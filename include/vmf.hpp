#ifndef __DSO_VMF3_TROPO_MODELING_HPP__
#define __DSO_VMF3_TROPO_MODELING_HPP__

#include "vmf3.hpp"
#include "vmf3_grid_data.hpp"

namespace dso {

namespace vmf3 {
/* Contains a Site and grid data from interpolation (for some epoch). */
struct SiteBlock {
  /* site name +  domes */
  vmf3::Site msite;
  /* ellispoidal coordinates */
  dso::GeodeticCrd mcrd;
  /* t0 epoch, corresponding to mdata_t0 */
  MjdEpoch mt0;
  /* grid data for surrounding nodes in the order: [bl, br, tl, tr] */
  vmf3::GridVmf3Data::Data mdata_t0[4];
  /* t1 epoch, corresponding to mdata_t1 */
  MjdEpoch mt1;
  /* grid data for surrounding nodes in the order: [bl, br, tl, tr] */
  vmf3::GridVmf3Data::Data mdata_t1[4];
  /* coeffs for each surrounding node, [bl, br, tl, tr] */
  vmf3::Vmf3FullCoeffs msitebc[4];
}; /* struct SiteBlock */

struct Vmf3Result {
  double zhd, zwd, mfh, mfw;
  double zhd() const noexcept { return zhd; }
  double zwd() const noexcept { return zwd; }
  double mfh() const noexcept { return mfh; }
  double mfw() const noexcept { return mfw; }
  double &zhd() noexcept { return zhd; }
  double &zwd() noexcept { return zwd; }
  double &mfh() noexcept { return mfh; }
  double &mfw() noexcept { return mfw; }
};
} /* namespace vmf3 */

class Vmf3SiteHandler {
 private:
  std::vector<vmf3::SiteBlock> msites;
  std::string mdata_dir;
  Vmf3 mvmf3;

  /** @brief Load a VMF3GR data grid and populate the instance's  */
  [[no_discard]]
  int dso::Vmf3SiteHandler::load_sites_for_epoch(const ymd_date &ymd,
                                                 int day_hours) noexcept;

  auto vmf3(const char *site, const MjdEpoch &t, double el,
            vmf3::Vmf3Result &result) {
    /* find the site of interest */
    const auto it = std::find_if(msites.begin(), msites.end(),
                                 [&site](const vmf3::SiteBlock &s) {
                                   return !std::strcmp(s.msite.name(), site);
                                 });
    if (it == msites.end()) {
      fprintf(stderr,
              "[ERROR] Site %s is not included in the Vmf3SiteHandler! "
              "(traceback: %s)\n",
              __func__);
      fprintf(stderr,
              "[ERROR] C'nued Cannot compute VMF3 mapping function (traceback: "
              "%s)\n",
              __func__);
      return 1;
    }

    /* make sure we have the correct time interval, else try to collect it */
    constexpr const double ThreeHoursInSec = 3 * 60 * 60e0;
    const auto dt0 =
        t.diff<dso::DateTimeDifferenceType::FractionalSeconds>(it->mt0)
            .seconds();
    const auto dt1 =
        it->mt1.diff<dso::DateTimeDifferenceType::FractionalSeconds>(t)
            .seconds();

    // if ((dt0 >= 0 && dt0 < ThreeHoursInSec) &&
    //     (dt1 >= 0 && dt1 < ThreeHoursInSec)) {

    vmf3::Vmf3Result res[4];

    /* everything in place, got the correct epochs! temporal interpolation */
    dso::vmf3::GridVmf3Data::Data d[4];  // i.e. [bl, br, tl, tr]
    const double dt10 =
        it->mt1.diff<dso::DateTimeDifferenceType::FractionalSeconds>(it->mt0)
            .seconds();
    const double w = dt0 / dt10;
    for (int node = 0; node < 4; node++) {
      const dso::vmf3::GridVmf3Data::Data *__restrict__ d0 = it->mdata_t0[node];
      const dso::vmf3::GridVmf3Data::Data *__restrict__ d1 = it->mdata_t1[node];
      for (int k = 0; k < dso::vmf3::GridVmf3Data::Data::NUM_ELEMENTS; k++) {
        d[node].data()[k] = d0->data()[k] + w * (d1->data()[k] - d0->data()[k]);
      }
    }

    /* d[4] holds data values for the epoch of request (i.e. t) after temporal
     * linear interpolation between t0 and t1
     *
     * Next step: Compute ZHD, ZWD, MFH, and MFW at each node, after applying
     * corrections.
     *
     * lift zhd/zwd to station height:
     * (1) convert the hydrostatic zenith delay at grid height to the respective
     * pressure value
     */
    for (int node = 0; node < 4; node++) {
      const dso::vmf3::GridVmf3Data::Data *__restrict__ dn = d[node];
      const double cos2lat = std::cos(2e0 * dn->lat());

      const double pg = (dn->zhd() / 0.0022768) *
                        (1. - 0.00266 * cos2lat - 0.28e-6 * it->oro_ell);
      /* (2) lift the pressure each from grid height to site height */
      const double ps =
          pg * std::pow(1. - 2.26e-5 * (it->mcrd.hgt() - it->oro_ell), 5.225);
      /* (3) convert the lifted pressure to zhd again (as proposed by Kouba,
       * 2008)
       */
      res[node].zhd() =
          0.0022768 * ps / (1. - 0.00266 * cos2lat - 0.28e-6 * it->mcrd.hgt());
      /* zwd to site height */
      res[node].zwd() *= std::exp(-(it->mcrd.hgt() - it->oro_ell) / 2e3);

      double mfh, mfw;
      mvmf3.mf(t, el, dn->ah(), dn->aw(), it->msitebc[node], res[node].mfh(),
               res[node].mfw());

      /* Apply Niell (1996) height correction to the hydrostatic mapping */
      const double aht = 2.53e0 - 5;
      const double bht = 5.49e0 - 3;
      const double cht = 1.14e0 - 3;
      const double sel = std::sin(el);
      const double nmfh = (1 + (aht / (1 + bht / (1 + cht)))) /
                          (sel + (aht / (sel + bht / (sel + cht))));
      res[node].mfh() += (1e0 / sel - nmfh) * it->mcrd.hgt() / 1e3;
    }

    /* done for all nodes, now bilinear interpolation */
    result = res[0];
    return 0;
  }
}; /*class Vmf3SiteHandler */
} /* namespace dso  */

#endif
