#ifndef __DSO_VMF3_TROPO_MODELING_HPP__
#define __DSO_VMF3_TROPO_MODELING_HPP__

#include "vmf3.hpp"
#include "vmf3_grid_stream.hpp"

namespace dso {

/* Contains a Site and grid data from interpolation (for some epoch). */
struct SiteBlock {
  vmf3::Site msite;
  dso::GeodeticCrd mcrd;
  double mcos2lat;
  vmf3::GridVmf3Data::Data mdata_t0;
  vmf3::GridVmf3Data::Data mdata_t1;
  vmf3::Vmf3FullCoeffs msitebc;
}; /* struct SiteBlock */

class Vmf3SiteHandler {
 private:
  MjdEpoch mt0{MjdEpoch::max()};
  MjdEpoch mt1{MjdEpoch::min()};
  std::vector<vmf3::SiteBlock> msites;
  std::string mdata_dir;
  Vmf3 mvmf3;

  [[no_discard]]
  int dso::Vmf3SiteHandler::load_sites_for_epoch(const ymd_date &ymd,
                                                 int day_hours) noexcept;

  auto vmf3(const char *site, const MjdEpoch &t, double el) {
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
        t.diff<dso::DateTimeDifferenceType::FractionalSeconds>(mt0).seconds();
    const auto dt1 =
        mt1.diff<dso::DateTimeDifferenceType::FractionalSeconds>(t).seconds();

    if ((dt0 >= 0 && dt0 < ThreeHoursInSec) &&
        (dt1 >= 0 && dt1 < ThreeHoursInSec)) {
      /* everything in place, got the correct epochs! temporal interpolation */
      dso::vmf3::GridVmf3Data::Data d;
      const double dt10 =
          mt1.diff<dso::DateTimeDifferenceType::FractionalSeconds>(mt0)
              .seconds();
      const double w = dt0 / dt10;
      for (int k = 0; k < dso::vmf3::GridVmf3Data::Data::NUM_ELEMENTS; k++) {
        d.data()[k] = it->mdata_t0.data()[k] +
                      w * (it->mdata_t1.data()[k] - it->mdata_t0.data()[k]);
      }
      /* interpolated coeffs are now in d, for the current site */

    } /* interval ok */

    /* lift zhd/zwd to station height:
     * (1) convert the hydrostatic zenith delay at grid height to the respective
     * pressure value
     */
    const double pg = (d.zhd() / 0.0022768) *
                      (1. - 0.00266 * it->mcos2lat - 0.28e-6 * it->oro_ell);
    /* (2) lift the pressure each from grid height to site height */
    const double ps =
        pg * std::pow(1. - 2.26e-5 * (it->mcrd.hgt() - it->oro_ell), 5.225);
    /* (3) convert the lifted pressure to zhd again (as proposed by Kouba, 2008)
     */
    d.zhd() = 0.0022768 * ps /
              (1. - 0.00266 * it->mcos2lat - 0.28e-6 * it->mcrd.hgt());
    /* zwd to site height */
    d.zwd() *= std::exp(-(it->mcrd.hgt() - it->oro_ell) / 2e3);

    /* compute mapping functions */
    mvmf3.mf(t, el, d.ah(), d.aw(), it->msitebc, mfh, mfw);
  }
}; /*class Vmf3SiteHandler */
} /* namespace dso  */

#endif
