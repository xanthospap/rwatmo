#include "geodesy/units.hpp"
#include "vmf3.hpp"

int dso::Vmf3SiteHandler::vmf3_impl(const char *site, const dso::MjdEpoch &t,
                                    double el,
                                    vmf3::Vmf3Result &result) noexcept {
  /* check we are at the correct time interval */
  if (!(t >= mt0 && t < mt1)) {
    fprintf(stderr,
            "[ERROR] Cannot compute vmf3! Wrong time interval, given t=%.9f "
            "while range is [%.9f to %.9f] (traceback: %s)\n",
            t.as_mjd(), mt0.as_mjd(), mt1.as_mjd(), __func__);
  }

  /* find the site of interest */
  const auto it = std::find_if(msites.begin(), msites.end(),
                               [&site](const vmf3::SiteBlock &s) {
                                 return !std::strcmp(s.msite.name(), site);
                               });
  if (it == msites.end()) {
    fprintf(stderr,
            "[ERROR] Site %s is not included in the Vmf3SiteHandler! "
            "(traceback: %s)\n",
            site, __func__);
    fprintf(stderr,
            "[ERROR] C'nued Cannot compute VMF3 mapping function (traceback: "
            "%s)\n",
            __func__);
    return 1;
  }

  const auto dt0 =
      t.diff<dso::DateTimeDifferenceType::FractionalSeconds>(mt0).seconds();
  [[maybe_unused]] const auto dt1 =
      mt1.diff<dso::DateTimeDifferenceType::FractionalSeconds>(t).seconds();
  const double dt10 =
      mt1.diff<dso::DateTimeDifferenceType::FractionalSeconds>(mt0).seconds();

  /* store results per cell (4 nodes) */
  vmf3::Vmf3Result res[4];
  dso::vmf3::GridVmf3Data::Data d[4]; /* i.e. [bl, br, tl, tr] */
  const double w = dt0 / dt10;        /* for linear interpolation */

  for (int node = 0; node < 4; node++) {
    /* step 1: Linear interpolation in time */
    const dso::vmf3::GridVmf3Data::Data *__restrict__ d0 =
        &(it->mdata_t0[node]);
    const dso::vmf3::GridVmf3Data::Data *__restrict__ d1 =
        &(it->mdata_t1[node]);
    for (int k = 0; k < dso::vmf3::GridVmf3Data::Data::NUM_DATA_ELEMENTS; k++) {
      d[node].data()[k] = d0->data()[k] + w * (d1->data()[k] - d0->data()[k]);
    }

    /* Step 2: Compute ZHD, ZWD, MFH, and MFW at each node */

    /* lift zhd/zwd to station height:
     * (1) convert zhd at grid height to the respective pressure value */
    const dso::vmf3::GridVmf3Data::Data *__restrict__ dn = &(d[node]);
    const double cos2lat = std::cos(2e0 * dso::deg2rad(dn->lat_deg()));
    const double pg = (dn->zhd() / 0.0022768) *
                      (1. - 0.00266 * cos2lat - 0.28e-6 * dn->oro_ell());

    /* (2) lift the pressure each from grid height to site height */
    const double ps =
        pg * std::pow(1. - 2.26e-5 * (it->mcrd.hgt() - dn->oro_ell()), 5.225);

    /* (3) convert the lifted pressure to zhd again (Kouba 2008) */
    res[node].zhd() =
        0.0022768 * ps / (1. - 0.00266 * cos2lat - 0.28e-6 * it->mcrd.hgt());

    /* zwd to site height */
    res[node].zwd() =
        dn->zwd() * std::exp(-(it->mcrd.hgt() - dn->oro_ell()) / 2e3);

    /* mapping functions */
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