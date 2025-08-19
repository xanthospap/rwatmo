#ifndef __DSO_VMF3_TROPO_MODELING_HPP__
#define __DSO_VMF3_TROPO_MODELING_HPP__

#include "vmf3_grid_stream.hpp"

namespace dso {
class Vmf3 {}; /* class Vmf3 */

class Vmf3SiteHandler {
 private:
  MjdEpoch mt0{MjdEpoch::max()};
  MjdEpoch mt1{MjdEpoch::min()};
  std::vector<vmf3::SiteBlock> msites;
  std::string mdata_dir;

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
      /* interpolated coeffs are not in d, for the current site */

    } /* interval ok */
  }
}; /*class Vmf3SiteHandler */
} /* namespace dso  */

#endif
