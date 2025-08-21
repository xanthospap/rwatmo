#include "vmf3.hpp"

int dso::Vmf3SiteHandler::initialize(const dso::MjdEpoch &t) noexcept {
  /* do we also need to initialize the Vmf3 b,c coeffs ? */
  bool initialize_coeffs =
      (mt0 == dso::MjdEpoch::max() && mt1 == dso::MjdEpoch::min()) ? 1 : 0;

  /* find surounding epochs for t */
  const double hours = t.seconds().seconds() / 3600e0;
  const int intrv = (int)(hours / (double)FILE_TIME_INTERVAL);
  const int hp = intrv * FILE_TIME_INTERVAL;
  const int hn =
      ((hp + FILE_TIME_INTERVAL) >= 24) ? 0 : (hp + FILE_TIME_INTERVAL);
  if (hp < 0 || hn > 18) {
    fprintf(stderr,
            "[ERROR] Failed locating a time interval for given date "
            "(traceback: %s)\n",
            __func__);
    return 1;
  }

  /* load the file prior to t. if this is the first loading/initialization, we
   * should also probably load the orography file
   */
  bool load_orography = (initialize_coeffs) ? true : false;
  if (load_sites_for_epoch(t.to_ymd(), hp, load_orography)) {
    fprintf(stderr,
            "[ERROR] Failed loading VMF3(GR) grid data file (traceback: %s)\n",
            __func__);
    return 1;
  }
  /* good, but t0 is now loaded as t1; left shift the data */
  this->left_shift();

  /* load the file next to t */
  auto ymdn = t.to_ymd();
  if (hn == 0) ymdn = t.add_seconds(FractionalSeconds(3599e0)).to_ymd();
  if (load_sites_for_epoch(ymdn, hn)) {
    fprintf(stderr,
            "[ERROR] Failed loading VMF3(GR) grid data file (traceback: %s)\n",
            __func__);
    return 1;
  }
  /* good, t1 now holds the data immidiate after t */
  this->left_shift();

  /* if needed, compute empirical vmf3 coeffs for all collected cells */
  if (initialize_coeffs) {
    if (compute_spatial_vmf3_coeffs()) {
      fprintf(stderr,
              "[ERROR] Failed computing empirical (b,c) coefficients for "
              "VMF3 (traceback: %s)\n",
              __func__);
      return 1;
    }
  }

  return (!grid_matches_between_epochs());
}

int dso::Vmf3SiteHandler::load_correct_interval(
    const dso::MjdEpoch &t) noexcept {
  constexpr const double INTRV = FILE_TIME_INTERVAL * 3600e0;
  if (t >= mt0 && t < mt1) return 0;

  /* the most probable case would be that we overrun the interval [t0, t1]
   * and we now need [t1, t1+1]
   */
  if ((mt1 != MjdEpoch::min()) &&
      (mt1 <= t &&
       t.diff<dso::DateTimeDifferenceType::FractionalSeconds>(mt1).seconds() <
           INTRV)) {
    /* first move t1 to t0 */
    this->left_shift();
    /* find date and hours of day for next file */
    auto tn = mt0.add_seconds(FractionalSeconds(INTRV));
    const double hours = tn.seconds().seconds() / 3600e0;
    const int intrv = (int)(hours / (double)FILE_TIME_INTERVAL);
    const int hn = intrv * FILE_TIME_INTERVAL;
    /* load the file prior to t */
    if (load_sites_for_epoch(tn.to_ymd(), hn)) {
      fprintf(
          stderr,
          "[ERROR] Failed loading VMF3(GR) grid data file (traceback: %s)\n",
          __func__);
      return 1;
    }
  }

  /* if this is not the case fuck it, re-initialize */
  if (initialize(t)) {
    fprintf(
        stderr,
        "[ERROR] Failed loading VMF3(GR) grid data file(s) (traceback: %s)\n",
        __func__);
    return 1;
  }

  return (!grid_matches_between_epochs());
}