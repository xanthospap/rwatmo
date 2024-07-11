#include "vmf3.hpp"

dso::detail::Vmf3EmpiricalCoeffs dso::detail::Vmf3FullCoeffs::computeCoeffs(
    const dso::MjdEpoch &t) const noexcept {
  /* fractional day of year */
  const dso::ydoy_date yd(t.to_ydoy());
  const double doy = yd.dy().as_underlying_type() + t.fractional_days();

  constexpr const double pi = M_PI;
  dso::Vmf3EmpiricalCoeffs vc;

  /* trigs */
  const double sy2 = std::sin(doy / 365.25 * 2 * pi);
  const double cy2 = std::cos(doy / 365.25 * 2 * pi);
  const double sy4 = std::sin(doy / 365.25 * 4 * pi);
  const double cy4 = std::cos(doy / 365.25 * 4 * pi);

  /* add seasonal amplitudes for the specified doy to the mean values */
  vc.bh = bh_A0 + bh_A1 * cy2 + bh_B1 * sy2 + bh_A2 * cy4 + bh_B2 * sy4;
  vc.bw = bw_A0 + bw_A1 * cy2 + bw_B1 * sy2 + bw_A2 * cy4 + bw_B2 * sy4;
  vc.ch = ch_A0 + ch_A1 * cy2 + ch_B1 * sy2 + ch_A2 * cy4 + ch_B2 * sy4;
  vc.cw = cw_A0 + cw_A1 * cy2 + cw_B1 * sy2 + cw_A2 * cy4 + cw_B2 * sy4;

  return vc;
}
