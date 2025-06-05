#include "nrlmsise.hpp"

double dso::Nrlmsise00::densu(double alt, double dlb, double tinf, double tlb, double xm,
    double alpha, double &tz, double zlb, double s2, int mn1,
    const double *zn1, double *tn1, double *tgn1, double *u) {

  /* Calculate Temperature and Density Profiles for MSIS models
   * New lower thermo polynomial
   */
  constexpr const double rgas = 831.4;
  double xs[5], ys[5], y2out[5];

  /* joining altitudes of Bates and spline */
  const double z = (alt > zn1[0]) ? alt : zn1[0];

  /* geopotential altitude difference from ZLB */
  const double zg2 = zeta_(z, zlb);

  const double tt = tinf - (tinf - tlb) * std::exp(-s2 * zg2);

  /* Bates temperature */
  tz = tt;
  if (alt < zn1[0]) {
    /* calculate temperature below ZA
     * temperature gradient at ZA from Bates profile
     */
    const double dta =
      (tinf - ta) * s2 * std::pow(((re + zlb) / (re + za)), 2.0);
    tgn1[0] = dta;
    tn1[0] = tt;
    const double z = (alt > zn1[mn1 - 1]) ? alt : zn1[mn1 - 1];
    /* geopotental difference from z1 */
    zg = zeta_(z, zn1[0]);
    zgdif = zeta_(zn1[mn1 - 1], zn1[0]);
    /* set up spline nodes */
    for (k = 0; k < mn1; k++) {
      xs[k] = zeta_(zn1[k], zn1[0]) / zgdif;
      ys[k] = 1.0 / tn1[k];
    }
    /* end node derivatives */
    const double yd1 = -tgn1[0] / (tn1[0] * tn1[0]) * zgdif;
    const double yd2 =
      -tgn1[1] / (tn1[mn1 - 1] * tn1[mn1 - 1]) * zgdif * std::pow(((re + zn1[mn1 - 1]) / (re + zn1[0])), 2.0);
    /* calculate spline coefficients */
    this->spline(xs, ys, mn1, yd1, yd2, y2out, u);
    double x = zg / zgdif;
    const double y = this->splint(xs, ys, y2out, mn1, x);
    /* temperature at altitude */
    tz = 1.0 / y;
  }
  if (xm == 0) [[unlikely]]
    return tz;

  /* calculate density above za */
  double glb = gsurf / std::pow((1.0 + zlb / re), 2.0);
  double gamma = xm * glb / (s2 * rgas * tinf);
  double expl = std::exp(-s2 * gamma * zg2);
  if (expl > 50.0)
    expl = 50.0;
  if (tt <= 0)
    expl = 50.0;

  /* density at altitude */
  const double densa =
    dlb * std::pow((tlb / tt), ((1.0 + alpha + gamma))) * expl;
  if (alt >= zn1[0])
    return densa;

  /* calculate density below za */
  glb = gsurf / std::pow((1.0 + z1 / re), 2.0);
  gamm = xm * glb * zgdif / rgas;

  /* integrate spline temperatures */
  const double yi = this->splini(xs, ys, y2out, mn1, x);
  expl = gamm * yi;
  if (expl > 50.0)
    expl = 50.0;
  if (tz <= 0)
    expl = 50.0;

  /* density at altitude */
  return densa * std::pow((t1 / tz), (1.0 + alpha)) * std::exp(-expl);
}
