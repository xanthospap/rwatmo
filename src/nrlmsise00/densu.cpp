#include "nrlmsise.hpp"
#include "rwatmo_math.hpp"

double dso::Nrlmsise00::densu(double lat, double alt, double dlb, double tinf,
                              double tlb, double xm, double alpha, double &tz,
                              double zlb, double s2, int mn1, const double *zn1,
                              double *tn1, double *tgn1,
                              double *u) const noexcept {

  /* Calculate Temperature and Density Profiles for MSIS models
   * New lower thermo polynomial
   */
  double xs[5], ys[5], y2out[5];

  /* joining altitudes of Bates and spline */
  double z = (alt > zn1[0]) ? alt : zn1[0];
  
  double re;
  const double gsurf = glatf(lat, re);

  /* geopotential altitude difference from ZLB */
  const double zg2 = ((z - zlb) * (re + zlb) / (re + z));

  const double tt = tinf - (tinf - tlb) * std::exp(-s2 * zg2);

  /* Bates temperature */
  tz = tt;
  double z1 = 0e0;
  double zgdif = 0e0;
  double x = 0e0;
  int mn = 0;
  double t1 = 0e0;

  if (alt < zn1[0]) {
    /* calculate temperature below ZA temperature gradient at ZA from Bates
     * profile */
    z1 = zn1[0];
    mn = mn1;
    t1 = tn1[0];

    const double dta =
        (tinf - tt) * s2 * std::pow(((re + zlb) / (re + zn1[0])), 2.0);
    tgn1[0] = dta;
    tn1[0] = tt;
    z = (alt > zn1[mn1 - 1]) ? alt : zn1[mn1 - 1];
    /* geopotental difference from z1 */
    const double zg = ((z - zn1[0]) * (re + zn1[0]) / (re + z));
    zgdif = ((zn1[mn1-1] - zn1[0]) * (re + zn1[0]) / (re + zn1[mn1-1]));
    /* set up spline nodes */
    for (int i = 0; i < mn1; i++) {
      xs[i] = ((zn1[i] - zn1[0]) * (re + zn1[0]) / (re + zn1[i])) / zgdif;
      ys[i] = 1.0 / tn1[i];
    }
    /* end node derivatives */
    const double yd1 = -tgn1[0] / (tn1[0] * tn1[0]) * zgdif;
    const double yd2 = -tgn1[1] / (tn1[mn1 - 1] * tn1[mn1 - 1]) * zgdif *
                       std::pow(((re + zn1[mn1 - 1]) / (re + zn1[0])), 2.0);
    /* calculate spline coefficients */
    spline(xs, ys, mn1, yd1, yd2, y2out, u);
    x = zg / zgdif;
    const double y = splint(xs, ys, y2out, mn1, x);
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
  gamma = xm * glb * zgdif / rgas;

  /* integrate spline temperatures */
  const double yi = splini(xs, ys, y2out, mn, x);
  expl = gamma * yi;
  if (expl > 50.0)
    expl = 50.0;
  if (tz <= 0)
    expl = 50.0;

  /* density at altitude */
  return densa * std::pow((t1 / tz), (1.0 + alpha)) * std::exp(-expl);
}
