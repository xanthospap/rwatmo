#include <cmath>

namespace {
/* returns gv */
inline double glatf(double lat, double& reff)
{
  const double dgtr = 1.74533e-2e0;
  const double c2 = std::cos(2.0 * dgtr * lat);
  const double gv = 980.616 * (1e0 - 0.0026373 * c2);
  reff = 2e0 * (gv) / (3.085462e-6 + 2.27e-9 * c2) * 1e-5;
  return gv;
}
constexpr const double ZMIX = 62.5;
} /* anonymous namespace */

void dso::Nrlmsise00::gtd7(const dso::MjdEpoch& t, const Eigen::Vector3d& ecef,
    double f107A, double f107, const double* const ap, double* densities, double* temperatures)
{
  constexpr const double zn3[5] = { 32.5, 20.0, 15.0, 10.0, 0.0 };
  constexpr const double zn2[4] = { 72.5, 55.0, 45.0, 32.5 };
  constexpr const int mn3 = sizeof(zn3) / sizeof double;
  constexpr const int mn2 = sizeof(zn2) / sizeof double;
#ifdef DEBUG
  static_assert(mn3 == 5);
  static_assert(mn2 == 4);
#endif

  /* TODO ecef, cartesian to geographic/geodetic */
  double lat, lon, hgt;
  /* TODO compute latitude */
  double alt;
  assert(alt > zn2[0]);

  double re; /* radius */
  const double gsurf = glatf(xlat, re);

  /* thermosphere / mesosphere (above zn2[0]) */
  /* TODO gts7 computes and needs to return:
   * dm28
   * meso_tgn1[2]
   * meso_tn1[5]
   *
   *
   * gts7 computes temperatures and densities
   */
  this->gts7(t, alt, f107A, f107, ap, densities, temperatures);

  /*       LOWER MESOSPHERE/UPPER STRATOSPHERE (between zn3[0] and zn2[0])
   *         Temperature at nodes and gradients at end nodes
   *         Inverse temperature a linear function of spherical harmonics
   */
  meso_tgn2[0] = meso_tgn1[1];
  meso_tn2[0] = meso_tn1[4];
  meso_tn2[1] = pma[0][0] * pavgm[0] / (1e0 - glob7s(pma[0], input, flags));
  meso_tn2[2] = pma[1][0] * pavgm[1] / (1e0 - glob7s(pma[1], input, flags));
  meso_tn2[3] = pma[2][0] * pavgm[2] / (1e0 - glob7s(pma[2], input, flags));
  meso_tgn2[1] = pavgm[8] * pma[9][0] * (1e0 + glob7s(pma[9], input, flags)) * meso_tn2[3] * meso_tn2[3] / (pow((pma[2][0] * pavgm[2]), 2.0));
  meso_tn3[0] = meso_tn2[3];

  if (alt < zn3[0]) {
    /*       LOWER STRATOSPHERE AND TROPOSPHERE (below zn3[0])
     *         Temperature at nodes and gradients at end nodes
     *         Inverse temperature a linear function of spherical harmonics
     */
    meso_tgn3[0] = meso_tgn2[1];
    meso_tn3[1] = pma[3][0] * pavgm[3] / (1.0 - glob7s(pma[3], input, flags));
    meso_tn3[2] = pma[4][0] * pavgm[4] / (1.0 - glob7s(pma[4], input, flags));
    meso_tn3[3] = pma[5][0] * pavgm[5] / (1.0 - glob7s(pma[5], input, flags));
    meso_tn3[4] = pma[6][0] * pavgm[6] / (1.0 - glob7s(pma[6], input, flags));
    meso_tgn3[1] = pma[7][0] * pavgm[7] * (1.0 + glob7s(pma[7], input, flags)) * meso_tn3[4] * meso_tn3[4] / (pow((pma[6][0] * pavgm[6]), 2.0));
  } /* alt<zn3[0] */

  /* LINEAR TRANSITION TO FULL MIXING BELOW zn2[0] */
  dmc = 1e0 - (zn2[0] - alt) / (zn2[0] - ZMIX);
  dz28 = densities[2];

  /**** N2 density ****/
  dmr = dz28 / dm28m - 1.0;
  densities[2] = densm(alt, dm28m, xmm, &tz, mn3, zn3, meso_tn3, meso_tgn3, mn2, zn2, meso_tn2, meso_tgn2) * (1.0 + dmr * dmc);

  /**** HE density TODO ****/
  dmr = soutput.d[0] / (dz28 * pdm[0][1]) - 1.0;
  output->d[0] = output->d[2] * pdm[0][1] * (1.0 + dmr*dmc);
}
