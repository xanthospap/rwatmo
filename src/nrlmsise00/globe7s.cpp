#include "nrlmsise.hpp"

double dso::Nrlmsise00::glob7s(const double *p, dso::Nrlmsise00::DataTrigs &dt,
                                double doy, double sec, double glon,
                                double f107, double f107A,
                                const double *const apt, const double plg[4][9]
                                ) const noexcept {
  /*    VERSION OF GLOBE FOR LOWER ATMOSPHERE 10/26/99
   */
  constexpr const double pset = 2.0;

  /* initialize temperatures */
  double t[14];
  for (int j = 0; j < 14; j++)
    t[j] = 0.0;

  const double cd32 = std::cos(dr * (doy - p[31]));
  const double cd18 = std::cos(2.0 * dr * (doy - p[17]));
  const double cd14 = std::cos(dr * (doy - p[13]));
  const double cd39 = std::cos(2.0 * dr * (doy - p[38]));

  /* F10.7 */
  t[0] = p[21] * (f107A-150e0);

  /* time independent */
  t[1] = p[1] * plg[0][2] + p[2] * plg[0][4] + p[22] * plg[0][6] +
         p[26] * plg[0][1] + p[14] * plg[0][3] + p[59] * plg[0][5];

  /* SYMMETRICAL ANNUAL */
  t[2] = (p[18] + p[47] * plg[0][2] + p[29] * plg[0][4]) * cd32;

  /* SYMMETRICAL SEMIANNUAL */
  t[3] = (p[15] + p[16] * plg[0][2] + p[30] * plg[0][4]) * cd18;

  /* ASYMMETRICAL ANNUAL */
  t[4] = (p[9] * plg[0][1] + p[10] * plg[0][3] + p[20] * plg[0][5]) * cd14;

  /* ASYMMETRICAL SEMIANNUAL */
  t[5] = (p[37] * plg[0][1]) * cd39;

  /* DIURNAL */
  const double t71 = p[11] * plg[1][2] * cd14;
  const double t72 = p[12] * plg[1][2] * cd14;
  t[6] = ((p[3] * plg[1][1] + p[4] * plg[1][3] + t71) * dt.ctloc +
          (p[6] * plg[1][1] + p[7] * plg[1][3] + t72) * dt.stloc);

  /* SEMIDIURNAL */
  const double t81 = (p[23] * plg[2][3] + p[35] * plg[2][5]) * cd14;
  const double t82 = (p[33] * plg[2][3] + p[36] * plg[2][5]) * cd14;
  t[7] = ((p[5] * plg[2][2] + p[41] * plg[2][4] + t81) * dt.c2tloc +
          (p[8] * plg[2][2] + p[42] * plg[2][4] + t82) * dt.s2tloc);

  /* TERDIURNAL */
  t[13] = p[39] * plg[3][3] * dt.s3tloc + p[40] * plg[3][3] * dt.c3tloc;

  /* MAGNETIC ACTIVITY */
  t[8] = (p[50] * apt[0] + p[96] * plg[0][2] * apt[0]);

  /* LONGITUDINAL */
  if (!(dso::rad2deg(glon) <= -1000.0)) {
    t[10] = (1.0 +
             plg[0][1] * (p[80] * std::cos(dr * (doy - p[81])) +
                          p[85] * std::cos(2.0 * dr * (doy - p[86]))) +
             p[83] * std::cos(dr * (doy - p[84])) +
             p[87] * std::cos(2.0 * dr * (doy - p[88]))) *
            ((p[64] * plg[1][2] + p[65] * plg[1][4] + p[66] * plg[1][6] +
              p[74] * plg[1][1] + p[75] * plg[1][3] + p[76] * plg[1][5]) *
                 std::cos(glon) +
             (p[90] * plg[1][2] + p[91] * plg[1][4] + p[92] * plg[1][6] +
              p[77] * plg[1][1] + p[78] * plg[1][3] + p[79] * plg[1][5]) *
                 std::sin(glon));
  }

  double tt = 0;
  for (int i = 0; i < 14; i++)
    tt += t[i];
  return tt;
}
