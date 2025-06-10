#include "nrlmsise.hpp"

namespace {
/*    3hr Magnetic activity functions */
/*    Eq. A24d */
double g0(double a, const double *const p) {
  const double p24 = (p[24] < 1.0e-4) ? 1.0e-4 : p[24];
  return (a - 4.0 +
          (p[25] - 1.0) * (a - 4.0 +
                           (std::exp(-std::sqrt(p24 * p24) * (a - 4.0)) - 1.0) /
                               std::sqrt(p24 * p24)));
}
/*    Eq. A24c */
double sumex(double ex) {
  return (1.0 + (1.0 - std::pow(ex, 19.0)) / (1.0 - ex) * std::pow(ex, 0.5));
}

/*    Eq. A24a */
double sg0(double ex, const double *const p, const double *const ap) {
  return (g0(ap[1], p) + (g0(ap[2], p) * ex + g0(ap[3], p) * ex * ex +
                          g0(ap[4], p) * std::pow(ex, 3.0) +
                          (g0(ap[5], p) * std::pow(ex, 4.0) +
                           g0(ap[6], p) * std::pow(ex, 12.0)) *
                              (1.0 - std::pow(ex, 8.0)) / (1.0 - ex))) /
         sumex(ex);
}
} // namespace

double dso::Nrlmsise00::glob7(const double *const p,
                               const dso::Nrlmsise00::DataTrigs &dt,
                               /*double tloc,*/ double doy, double sec,
                               /*double glat,*/ double glon, double f107,
                               double f107A, const double *const ap,
                               double plg[4][9], double *apt) const noexcept {
  /* CALCULATE G(L) FUNCTION Upper Thermosphere Parameters */
  double t[15] = {0e0};

  /* calculate legendre polynomials */
  const double c = dt.c; // std::sin(glat);
  const double s = dt.s; // std::cos(glat);
  const double c2 = c * c;
  const double c4 = c2 * c2;
  const double s2 = s * s;

  if (plg[0][1] == 0e0) {
    plg[0][1] = c;
    plg[0][2] = 0.5 * (3.0 * c2 - 1.0);
    plg[0][3] = 0.5 * (5.0 * c * c2 - 3.0 * c);
    plg[0][4] = (35.0 * c4 - 30.0 * c2 + 3.0) / 8.0;
    plg[0][5] = (63.0 * c2 * c2 * c - 70.0 * c2 * c + 15.0 * c) / 8.0;
    plg[0][6] = (11.0 * c * plg[0][5] - 5.0 * plg[0][4]) / 6.0;
    /* plg[0][7] = (13.0*c*plg[0][6] - 6.0*plg[0][5])/7.0; */
    plg[1][1] = s;
    plg[1][2] = 3.0 * c * s;
    plg[1][3] = 1.5 * (5.0 * c2 - 1.0) * s;
    plg[1][4] = 2.5 * (7.0 * c2 * c - 3.0 * c) * s;
    plg[1][5] = 1.875 * (21.0 * c4 - 14.0 * c2 + 1.0) * s;
    plg[1][6] = (11.0 * c * plg[1][5] - 6.0 * plg[1][4]) / 5.0;
    /* plg[1][7] = (13.0*c*plg[1][6]-7.0*plg[1][5])/6.0; */
    /* plg[1][8] = (15.0*c*plg[1][7]-8.0*plg[1][6])/7.0; */
    plg[2][2] = 3.0 * s2;
    plg[2][3] = 15.0 * s2 * c;
    plg[2][4] = 7.5 * (7.0 * c2 - 1.0) * s2;
    plg[2][5] = 3.0 * c * plg[2][4] - 2.0 * plg[2][3];
    plg[2][6] = (11.0 * c * plg[2][5] - 7.0 * plg[2][4]) / 4.0;
    plg[2][7] = (13.0 * c * plg[2][6] - 8.0 * plg[2][5]) / 5.0;
    plg[3][3] = 15.0 * s2 * s;
    plg[3][4] = 105.0 * s2 * s * c;
    plg[3][5] = (9.0 * c * plg[3][4] - 7. * plg[3][3]) / 2.0;
    plg[3][6] = (11.0 * c * plg[3][5] - 8. * plg[3][4]) / 3.0;
  }

  const double stloc = dt.stloc;
  const double ctloc = dt.ctloc;
  const double s2tloc = dt.s2tloc;
  const double c2tloc = dt.c2tloc;
  const double s3tloc = dt.s3tloc;
  const double c3tloc = dt.c3tloc;

  const double cd32 = std::cos(dr * (doy - p[31]));
  const double cd18 = std::cos(2.0 * dr * (doy - p[17]));
  const double cd14 = std::cos(dr * (doy - p[13]));
  const double cd39 = std::cos(2.0 * dr * (doy - p[38]));

  /* init all temperatures to 0 */
  for (int i = 0; i < 14; i++)
    t[i] = 0e0;

  /* F10.7 EFFECT */
  const double df = f107 - f107A;
  const double dfa = f107A - 150.0;
  t[0] = p[19] * df * (1.0 + p[59] * dfa) + p[20] * df * df + p[21] * dfa +
         p[29] * std::pow(dfa, 2.0);
  const double f1 = 1.0 + (p[47] * dfa + p[19] * df + p[20] * df * df);
  const double f2 = 1.0 + (p[49] * dfa + p[19] * df + p[20] * df * df);

  /*  TIME INDEPENDENT */
  t[1] = (p[1] * plg[0][2] + p[2] * plg[0][4] + p[22] * plg[0][6]) +
         (p[14] * plg[0][2]) * dfa + p[26] * plg[0][1];

  /*  SYMMETRICAL ANNUAL */
  t[2] = p[18] * cd32;

  /*  SYMMETRICAL SEMIANNUAL */
  t[3] = (p[15] + p[16] * plg[0][2]) * cd18;

  /*  ASYMMETRICAL ANNUAL */
  t[4] = f1 * (p[9] * plg[0][1] + p[10] * plg[0][3]) * cd14;

  /*  ASYMMETRICAL SEMIANNUAL */
  t[5] = p[37] * plg[0][1] * cd39;

  /* DIURNAL */
  const double t71 = (p[11] * plg[1][2]) * cd14;
  const double t72 = (p[12] * plg[1][2]) * cd14;
  t[6] =
      f2 *
      ((p[3] * plg[1][1] + p[4] * plg[1][3] + p[27] * plg[1][5] + t71) * ctloc +
       (p[6] * plg[1][1] + p[7] * plg[1][3] + p[28] * plg[1][5] + t72) * stloc);

  /* SEMIDIURNAL */
  const double t81 = (p[23] * plg[2][3] + p[35] * plg[2][5]) * cd14;
  const double t82 = (p[33] * plg[2][3] + p[36] * plg[2][5]) * cd14;
  t[7] = f2 * ((p[5] * plg[2][2] + p[41] * plg[2][4] + t81) * c2tloc +
               (p[8] * plg[2][2] + p[42] * plg[2][4] + t82) * s2tloc);

  /* TERDIURNAL */
  t[13] =
      f2 *
      ((p[39] * plg[3][3] + (p[93] * plg[3][4] + p[46] * plg[3][6]) * cd14) *
           s3tloc +
       (p[40] * plg[3][3] + (p[94] * plg[3][4] + p[48] * plg[3][6]) * cd14) *
           c3tloc);

  /* magnetic activity based on daily ap */
  if (p[51] != 0) {
    double exp1 =
        std::exp(-10800.0 * std::sqrt(p[51] * p[51]) /
                 (1.0 + p[138] * (45.0 - std::abs(dso::rad2deg(dt.glat)))));
    if (exp1 > 0.99999)
      exp1 = 0.99999;
    apt[0] = sg0(exp1, p, ap);
    /* apt[1]=sg2(exp1,p,ap->a);
       apt[2]=sg0(exp2,p,ap->a);
       apt[3]=sg2(exp2,p,ap->a);
       */
    t[8] =
        apt[0] *
        (p[50] + p[96] * plg[0][2] + p[54] * plg[0][4] +
         (p[125] * plg[0][1] + p[126] * plg[0][3] + p[127] * plg[0][5]) * cd14 +
         (p[128] * plg[1][1] + p[129] * plg[1][3] + p[130] * plg[1][5]) *
             std::cos(dso::hours2rad(dt.tloc - p[131])));
  }

  if (dso::rad2deg(glon) > -1000.0) {
    /* longitudinal */
    t[10] = (1.0 + p[80] * dfa) *
            ((p[64] * plg[1][2] + p[65] * plg[1][4] + p[66] * plg[1][6] +
              p[103] * plg[1][1] + p[104] * plg[1][3] + p[105] * plg[1][5] +

              (p[109] * plg[1][1] + p[110] * plg[1][3] + p[111] * plg[1][5]) *
                  cd14) *
                 std::cos(glon) +
             (p[90] * plg[1][2] + p[91] * plg[1][4] + p[92] * plg[1][6] +
              p[106] * plg[1][1] + p[107] * plg[1][3] + p[108] * plg[1][5] +

              (p[112] * plg[1][1] + p[113] * plg[1][3] + p[114] * plg[1][5]) *
                  cd14) *
                 std::sin(glon));

    /* ut and mixed ut, longitude */
    t[11] = (1.0 + p[95] * plg[0][1]) * (1.0 + p[81] * dfa) *
            (1.0 + p[119] * plg[0][1] * cd14) *
            ((p[68] * plg[0][1] + p[69] * plg[0][3] + p[70] * plg[0][5]) *
             std::cos(dso::hsec2rad(sec - p[71])));
    t[11] += (p[76] * plg[2][3] + p[77] * plg[2][5] + p[78] * plg[2][7]) *
             std::cos(dso::hsec2rad(sec - p[79]) + 2.0 * glon) *
             (1.0 + p[137] * dfa);

    /* ut, longitude magnetic activity */
    if (p[51]) {
      t[12] =
          apt[0] * (1. + p[132] * plg[0][1]) *
              ((p[52] * plg[1][2] + p[98] * plg[1][4] + p[67] * plg[1][6]) *
               std::cos(glon - p[97])) +
          apt[0] *
              (p[133] * plg[1][1] + p[134] * plg[1][3] + p[135] * plg[1][5]) *
              cd14 * std::cos(glon - p[136]) +
          apt[0] * (p[55] * plg[0][1] + p[56] * plg[0][3] + p[57] * plg[0][5]) *
              std::cos(dso::hsec2rad(sec - p[58]));
    }
  }

  /* parms not used: 82, 89, 99, 139-149 */
  double tinf = p[30];
  for (int i = 0; i < 14; i++)
    tinf += t[i];
  return tinf;
}
