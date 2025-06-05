#include "geodesy/transformations.hpp"
#include "nrlmsise.hpp"
#include <cmath>
#include <stdexcept>

/* lat = geodetic latitude
 *
 * Computes:
 * densities[0-8]
 * temperatures[0]
 *
 * */
double dso::Nrlmsise00::gts7(const dso::MjdEpoch &t,
                             const dso::GeodeticCrd &flh, double f107A,
                             double f107, const double *const ap,
                             double *densities, double *temperatures) {

  const double sec = t.seconds.seconds();
  dso::ydoy_date yd = t.to_ydoy();
  int doy = yd.day_of_year().as_underlying_type();
  const double tloc = sec/3600e0 + dso::rad2deg(flh.lon())/15.e0;

  /* altitude in km */
  const double altitude =
      (flh.h() - dso::mean_earth_radius<dso::Ellipsoid::grs80>()) * 1e-3;
  if (altitude<72.5) {
    fprintf(stderr, "[ERROR] Cannot compute density for altitude < 72.5 km. Need to append to implementation! (traceback: %s)\n", __func__);
    std::string ermsg = "Cannot compute density for altitude < 72.5 km. Need to append to implementation!\n";
    throw std::runtime_error(ermsg);
  }

  double zn1[5] = {120.0, 110.0, 100.0, 90.0, 72.5};
  const double za = pdl[1][15];
  zn1[0] = za;

  /* initialize all densities to 0 */
  for (int j = 0; j < 9; j++)
    densities[j] = 0e0;

  /* legendre polynomials */
  double plg[4][9];
  plg[0][1] = 0e0; /* signal computation of polynomials */

  /* tinf variations not important below za or zn1(1) */
  temperatures[0] =
      ptm[0] * pt[0] +
      (altitude > zn1[0]) *
          (1.0 + dso::nrlmsise_impl::globe7(pt, tloc, doy, sec, flh.lat(),
                                            flh.lon(), f107, f107A, ap, plg));

  /*  gradient variations not important below zn1(5) */
  double g0;
  if (altitude > zn1[4]) {
    g0 = ptm[3] * ps[0] *
         (1.0 + dso::nrlmsise_impl::globe7(ps, tloc, doy, sec, flh.lat(),
                                           flh.lon(), f107, f107A, ap, plg));
  } else {
    g0 = ptm[3] * ps[0];
  }
  const double tlb =
      ptm[1] *
      (1.0 + dso::nrlmsise_impl::globe7(pd[3], tloc, doy, sec, flh.lat(),
                                        flh.lon(), f107, f107A, ap, plg)) *
      pd[3][0];
  const double s = g0 / (tinf - tlb);

  double meso_tn1[5];
  double meso_tn2[4];
  double meso_tn3[5];
  double meso_tgn1[2];
  double meso_tgn2[2];
  double meso_tgn3[2];

  /* Lower thermosphere temp variations not significant for density above 300 km
   */
  if (altitude < 300.0) {
    meso_tn1[1] = ptm[6] * ptl[0][0] /
                  (1.0 - glob7s(ptl[0], tloc, doy, sec, flh.lat(), flh.lon(),
                                f107, f107A, ap, plg));
    meso_tn1[2] = ptm[2] * ptl[1][0] /
                  (1.0 - glob7s(ptl[1], tloc, doy, sec, flh.lat(), flh.lon(),
                                f107, f107A, ap, plg));
    meso_tn1[3] = ptm[7] * ptl[2][0] /
                  (1.0 - glob7s(ptl[2], tloc, doy, sec, flh.lat(), flh.lon(),
                                f107, f107A, ap, plg));
    meso_tn1[4] = ptm[4] * ptl[3][0] /
                  (1.0 - glob7s(ptl[3], tloc, doy, sec, flh.lat(), flh.lon(),
                                f107, f107A, ap, plg));
    meso_tgn1[1] = ptm[8] * pma[8][0] *
                   (1.0 + glob7s(pma[8], tloc, doy, sec, flh.lat(), flh.lon(),
                                 f107, f107A, ap, plg)) *
                   meso_tn1[4] * meso_tn1[4] /
                   (std::pow((ptm[4] * ptl[3][0]), 2.0));
  } else {
    meso_tn1[1] = ptm[6] * ptl[0][0];
    meso_tn1[2] = ptm[2] * ptl[1][0];
    meso_tn1[3] = ptm[7] * ptl[2][0];
    meso_tn1[4] = ptm[4] * ptl[3][0];
    meso_tgn1[1] = ptm[8] * pma[8][0] * meso_tn1[4] * meso_tn1[4] /
                   (std::pow((ptm[4] * ptl[3][0]), 2.0));
  }

  /* N2 variation factor at Zlb */
  const double g28 = dso::nrlmsise_impl::globe7(pd[2], tloc, doy, sec, flh.lat(), flh.lon(), f107, f107A, ap, plg);

  /* VARIATION OF TURBOPAUSE HEIGHT */
  zhf = pdl[1][24] * (1.0 + pdl[0][24] * std::sin(flh.lat()) *
                                std::cos(dr * (doy - pt[13])));
  const double xmm = pdm[2][4];
  const double z = altitude;
  double tz;

  /**** N2 DENSITY ****/
  {
    /* Diffusive density at Zlb */
    const double db28 = pdm[2][0] * std::exp(g28) * pd[2][0];
    /* Diffusive density at Alt */
    densities[2] = densu(z, db28, tinf, tlb, 28.0, alpha[2], temperatures[1],
                         ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
    if (z <= altl[2]) {
      /* Turbopause */
      const double zh28 = pdm[2][2] * zhf;
      const double zhm28 = pdm[2][3] * pdl[1][5];
      const double xmd = 28.0 - xmm;
      /* Mixed density at Zlb */
      const double b28 = densu(zh28, db28, tinf, tlb, xmd, (alpha[2] - 1.0), tz,
                               ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
      /*  Mixed density at Alt */
      const double dm28 = densu(z, b28, tinf, tlb, xmm, alpha[2], tz, ptm[5],
                                s, mn1, zn1, meso_tn1, meso_tgn1);
      /*  Net density at Alt */
      densities[2] = dnet(densities[2], dm28, zhm28, xmm, 28.0);
    }
  }

  /**** HE DENSITY ****/
  {
    /*   Density variation factor at Zlb */
    const double g4 = globe7(pd[0], tloc, doy, sec, flh.lat(), flh.lon(), f107, f107A, ap, plg);
    /*  Diffusive density at Zlb */
    const double db04 = pdm[0][0] * std::exp(g4) * pd[0][0];
    /*  Diffusive density at Alt */
    densities[0] = densu(z, db04, tinf, tlb, 4., alpha[0], temperatures[1],
                         ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
    if (z < altl[0]) {
      /*  Turbopause */
      const double zh04 = pdm[0][2];
      /*  Mixed density at Zlb */
      const double b04 =
          densu(zh04, db04, tinf, tlb, 4. - xmm, alpha[0] - 1., temperatures[1],
                ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
      /*  Mixed density at Alt */
      const double dm04 = densu(z, b04, tinf, tlb, xmm, 0., temperatures[1],
                                ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
      const double zhm04 = zhm28;
      /*  Net density at Alt */
      densities[0] = dnet(densities[0], dm04, zhm04, xmm, 4.);
      /*  Correction to specified mixing ratio at ground */
      const double rl = log(b28 * pdm[0][1] / b04);
      const double zc04 = pdm[0][4] * pdl[1][0];
      const double hc04 = pdm[0][5] * pdl[1][1];
      /*  Net density corrected at Alt */
      densities[0] *= ccor(z, rl, hc04, zc04);
    }
  }

  /**** O DENSITY ****/
  {
    /*  Density variation factor at Zlb */
    const double g16 = globe7(pd[1], tloc, doy, sec, flh.lat(), flh.lon(), f107, f107A, ap, plg);
    /*  Diffusive density at Zlb */
    const double db16 = pdm[1][0] * std::exp(g16) * pd[1][0];
    /*   Diffusive density at Alt */
    densities[1] = densu(z, db16, tinf, tlb, 16., alpha[1], temperatures[1],
                         ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
    if (z <= altl[1]) {
      /*   Turbopause */
      const double zh16 = pdm[1][2];
      /*  Mixed density at Zlb */
      const double b16 =
          densu(zh16, db16, tinf, tlb, 16.0 - xmm, (alpha[1] - 1.0),
                temperatures[1], ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
      /*  Mixed density at Alt */
      const double dm16 = densu(z, b16, tinf, tlb, xmm, 0., temperatures[1],
                                ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
      const double zhm16 = zhm28;
      /*  Net density at Alt */
      densities[1] = dnet(densities[1], dm16, zhm16, xmm, 16.);
      const double rl =
          pdm[1][1] * pdl[1][16] * (1.0 + pdl[0][23] * (f107A - 150.0));
      const double hc16 = pdm[1][5] * pdl[1][3];
      const double zc16 = pdm[1][4] * pdl[1][2];
      const double hc216 = pdm[1][5] * pdl[1][4];
      densities[1] *= ccor(z, rl, hc16, zc16, hc216);
      /*   Chemistry correction */
      const double hcc16 = pdm[1][7] * pdl[1][13];
      const double zcc16 = pdm[1][6] * pdl[1][12];
      const double rc16 = pdm[1][3] * pdl[1][14];
      /*  Net density corrected at Alt */
      densities[1] *= ccor(z, rc16, hcc16, zcc16);
    }
  }

  /**** O2 DENSITY ****/
  {
    /*   Density variation factor at Zlb */
    const double g32 = globe7(pd[4], tloc, doy, sec, flh.lat(), flh.lon(), f107, f107A, ap, plg);
    /*  Diffusive density at Zlb */
    const double db32 = pdm[3][0] * std::exp(g32) * pd[4][0];
    /*   Diffusive density at Alt */
    densities[3] = densu(z, db32, tinf, tlb, 32., alpha[3], temperatures[1],
                         ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
    if (z <= altl[3]) {
      /*   Turbopause */
      const double zh32 = pdm[3][2];
      /*  Mixed density at Zlb */
      const double b32 =
          densu(zh32, db32, tinf, tlb, 32. - xmm, alpha[3] - 1.,
                temperatures[1], ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
      /*  Mixed density at Alt */
      const double dm32 = densu(z, b32, tinf, tlb, xmm, 0., temperatures[1],
                                ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
      const double zhm32 = zhm28;
      /*  Net density at Alt */
      densities[3] = dnet(densities[3], dm32, zhm32, xmm, 32.);
      /*   Correction to specified mixing ratio at ground */
      const double rl = std::log(b28 * pdm[3][1] / b32);
      const double hc32 = pdm[3][5] * pdl[1][7];
      const double zc32 = pdm[3][4] * pdl[1][6];
      densities[3] *= ccor(z, rl, hc32, zc32);
    }
    /*  Correction for general departure from diffusive equilibrium above Zlb */
    const double hcc32 = pdm[3][7] * pdl[1][22];
    const double hcc232 = pdm[3][7] * pdl[0][22];
    const double zcc32 = pdm[3][6] * pdl[1][21];
    const double rc32 =
        pdm[3][3] * pdl[1][23] * (1. + pdl[0][23] * (f107A - 150.));
    /*  Net density corrected at Alt */
    densities[3] *= ccor(z, rc32, hcc32, zcc32, hcc232);
  }

  /**** AR DENSITY ****/
  {
    /*   Density variation factor at Zlb */
    const double g40 = globe7(pd[5], tloc, doy, sec, flh.lat(), flh.lon(), f107, f107A, ap, plg);
    /*  Diffusive density at Zlb */
    const double db40 = pdm[4][0] * std::exp(g40) * pd[5][0];
    /*   Diffusive density at Alt */
    densities[4] = densu(z, db40, tinf, tlb, 40., alpha[4], temperatures[1],
                         ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
    if (z <= altl[4]) {
      /*   Turbopause */
      const double zh40 = pdm[4][2];
      /*  Mixed density at Zlb */
      const double b40 =
          densu(zh40, db40, tinf, tlb, 40. - xmm, alpha[4] - 1.,
                temperatures[1], ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
      /*  Mixed density at Alt */
      const double dm40 = densu(z, b40, tinf, tlb, xmm, 0., temperatures[1],
                                ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
      const double zhm40 = zhm28;
      /*  Net density at Alt */
      densities[4] = dnet(densities[4], dm40, zhm40, xmm, 40.);
      /*   Correction to specified mixing ratio at ground */
      const double rl = log(b28 * pdm[4][1] / b40);
      const double hc40 = pdm[4][5] * pdl[1][9];
      const double zc40 = pdm[4][4] * pdl[1][8];
      /*  Net density corrected at Alt */
      densities[4] *= ccor(z, rl, hc40, zc40);
    }
  }

  /**** HYDROGEN DENSITY ****/
  {
    /*   Density variation factor at Zlb */
    const double g1 = globe7(pd[6], tloc, doy, sec, flh.lat(), flh.lon(), f107, f107A, ap, plg);
    /*  Diffusive density at Zlb */
    const double db01 = pdm[5][0] * std::exp(g1) * pd[6][0];
    /*   Diffusive density at Alt */
    densities[6] = densu(z, db01, tinf, tlb, 1., alpha[6], temperatures[1],
                         ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
    if (z <= altl[6]) {
      /*   Turbopause */
      const double zh01 = pdm[5][2];
      /*  Mixed density at Zlb */
      const double b01 =
          densu(zh01, db01, tinf, tlb, 1. - xmm, alpha[6] - 1., temperatures[1],
                ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
      /*  Mixed density at Alt */
      const double dm01 = densu(z, b01, tinf, tlb, xmm, 0., temperatures[1],
                                ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
      const double zhm01 = zhm28;
      /*  Net density at Alt */
      densities[6] = dnet(densities[6], dm01, zhm01, xmm, 1.);
      /*   Correction to specified mixing ratio at ground */
      const double rl =
          std::log(b28 * pdm[5][1] * sqrt(pdl[1][17] * pdl[1][17]) / b01);
      const double hc01 = pdm[5][5] * pdl[1][11];
      const double zc01 = pdm[5][4] * pdl[1][10];
      densities[6] *= ccor(z, rl, hc01, zc01);
      /*   Chemistry correction */
      const double hcc01 = pdm[5][7] * pdl[1][19];
      const double zcc01 = pdm[5][6] * pdl[1][18];
      const double rc01 = pdm[5][3] * pdl[1][20];
      /*  Net density corrected at Alt */
      densities[6] *= ccor(z, rc01, hcc01, zcc01);
    }
  }

  /**** ATOMIC NITROGEN DENSITY ****/
  {
    /*   Density variation factor at Zlb */
    const double g14 = globe7(pd[7], tloc, doy, sec, flh.lat(), flh.lon(), f107, f107A, ap, plg);
    /*  Diffusive density at Zlb */
    const double db14 = pdm[6][0] * std::exp(g14) * pd[7][0];
    /*   Diffusive density at Alt */
    densities[7] = densu(z, db14, tinf, tlb, 14., alpha[7], temperatures[1],
                         ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
    if (z <= altl[7]) {
      /*   Turbopause */
      const double zh14 = pdm[6][2];
      /*  Mixed density at Zlb */
      const double b14 =
          densu(zh14, db14, tinf, tlb, 14. - xmm, alpha[7] - 1.,
                temperatures[1], ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
      /*  Mixed density at Alt */
      const double dm14 = densu(z, b14, tinf, tlb, xmm, 0., temperatures[1],
                                ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
      const double zhm14 = zhm28;
      /*  Net density at Alt */
      densities[7] = dnet(densities[7], dm14, zhm14, xmm, 14.);
      /*   Correction to specified mixing ratio at ground */
      const double rl =
          std::log(b28 * pdm[6][1] * std::sqrt(pdl[0][2] * pdl[0][2]) / b14);
      const double hc14 = pdm[6][5] * pdl[0][1];
      const double zc14 = pdm[6][4] * pdl[0][0];
      densities[7] *= ccor(z, rl, hc14, zc14);
      /*   Chemistry correction */
      cosnt double hcc14 = pdm[6][7] * pdl[0][4];
      cosnt double zcc14 = pdm[6][6] * pdl[0][3];
      cosnt double rc14 = pdm[6][3] * pdl[0][5];
      /*  Net density corrected at Alt */
      densities[7] *= ccor(z, rc14, hcc14, zcc14);
    }
  }

  /**** Anomalous OXYGEN DENSITY ****/
  {
    const double g16h = globe7(pd[8], tloc, doy, sec, flh.lat(), flh.lon(), f107, f107A, ap, plg);
    const double db16h = pdm[7][0] * std::exp(g16h) * pd[8][0];
    const double tho = pdm[7][9] * pdl[0][6];
    const double dd = densu(z, db16h, tho, tho, 16., alpha[8], temperatures[1],
                            ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
    const double zsht = pdm[7][5];
    const double zmho = pdm[7][4];
    const double zsho = scalh(zmho, 16.0, tho);
    densities[8] = dd * std::exp(-zsht / zsho * (std::exp(-(z - zmho) / zsht) - 1.));
  }

  /* total mass density */
  densities[5] =
      1.66e-24 * (4.0 * densities[0] + 16.0 * densities[1] +
                  28.0 * densities[2] + 32.0 * densities[3] +
                  40.0 * densities[4] + densities[6] + 14.0 * densities[7]);

  /* temperature */
  return densu(altitude, 1.0, tinf, tlb, 0.0, 0.0, temperatures[1], ptm[5], s,
               mn1, zn1, meso_tn1, meso_tgn1);
}
