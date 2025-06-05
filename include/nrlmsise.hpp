#ifndef __DSO_NRLMSISE00_UPPERATMO_HPP__
#define __DSO_NRLMSISE00_UPPERATMO_HPP__

#include <cmath>
#include <cassert>


namespace dso {

namespace nrlmsise_impl {
double globe7(double *p, double tloc, double doy, double sec, double glat, double glon, double f107, double f107A, const double *const ap) noexcept;
double globe7s(double *p, double tloc, double doy, double sec, double glat, double glon, double f107, double f107A, const double *const ap, const double* plg[]) noexcept;
} /* namespace nrlmsise_impl */

class SpaceWeatherData {
  double f107;
  double f107A;
  double Ap;
  double ap[7];
}; /* SpaceWeatherData */

class Nrlmsise00 {

  double density(Eigen::Vector3d rsat, dso::MjdEpoch &tt, const SpaceWeatherData &data);

  /* returns gv
   * computes rref
   */
  static double glatf(double lat, double &reff) noexcept {
    const double dgtr = 1.74533e-2;
    const double c2 = std::cos(2e0 * dgtr * lat);
    const double gv = 980.616 * (1e0 - 0.0026373 * c2);
    reff = 2e0 * gv / (3.085462e-6 + 2.27e-9 * c2) * 1.0e-5;
    return gv;
  }

  static double ccor(double alt, double r, double h1, double zh) noexcept {
    /* chemistry/dissociation correction for msis models
     *  alt - altitude
     *  r - target ratio
     *  h1 - transition scale length
     *  zh - altitude of 1/2 r
     */
    const double e = (alt - zh) / h1;
    if (e > 70) {
      return 1e0;
    }
    if (e < -70) {
      return std::exp(r);
    }
    return std::exp(r / (1.0 + std::exp(e)));
  }

  static double ccor(double alt, double r, double h1, double zh,
                     double h2) noexcept {
    /* chemistry/dissociation correction for msis models
     *  alt - altitude
     *  r - target ratio
     *  h1 - transition scale length
     *  zh - altitude of 1/2 r
     *  h2 - transition scale length #2 ?
     */
    const double e1 = (alt - zh) / h1;
    const double e2 = (alt - zh) / h2;
    if ((e1 > 70) || (e2 > 70)) {
      return 1e0;
    }
    if ((e1 < -70) && (e2 < -70)) {
      return std::exp(r);
    }
    return std::exp(r / (1e0 + 5e-1 * (std::exp(e1) + std::exp(e2))));
  }

  static double scalh(double alt, double xm, double temp) noexcept {
    constexpr const double rgas = 831.4;
    return rgas * temp / ((gsurf / (std::pow((1e0 + alt / re), 2e0))) * xm);
  }

  static double dnet(double dd, double dm, double zhm, double xmm,
                     double xm) noexcept {
    /*       TURBOPAUSE CORRECTION FOR MSIS MODELS
     *        Root mean density
     *         DD - diffusive density
     *         DM - full mixed density
     *         ZHM - transition scale length
     *         XMM - full mixed molecular weight
     *         XM  - species molecular weight
     *         DNET - combined density
     */
    assert((dm > 0) && (dd > 0));
    const double a = zhm / (xmm - xm);
    const double ylog = a * log(dm / dd);
    if (ylog < -10)
      return dd;
    if (ylog > 10)
      return dm;
    return dd * std::pow((1.0 + std::exp(ylog)), (1.0 / a));
  }

  double splini(const double *__restrict__ xa, const double *__restrict__ ya,
                const double *__restrict__ y2a, int n, double x) noexcept;
  double splint(const double *__restrict__ xa, const double *__restrict__ ya,
                const double *__restrict__ y2a, int n, double x) noexcept;
  void spline(const double *__restrict__ x, const double *__restrict__ y, int n,
              double yp1, double ypn, const double *__restrict__ y2,
              double *__restrict__ u) noexcept;
  double densu(double alt, double dlb, double tinf, double tlb, double xm,
               double alpha, double &tz, double zlb, double s2, int mn1,
               const double *zn1, double *tn1, double *tgn1);
  double densm(double alt, double d0, double xm, double &tz, int mn3,
               const double *zn3, const double *tn3, const double *tgn3,
               int mn2, const double *const zn2, const double *tn2,
               const double *tgn2, double *__restrict__ u) noexcept;
               
  double gts7(const dso::MjdEpoch &t, double lon, double lat, double altitude, double f107A, double f107, const double *const ap, double *densities, double *temperatures);
  void gtd7(const dso::MjdEpoch& tt, const dso::GeodeticCrd &flh,double f107A, double f107, const double* const ap, double* densities, double* temperatures);

  /* POWER7 */
  static constexpr double pt[150]; // used by glob7
  static constexpr double pd[9][150]; // 
  static constexpr double ps[150]; // used by glob7
  static constexpr double pdl[2][25]; // used by gts7
  static constexpr double ptl[4][100];
  static constexpr double pma[10][100];
  static constexpr double sam[100];
  /* LOWER7 */
  static constexpr double ptm[10];
  static constexpr double pdm[8][10];
  static constexpr double pavgm[10];

  /* TODO why are these here? */
  double gsurf;
  double re;
  /* GTS3C */
  double dd;
  /* DMIX */
  double dm04, dm16, dm28, dm32, dm40, dm01, dm14;
  /* LPOLY */
  double dfa;
  double plg[4][9];
  double ctloc, stloc;
  double c2tloc, s2tloc;
  double s3tloc, c3tloc;
  double apdf, apt[4];
}; /* Nrlmsise00 */
} /* namespace dso */

#endif
