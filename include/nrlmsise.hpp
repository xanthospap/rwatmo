#ifndef __DSO_NRLMSISE00_UPPERATMO_HPP__
#define __DSO_NRLMSISE00_UPPERATMO_HPP__

#include "datetime/calendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "geodesy/transformations.hpp"
#include "geodesy/units.hpp"
#include <cassert>
#include <cmath>

namespace dso {

struct SpaceWeatherData {
  double f107;
  double f107A;
  double Ap;
  double ap[7];
}; /* SpaceWeatherData */

class Nrlmsise00 {

public:
  double density(Eigen::Vector3d rsat, dso::MjdEpoch &tt,
                 const SpaceWeatherData &data);

private:
  static constexpr const double rgas = 831.4;
  static constexpr const double dr = dso::D2PI / 365e0; // days to radians

  struct DataTrigs {
    /** @brief Holds trigs to avoid recomputing in glob7 and glob7s.
     *
     * @param[in] glat Geodetic Latitude in [rad]
     * @param[in] tloc Local solar time in hours (of day)
     */
    DataTrigs(double glat_, double tloc_) noexcept
        : glat(glat_), tloc(tloc_), c(std::sin(glat)), s(std::cos(glat)),
          stloc(std::sin(hours2rad(tloc))), ctloc(std::cos(hours2rad(tloc))),
          s2tloc(std::sin(2e0 * hours2rad(tloc))),
          c2tloc(std::cos(2e0 * hours2rad(tloc))),
          s3tloc(std::sin(3e0 * hours2rad(tloc))),
          c3tloc(std::cos(3e0 * hours2rad(tloc))) {}
    double glat, tloc, c, s, stloc, ctloc, s2tloc, c2tloc, s3tloc, c3tloc;
  }; /* DataTrigs */

  /* returns gv
   * computes rref
   */
  double glatf(double lat, double &reff) const noexcept {
    const double c2 = std::cos(2e0 * lat);
    const double gv = 980.616 * (1e0 - 0.0026373 * c2);
    reff = 2e0 * gv / (3.085462e-6 + 2.27e-9 * c2) * 1.0e-5;
    return gv;
  }

  double ccor(double alt, double r, double h1, double zh) const noexcept {
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

  double ccor(double alt, double r, double h1, double zh,
              double h2) const noexcept {
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

  double dnet(double dd, double dm, double zhm, double xmm,
              double xm) const noexcept {
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
    const double ylog = a * std::log(dm / dd);
    if (ylog < -10)
      return dd;
    if (ylog > 10)
      return dm;
    return dd * std::pow((1.0 + std::exp(ylog)), (1.0 / a));
  }

  double densu(double lat, double alt, double dlb, double tinf, double tlb,
               double xm, double alpha, double &tz, double zlb, double s2,
               int mn1, const double *zn1, double *tn1, double *tgn1,
               double *u) const noexcept;
  double gts7(const dso::MjdEpoch &t, const dso::GeodeticCrd &flh, double f107A,
              double f107, const double *const ap, double *densities,
              double *temperatures);
  double glob7(const double *const p, const DataTrigs &dt, double doy,
                double sec, double glon, double f107, double f107A,
                const double *const ap, double plg[4][9],
                double *apt) const noexcept;
  double glob7s(const double *p, DataTrigs &dt, double doy, double sec,
                 double glon, double f107, double f107A, const double *const apt,
                 const double plg[4][9]) const noexcept;

  /* POWER7 */
  static const double pt[150];
  static const double pd[9][150];
  static const double ps[150];
  static const double pdl[2][25];
  static const double ptl[4][100];
  static const double pma[10][100];
  static const double sam[100];
  /* LOWER7 */
  static const double ptm[10];
  static const double pdm[8][10];
  static const double pavgm[10];

}; /* Nrlmsise00 */
} /* namespace dso */

#endif
