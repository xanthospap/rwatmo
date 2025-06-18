#ifndef __DSO_NRLMSISE00_UPPERATMO_HPP__
#define __DSO_NRLMSISE00_UPPERATMO_HPP__
#include "datetime/calendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "geodesy/transformations.hpp"
#include "geodesy/units.hpp"
#include "space_weather.hpp"
#include <cassert>
#include <cmath>

namespace dso {

struct Msise00Data {
  double f107A; /* 81 day average of F10.7 flux (centered on doy) */
  double f107;  /* daily F10.7 flux for previous day */
  double ap[7];
  // 0 : daily AP
  // 1 : 3 hr AP index for current time
  // 2 : 3 hr AP index for 3 hrs before current time
  // 3 : 3 hr AP index for 6 hrs before current time
  // 4 : 3 hr AP index for 9 hrs before current time
  // 5 : Average of eight 3 hr AP indicies from 12 to 33 hrs
  //         prior to current time
  // 6 : Average of eight 3 hr AP indicies from 36 to 57 hrs
  //         prior to current time
};

class NrlmsiseDataHunter {
public:
  const Msise00Data &get_data(const MjdEpoch &tt) {
    if (this->hunt(tt)) {
      throw std::runtime_error(
          "[ERROR] Failed collecting MSISE space weather data.\n");
    }
    return _lastdata;
  }

  void set_data_ptr(const std::vector<SpaceWeatherData> &data_ptr) noexcept {
    _swdata = &data_ptr;
    _lastepoch = MjdEpoch::min();
  }

  const Msise00Data &msise_data() const noexcept { return _lastdata; }

private:
  /* seconds in a3-hour interval */
  static constexpr const long SEC_IN_3H = 60 * 60 * 3;
  /* ptr to a SpaceWeatherData vector where we shall extract data from */
  const std::vector<SpaceWeatherData> *_swdata{nullptr};
  /* current index in the SpaceWeatherData vector */
  int _ci = 0;
  /* this is the seconds (in day) of the epoch of which _ci referes to. Why do
   * we need this? Well, we are transforming epochs from UTC to TT, but we have
   * to deal with 3-hour intervals that start at 0H UTC. So, to find the 3-hour
   * interval an epoch belongs to, give that the reference epoch start at XX
   * seconds, we must subtract XX seconds.
   *
   * See the function _3hidx
   */
  double _ci_start_of_day_sec = 0e0;
  /* last epoch used to hunt data */
  MjdEpoch _lastepoch{MjdEpoch::min()};
  /* last Msise00Data acquired */
  Msise00Data _lastdata{};

  double start_of_day_sec() const noexcept {
    return (_swdata->cbegin() + _ci)->tt().seconds().seconds();
  }

  /* find index of 3-hour interval of given date, in range [0-8) */
  int _3hidx(const MjdEpoch &tt) const noexcept {
    const double sec_in_day = tt.seconds().seconds() - _ci_start_of_day_sec;
    const int idx = static_cast<long>(sec_in_day) / SEC_IN_3H;
    return idx;
  }

  /** @brief Collect Msise00Data for a given date, considering instance's state.
   *
   * By "considering instance's state" we mean that we will try to see if we
   * have already called the function requesting data for an epoch close to
   * the one we request here. If we can take advantage of that (re-use data
   * and/or indexes) we will, else we will eventually just call get_data to
   * do the job.
   *
   * This function virtually does the same as get_data but is optimized when
   * requesting often data for epochs that are close.
   *
   * Note that after the function call, _ci will be set to the index of the
   * record (within _swdata) that corresponds to the date tt, and _lastepoch
   * will be set to tt. The collected data will be stores in the instance's
   * _lastdata member.
   *
   * @param[in] tt The epoch of interest in [TT]
   * @return Anything other than zero denotes an error.
   */
  [[nodiscard]] int hunt(const dso::MjdEpoch &tt);

  /** @brief Fill in Msise00Data from an vector<SpaceWeatherData> given the
   * needed indexes.
   *
   * This function will fill in the values in the data instance (Msise00Data)
   * given that we wave a vector of SpaceWeatherData AND we know the index of
   * the current epoch (it) within this vector and the index of the current
   * 3-hour interval (apidx).
   *
   * Note that the data required may be up to 3 days prior to current date!
   * Make sure that the calling function that the passed in iterator is at
   * least >= std::vector<SpaceWeatherData>::cbegin()+3
   *
   * @param[in] it Const iterator to the entry (in an vector<SpaceWeatherData>,
   * normally _swdata) for the date of interest. You must have a very good
   * reason for this iterator to be outside the intsance's _swdata vector!
   * @param[in] apidx Index of the 3hour interval [0,7) for the time of
   * interest.
   * @param[out] data A pointer to an Msise00Data instance, where the
   * exctracted data will be stored at. If the pointer is NULL, then the
   * instance's _lastdata member will be used.
   * @return Always 0.
   */
  int fill_data(std::vector<SpaceWeatherData>::const_iterator it, int apidx,
                Msise00Data *data = nullptr);

  /** @brief Collect Msise00Data for a random date.
   *
   * This function will search through all of the elements of the _swdata
   * vector to find a matching date, and collect the relevant data.
   *
   * Note that after the function call, _ci will be set to the index of the
   * record (within _swdata) that corresponds to the date tt, and _lastepoch
   * will be set to tt.
   *
   * @param[in] tt The epoch of interest in [TT]
   * @param[in] mdata A pointer to an Msise00Data data, where the exctracted
   * data will be stored at. If the pointer is NULL, then the instance's
   * _lastdata member will hold the data.
   * @return Anything other than 0 denotes an error.
   */
  [[nodiscard]] int pr_get_data(const MjdEpoch &tt,
                                Msise00Data *mdata = nullptr);
}; /* NrlmsiseHunter */

class Nrlmsise00 {

public:
  double density(Eigen::Vector3d rsat_ecef, const dso::MjdEpoch &tt,
                 const Msise00Data &data);
  double density(Eigen::Vector3d rsat_ecef, const dso::MjdEpoch &tt) {
    if (mhunter)
      return density(rsat_ecef, tt, mhunter->get_data(tt));
    else
      throw std::runtime_error(
          "[ERROR] NrlmsiseDataHunter not setup, cannot compute density!\n");
  }

  NrlmsiseDataHunter &
  setup_data_hunter(const std::vector<SpaceWeatherData> &data) {
    if (mhunter)
      delete mhunter;
    mhunter = new NrlmsiseDataHunter();
    mhunter->set_data_ptr(data);
    return *mhunter;
  }

  ~Nrlmsise00() noexcept {
    if (mhunter)
      delete mhunter;
  }

private:
  static constexpr const double rgas = 831.4;
  static constexpr const double dr = dso::D2PI / 365e0; // days to radians
  NrlmsiseDataHunter *mhunter = nullptr;

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
  double glob7s(const double *p, DataTrigs &dt, double doy,
                double glon, double f107A, const double *const apt,
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
