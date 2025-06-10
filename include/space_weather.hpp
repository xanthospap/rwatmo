#ifndef __DSO_RWATMO_SPACEWEATHER_VAR_HPP__
#define __DSO_RWATMO_SPACEWEATHER_VAR_HPP__

#include "datetime/calendar.hpp"
#include <vector>

namespace dso {
/* ---------------------------------------------------------------------------
 *[ 0| -| -] DATE	Year-Month-Day (ISO 8601)
 *[ 1| 0| -] BSRN	Bartels Solar Rotation Number. A sequence of 27-day
 *                  intervals counted continuously from 1832 Feb 8.
 *[ 2| 1| -] ND	    Number of Day within the Bartels 27-day cycle (01-27).
 *[ 3| 2| -] KP1	Planetary 3-hour Range Index (Kp) for 0000-0300 UT.
 *[ 4| 3| -] KP2	Planetary 3-hour Range Index (Kp) for 0300-0600 UT.
 *[ 5| 4| -] KP3	Planetary 3-hour Range Index (Kp) for 0600-0900 UT.
 *[ 6| 5| -] KP4	Planetary 3-hour Range Index (Kp) for 0900-1200 UT.
 *[ 7| 6| -] KP5	Planetary 3-hour Range Index (Kp) for 1200-1500 UT.
 *[ 8| 7| -] KP6	Planetary 3-hour Range Index (Kp) for 1500-1800 UT.
 *[ 9| 8| -] KP7	Planetary 3-hour Range Index (Kp) for 1800-2100 UT.
 *[10| 9| -] KP8	Planetary 3-hour Range Index (Kp) for 2100-0000 UT.
 *[11|10| -] KP_SUM	Sum of the 8 Kp indices for the day. Kp has values of
 *                  0o, 0+, 1-, 1o, 1+, 2-, 2o, 2+, ... , 8o, 8+, 9-, 9o,
 *                  which are expressed in steps of one third unit. These
 *                  values are multiplied by 10 and rounded to an integer value.
 *[12|11| -] AP1	Planetary Equivalent Amplitude (Ap) for 0000-0300 UT.
 *[13|12| -] AP2	Planetary Equivalent Amplitude (Ap) for 0300-0600 UT.
 *[14|13| -] AP3	Planetary Equivalent Amplitude (Ap) for 0600-0900 UT.
 *[15|14| -] AP4	Planetary Equivalent Amplitude (Ap) for 0900-1200 UT.
 *[16|15| -] AP5	Planetary Equivalent Amplitude (Ap) for 1200-1500 UT.
 *[17|16| -] AP6	Planetary Equivalent Amplitude (Ap) for 1500-1800 UT.
 *[18|17| -] AP7	Planetary Equivalent Amplitude (Ap) for 1800-2100 UT.
 *[19|18| -] AP8	Planetary Equivalent Amplitude (Ap) for 2100-0000 UT.
 *[20|19| -] AP_AVG	Arithmetic average of the 8 Ap indices for the day.
 *[21| -| 0] CP	    Cp or Planetary Daily Character Figure. A qualitative
 *                  estimate of overall level of magnetic activity for the day
 *                  determined from the sum of the 8 Ap indices. Cp ranges, in
 *                  steps of one-tenth, from 0 (quiet) to 2.5 (highly
 *                  disturbed).
 *[22|20| -] C9	    C9. A conversion of the 0-to-2.5 range of the Cp index to
 *                  one digit between 0 and 9.
 *[23|21| -] ISN	International Sunspot Number. Records contain the Zurich
 *                  number through 1980 Dec 31 and the International Brussels
 *                  number thereafter.
 *[24| -| 1] F10.7_OBS	Observed 10.7-cm Solar Radio Flux (F10.7). Measured at
 *                  Ottawa at 1700 UT daily from 1947 Feb 14 until 1991 May 31
 *                  and measured at Penticton at 2000 UT from 1991 Jun 01 on.
 *                  Expressed in units of 10-22 W/m2/Hz.
 *[25| -| 2] F10.7_ADJ	10.7-cm Solar Radio Flux (F10.7) adjusted to 1 AU.
 *[26| -| -] F10.7_DATA_TYPE	Flux Qualifier.
 *                  OBS: Observed flux measurement
 *                  INT: CelesTrak linear interpolation of missing data
 *                  PRD: 45-Day predicted flux
 *                  PRM: Monthly predicted flux
 *[27| -| 3] F10.7_OBS_CENTER81	Centered 81-day arithmetic average of F10.7
 *                  (observed).
 *[28| -| 4] F10.7_OBS_LAST81	Last 81-day arithmetic average of F10.7
 *                  (observed).
 *[29| -| 5] F10.7_ADJ_CENTER81	Centered 81-day arithmetic average of F10.7
 *                  (adjusted).
 *[30| -| 6] F10.7_ADJ_LAST81	Last 81-day arithmetic average of F10.7
 *                 (adjusted).
---------------------------------------------------------------------------- */
class SpaceWeatherData {
  MjdEpoch mtt;
  char mfq[4];
  int miar[22];
  double mdar[7];

public:
  /* 81 day average of F10.7 flux (centered on doy) */
  double f107A() const noexcept { return mdar[3]; }
  /* daily F10.7 flux for current day */
  double f107() const noexcept { return mdar[0]; }
  static constexpr const int _3hapindex_start = 11;

public:
  const MjdEpoch &tt() const noexcept { return mtt; }
  MjdEpoch &tt() noexcept { return mtt; }
  const int *int_array() const noexcept { return miar + 0; }
  int *int_array() noexcept { return miar + 0; }
  const double *flt_array() const noexcept { return mdar + 0; }
  double *flt_array() noexcept { return mdar + 0; }
  const char *flux_qualifier() const noexcept { return mfq + 0; }
  char *flux_qualifier() noexcept { return mfq + 0; }

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
  class NrlmsiseHunter {
    /* current index in vector */
    int ci = 0;

    const long SEC_IN_3H = 60 * 60 * 3;
    int _3hi = -1;
    MjdEpoch _lastepoch;
    Msise00Data _lastdata;

    int update_by_3h(const MjdEpoch &tt,
                     const std::vector<SpaceWeatherData> &sdata,
                     Msise00Data &mdata) {
      if (tt.imjd() != _lastepoch.imjd()) {
        ++ci;
        auto it = sdata.cbegin() + ci;
        if ((tt.imjd() != it->tt().imjd()) || (_3hidx(tt) != 0))
      }

      return 1;
    }

    int get_msise00_data(const MjdEpoch &tt,
                         const std::vector<SpaceWeatherData> &sdata,
                         Msise00Data &mdata) {
      /* search for day of interest in array */
      auto it = sdata.cbegin();
      for (it; it != sdata.cend(); ++it) {
        if (tt.diff<DateTimeDifferenceType::FractionalDays>(it->tt()).days() <
            1.) {
          /* we must have data for two days prior to current */
          if (std::distance(sdata.cbegin(), it) < 2) {
            fprintf(stderr,
                    "[ERROR] Not enough space weather data (prior to current "
                    "date) for NRLMSISE00 model! (traceback: %s)\n",
                    __func__);
            return 1;
          }
          /* update current index */
          ci = std::distance(sdata.cbegin(), it);
          mdata.f107A = it->f107A();
          mdata.f107 = (it - 1)->f107();
          mdata.ap[0] = it->int_array()[19];
          /* current 3-hour interval within the day */
          int idx = _3hidx(tt);
          /* next three are the prior 3hour intervals */
          mdata.ap[1] = it->int_array()[_3hapindex_start + idx];
          for (int i = 0; i < 3; i++) {
            if (idx == 0) {
              --it;
              idx = 7;
            } else {
              --idx;
            }
            mdata.ap[1 + i] = it->int_array()[_3hapindex_start + idx];
          }
          /* next, Average of eight 3 hr AP indicies from 12 to 33 hrs prior to
           * current time */
          double average = 0;
          for (int i = 0; i < 8; i++) {
            if (idx == 0) {
              --it;
              idx = 7;
            } else {
              --idx;
            }
            average += it->int_array()[_3hapindex_start + idx];
          }
          mdata.ap[5] = average / 8.;
          /* next, Average of eight 3 hr AP indicies from 36 to 57 hrs prior to
           * current time */
          for (int i = 0; i < 8; i++) {
            if (idx == 0) {
              --it;
              idx = 7;
            } else {
              --idx;
            }
            average += it->int_array()[_3hapindex_start + idx];
          }
          mdata.ap[6] = average / 8.;

          return 0;
        }
      }

      return 1;
    }

    /* find index of 3-hour interval of given date, in range [0-8) */
    int _3hidx(const MjdEpoch &tt) noexcept {
      const double sec_in_day = tt.seconds().seconds();
      const int idx = static_cast<long>(sec_in_day) / SEC_IN_3H;
      return idx;
    }

    int get_data_for(const MjdEpoch &tt,
                     const std::vector<SpaceWeatherData> &data) {
      if ((_3hi >= 0) &&
          ((tt.imjd() == _lastepoch.imjd()) && (_3hi == index(tt)))) {
        return _lastdata;
      } else if ((_3hi >= 0) &&
                 ((tt.imjd() == _lastepoch.imjd()) &&
                  (tt.diff<DateTimeDifferenceType::FractionalSeconds>(
                         _lastepoch)
                       .seconds() < SEC_IN_3H))) {
      }
    }
  }; /* NrlmsiseHunter */
}; /* class SpaceWeatherData */

/** @brief Collect space weather data from Celestrak data file (new
 * non-legacy).
 *
 * Collect space weather data, as given in the Celestrak CSV new data file
 * (not tested with the legacy format).
 *
 * @param[in] fn Filename of the Celestrak CSV file (see
 * https://celestrak.org/SpaceData/SpaceWx-format.php)
 * @param[in] tstart Starting epoch, inclusive, in TT
 * @param[in] tend   End epoch, not inclusive, in TT
 * @return A vector of SpaceWeatherData instances, one record per day.
 */
std::vector<SpaceWeatherData>
load_celestrak_sw(const char *fn, const dso::MjdEpoch &tstart,
                  const dso::MjdEpoch &tend) noexcept;
} /* namespace dso */

#endif