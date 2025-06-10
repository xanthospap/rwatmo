#include "datetime/datetime_read.hpp"
#include <charconv>
#include <cstring>
#include <stdexcept>

/* ---------------------------------------------------------------------------
[ 0] DATE	Year-Month-Day (ISO 8601)
[ 1] BSRN	Bartels Solar Rotation Number. A sequence of 27-day intervals 
                counted continuously from 1832 Feb 8.
[ 2] ND	        Number of Day within the Bartels 27-day cycle (01-27).
[ 3] KP1	Planetary 3-hour Range Index (Kp) for 0000-0300 UT.
[ 4] KP2	Planetary 3-hour Range Index (Kp) for 0300-0600 UT.
[ 5] KP3	Planetary 3-hour Range Index (Kp) for 0600-0900 UT.
[ 6] KP4	Planetary 3-hour Range Index (Kp) for 0900-1200 UT.
[ 7] KP5	Planetary 3-hour Range Index (Kp) for 1200-1500 UT.
[ 8] KP6	Planetary 3-hour Range Index (Kp) for 1500-1800 UT.
[ 9] KP7	Planetary 3-hour Range Index (Kp) for 1800-2100 UT.
[10] KP8	Planetary 3-hour Range Index (Kp) for 2100-0000 UT.
[11] KP_SUM	Sum of the 8 Kp indices for the day.
                Kp has values of 0o, 0+, 1-, 1o, 1+, 2-, 2o, 2+, ... , 8o, 8+, 
                9-, 9o, which are expressed in steps of one third unit. These 
                values are multiplied by 10 and rounded to an integer value.
[12] AP1	Planetary Equivalent Amplitude (Ap) for 0000-0300 UT.
[13] AP2	Planetary Equivalent Amplitude (Ap) for 0300-0600 UT.
[14] AP3	Planetary Equivalent Amplitude (Ap) for 0600-0900 UT.
[15] AP4	Planetary Equivalent Amplitude (Ap) for 0900-1200 UT.
[16] AP5	Planetary Equivalent Amplitude (Ap) for 1200-1500 UT.
[17] AP6	Planetary Equivalent Amplitude (Ap) for 1500-1800 UT.
[18] AP7	Planetary Equivalent Amplitude (Ap) for 1800-2100 UT.
[19] AP8	Planetary Equivalent Amplitude (Ap) for 2100-0000 UT.
[20] AP_AVG	Arithmetic average of the 8 Ap indices for the day.
[21] CP	        Cp or Planetary Daily Character Figure. A qualitative estimate 
                of overall level of magnetic activity for the day determined 
                from the sum of the 8 Ap indices. Cp ranges, in steps of 
                one-tenth, from 0 (quiet) to 2.5 (highly disturbed).
[22] C9	        C9. A conversion of the 0-to-2.5 range of the Cp index to one 
                digit between 0 and 9.
[23] ISN	International Sunspot Number. Records contain the Zurich 
                number through 1980 Dec 31 and the International Brussels 
                number thereafter.
[24] F10.7_OBS	Observed 10.7-cm Solar Radio Flux (F10.7). Measured at Ottawa 
                at 1700 UT daily from 1947 Feb 14 until 1991 May 31 and 
                measured at Penticton at 2000 UT from 1991 Jun 01 on. 
                Expressed in units of 10-22 W/m2/Hz.
[25] F10.7_ADJ	10.7-cm Solar Radio Flux (F10.7) adjusted to 1 AU.
[26] F10.7_DATA_TYPE	Flux Qualifier.
                OBS: Observed flux measurement
                INT: CelesTrak linear interpolation of missing data
                PRD: 45-Day predicted flux
                PRM: Monthly predicted flux
[27] F10.7_OBS_CENTER81	Centered 81-day arithmetic average of F10.7 (observed).
[28] F10.7_OBS_LAST81	Last 81-day arithmetic average of F10.7 (observed).
[29] F10.7_ADJ_CENTER81	Centered 81-day arithmetic average of F10.7 (adjusted).
[30] F10.7_ADJ_LAST81	Last 81-day arithmetic average of F10.7 (adjusted).
---------------------------------------------------------------------------- */

// 27 f107A
// 28 f107 daily F10.7 flux for previous day
// 12-20

class SpaceWeatherData {
  MjdEpoch tt;
  int iar[];
  double dar[];
}; /* class SpaceWeatherData */


int parse_line(const char *line, dso::TwoPartDateUTC &utc, int *__restrict__ intar, double *__restrict__ fltar) noexcept {
  try {
    utc = dso::TwoPartDateUTC(dso::ReadInDate<dso::YMDFormat::YYYYMMDD>::read(line, nullptr));
  } catch (std::exception &e) {
    fprintf(stderr, "[ERROR] Failed resolvding date from Clestrak SpaceWaether data; line was %s (traceback: %s)\n", line, __func__);
    return 1;
  }

  int error = 0;
  int fltidx = 0;
  int intidx = 0;

  const char *str = line + 11;
  const char *last = line + std::strlen(line) - 1;
  for (int i=1; i<21; i++, intidx++) {
    const auto res = std::from_chars(str, last, intar[intidx]);
    if (res.ec != std::errc{}) ++error;
    str = last + 1;
  }
  
  {
    const auto res = std::from_chars(str, last, fltar[fltidx]);
    if (res.ec != std::errc{}) ++error;
    str = last + 1;
    ++fltidx;
  }
  
  for (int i=0; i<2; i++, intidx++) {
    const auto res = std::from_chars(str, last, intar[intidx]);
    if (res.ec != std::errc{}) ++error;
    str = last + 1;
  }
  
  for (int i=24; i<31; i++, fltidx++) {
    const auto res = std::from_chars(str, last, fltar[fltidx]);
    if (res.ec != std::errc{}) ++error;
    str = last + 1;
  }

  return error;
}
