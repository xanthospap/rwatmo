#include "datetime/datetime_read.hpp"
#include "space_weather.hpp"
#include <charconv>
#include <cstring>
#include <fstream>
#include <stdexcept>
#ifdef DEBUG
#include <cassert>
#endif

namespace {
dso::MjdEpoch parse_date(const char *line, int &error) noexcept {
  try {
    const auto utc = dso::TwoPartDateUTC(
        dso::ReadInDate<dso::YMDFormat::YYYYMMDD>::read(line, nullptr));
    /* important! this willset hours, minutes, seconds (all) to zero, which may
     * be misleading when transforming to other time scales. Set hours of day to
     * 12h to avoid this.*/
    return utc.add_seconds(dso::FractionalSeconds(86400e0 / 2)).utc2tt();
  } catch (std::exception &e) {
    fprintf(stderr,
            "[ERROR] Failed resolving date from Clestrak SpaceWeather data; "
            "line was %s (traceback: %s)\n",
            line, __func__);
    error += 1;
  }

  return dso::MjdEpoch::min();
}

int parse_line(const char *line, dso::SpaceWeatherData *data) noexcept {
  try {
    const auto utc = dso::TwoPartDateUTC(
        dso::ReadInDate<dso::YMDFormat::YYYYMMDD>::read(line, nullptr));
    data->tt() = utc.utc2tt();
  } catch (std::exception &e) {
    fprintf(stderr,
            "[ERROR] Failed resolvding date from Clestrak SpaceWaether data; "
            "line was %s (traceback: %s)\n",
            line, __func__);
    return 1;
  }

  int error = 0;
  int fltidx = 0;
  int intidx = 0;

  int *__restrict__ intar = data->int_array();
  double *__restrict__ fltar = data->flt_array();

  const char *str = line + 11;
  const char *last = line + std::strlen(line) - 1;
  for (int i = 1; i < 21; i++, intidx++) {
    const auto res = std::from_chars(str, last, intar[intidx]);
    if (res.ec != std::errc{})
      ++error;
    // printf("int[%2d]=%d\n",intidx,intar[intidx]);
    str = res.ptr + 1;
  }

  {
    const auto res = std::from_chars(str, last, fltar[fltidx]);
    if (res.ec != std::errc{})
      ++error;
    str = res.ptr + 1;
    ++fltidx;
  }

  for (int i = 0; i < 2; i++, intidx++) {
    const auto res = std::from_chars(str, last, intar[intidx]);
    if (res.ec != std::errc{})
      ++error;
    str = res.ptr + 1;
  }

  for (int i = 24; i < 26; i++, fltidx++) {
    const auto res = std::from_chars(str, last, fltar[fltidx]);
    if (res.ec != std::errc{})
      ++error;
    str = res.ptr + 1;
  }

  {
    std::strncpy(data->flux_qualifier(), str, 3);
    data->flux_qualifier()[3] = '\0';
    error += !(str[3] == ',');
    str += 4;
  }

  for (int i = 27; i < 31; i++, fltidx++) {
    const auto res = std::from_chars(str, last, fltar[fltidx]);
    if (res.ec != std::errc{})
      ++error;
    str = res.ptr + 1;
  }

#ifdef DEBUG
  assert(fltidx) == 7;
  assert(intidx) == 22;
#endif
  return error;
}
} /* anonymous namespace */

std::vector<dso::SpaceWeatherData>
dso::load_celestrak_sw(const char *fn, const dso::MjdEpoch &tstart,
                       const dso::MjdEpoch &tend) noexcept {
  std::vector<dso::SpaceWeatherData> res;
  res.reserve(static_cast<int>(
      tend.diff<dso::DateTimeDifferenceType::FractionalDays>(tstart).days()));

  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr,
            "[ERROR] Failed opening Celestrak space weather data file %s "
            "(traceback: %s)\n",
            fn, __func__);
    return res;
  }

  char line[256];
  int error = 0;

  /* first line is just the header */
  fin.getline(line, 256);

  /* keep on looking untill we meet the first date of interest */
  while (fin.getline(line, 256) && (!error)) {
    if (parse_date(line, error).imjd() == tstart.imjd())
      break;
  }

  dso::SpaceWeatherData data;
  /* collect data for dates of interest */
  if ((!error) && fin.good()) {
    /* first store line in buffer */
    error += parse_line(line, &data);
    res.emplace_back(data);
    /* keep on storing data until we meet tend */
    while (fin.getline(line, 256) && (!error)) {
      error += parse_line(line, &data);
      if (data.tt() >= tend)
        break;
      res.emplace_back(data);
    }
  }

  if (error) {
    fprintf(stderr,
            "[ERROR] Something went wrong while parsing space weather data; "
            "error code: %d, last line read: %s (traceback: %s)\n",
            error, line, __func__);
    return std::vector<dso::SpaceWeatherData>{};
  }

  return res;
}
