#include "datetime/datetime_read.hpp"
#include "space_weather.hpp"
#include <charconv>
#include <cstring>
#include <fstream>
#include <stdexcept>

namespace {
dso::MjdEpoch parse_date(const char *line, int &error) noexcept {
  try {
    const auto utc = dso::TwoPartDateUTC(
        dso::ReadInDate<dso::YMDFormat::YYYYMMDD>::read(line, nullptr));
    return utc.utc2tt();
  } catch (std::exception &e) {
    fprintf(stderr,
            "[ERROR] Failed resolving date from Clestrak SpaceWaether data; "
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
    str = last + 1;
  }

  {
    const auto res = std::from_chars(str, last, fltar[fltidx]);
    if (res.ec != std::errc{})
      ++error;
    str = last + 1;
    ++fltidx;
  }

  for (int i = 0; i < 2; i++, intidx++) {
    const auto res = std::from_chars(str, last, intar[intidx]);
    if (res.ec != std::errc{})
      ++error;
    str = last + 1;
  }

  for (int i = 24; i < 26; i++, fltidx++) {
    const auto res = std::from_chars(str, last, fltar[fltidx]);
    if (res.ec != std::errc{})
      ++error;
    str = last + 1;
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
    str = last + 1;
  }

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
  dso::MjdEpoch tc;

  /* keep on looking untill we meet the first date of interest */
  while (fin.getline(line, 256) && (!error)) {
    tc = parse_date(line, error);
  }

  /* collect data for dates of interest */
  if ((!error) && fin.good()) {
    dso::SpaceWeatherData data;
    error += parse_line(line, &data);
    res.emplace_back(data);
    while (fin.getline(line, 256) && (!error)) {
      error += parse_line(line, &data);
      if (data.tt() >= tend)
        break;
      res.emplace_back(data);
    }
  }

  if (error && fin.eof()) {
    return std::vector<dso::SpaceWeatherData>{};
  }

  return res;
}