#include <cassert>
#include <charconv>
#include <cstdio>
#include <cstring>
#include <fstream>

#include "datetime/datetime_read.hpp"
#include "vmf3_grid_data.hpp"

namespace {
constexpr const int MAX_VMF3_CHARS = 256;
constexpr const char *Data_types =
    "V3GR (lat lon ah aw zhd zwd Gn_h Ge_h Gn_w Ge_w)";
const char *skipws(const char *line) noexcept {
  while (*line && *line == ' ') ++line;
  return line;
}
}  // namespace

int dso::load_vmfgr3_grid_map(const char *fn,
                              dso::vmf3::GridVmf3Data *grid) noexcept {
  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr,
            "[ERROR] Failed opening VMF3 grid data file %s (traceback: %s)\n",
            fn, __func__);
    return 1;
  }

  char line[MAX_VMF3_CHARS];
  int error = 0;
  dso::MjdEpoch t = dso::MjdEpoch::min();
  double scale_factor = 1e0;
  double range[6];

  /* read header
   * ! Version:            1.0
   * ! Source:             D. Landskron, TU Vienna (created: 2021-01-04)
   * ! Data_types:         V3GR (lat lon ah aw zhd zwd Gn_h Ge_h Gn_w Ge_w)
   * ! Epoch:              2021 01 03 06 00  0.0
   * ! Scale_factor:       1.e+00
   * ! Range/resolution:   -87.5 87.5 2.5 357.5 5 5
   * ! Comment: vmf.geo.tuwien.ac.at/trop_products/GRID/5x5/V3GR/V3GR_OP/2021/
   */
  while (fin.getline(line, MAX_VMF3_CHARS) && (!error)) {
    if (line[0] != '!') break;

    if (!std::strncmp(line + 2, "Version:", 8)) {
      ;
    } else if (!std::strncmp(line + 2, "Source:", 7)) {
      ;
    } else if (!std::strncmp(line + 2, "Data_types:", 11)) {
      if (std::strncmp(skipws(line + 2 + 11), Data_types,
                       std::strlen(Data_types))) {
        fprintf(stderr,
                "[ERROR] Invalid(?) header line encountered in VMF3 grid data "
                "file %s (traceback: %s)\n",
                fn, __func__);
        fprintf(stderr,
                "[ERROR] C'nued erronuous line is: %s (traceback: %s)\n", line,
                __func__);
        error = 1;
      }
    } else if (!std::strncmp(line + 2, "Epoch:", 6)) {
      try {
        t = dso::from_char<dso::YMDFormat::YYYYMMDD, dso::HMSFormat::HHMMSSF>(
            line + 2 + 6);
      } catch (std::exception &) {
        fprintf(stderr,
                "[ERROR] Invalid(?) header line encountered in VMF3 grid data "
                "file %s (traceback: %s)\n",
                fn, __func__);
        fprintf(stderr,
                "[ERROR] C'nued erronuous line is: %s (traceback: %s)\n", line,
                __func__);
        error = 1;
      }
    } else if (!std::strncmp(line + 2, "Scale_factor:", 13)) {
      auto ec = std::from_chars(skipws(line + 2 + 13), line + std::strlen(line),
                                scale_factor);
      if (!(ec.ec == std::errc{}) || scale_factor != 1e0) {
        fprintf(stderr,
                "[ERROR] Invalid(?) header line encountered in VMF3 grid data "
                "file %s (traceback: %s)\n",
                fn, __func__);
        fprintf(stderr,
                "[ERROR] C'nued erronuous line is: %s (traceback: %s)\n", line,
                __func__);
        error = 1;
      }
    } else if (!std::strncmp(line + 2, "Range/resolution:", 17)) {
      const char *snum = line + 2 + 17;
      for (int i = 0; i < 6; i++) {
        auto ec =
            std::from_chars(skipws(snum), line + std::strlen(line), range[i]);
        if (!(ec.ec == std::errc{})) {
          fprintf(stderr,
                  "[ERROR] Invalid(?) header line encountered in VMF3 grid "
                  "data file %s (traceback: %s)\n",
                  fn, __func__);
          fprintf(stderr,
                  "[ERROR] C'nued erronuous line is: %s (traceback: %s)\n",
                  line, __func__);
          error = 1;
        }
        snum = ec.ptr;
      }
    } else if (!std::strncmp(line + 2, "Comment:", 8)) {
      ;
    } else {
      fprintf(stderr,
              "[ERROR] Invalid(?) header line encountered in VMF3 grid data "
              "file %s (traceback: %s)\n",
              fn, __func__);
      fprintf(stderr, "[ERROR] C'nued erronuous line is: %s (traceback: %s)\n",
              line, __func__);
      error = 1;
    }
  } /* end reading header */

  /* I do not know why this happens, but it seems that in all grid data files,
   * while we expect a range from -89.5 to 89.5 with e.g a step of 2.5 deg.,
   * the data records actually span the inverse range, i.e 89.5 to -89.5 with
   * a step of -2.5 deg.
   * Fix that!
   */
  if (range[0] < 0 && range[1] > 0) {
    std::swap(range[0], range[1]);
    range[4] *= -1e0;
  }

  /* set correct dimensions for the grid data we shall load/store (whatch out fo
   * the order!) */
  grid->set_range(range[0], range[1], range[4], range[2], range[3], range[5]);

  /* ready to start reading data; note that the first data line is buffered
   * in line[]. */
  std::size_t idx = 0;
  while (!(error) && ((idx < (std::size_t)grid->num_pts() && (fin.good())))) {
    const char *num = line;
    /* latitude and longitude */
    auto ec = std::from_chars(skipws(num), line + std::strlen(line),
                              grid->data(idx)->lat_deg());
    if (!(ec.ec == std::errc{})) ++error;
    num = ec.ptr;
    ec = std::from_chars(skipws(num), line + std::strlen(line),
                         grid->data(idx)->lon_deg());
    if (!(ec.ec == std::errc{})) ++error;
    num = ec.ptr;

    /* VMF3 grid data */
    for (int i = 0; i < dso::vmf3::GridVmf3Data::Data::NUM_DATA_ELEMENTS; i++) {
      ec = std::from_chars(skipws(num), line + std::strlen(line),
                           grid->data(idx)->data()[i]);
      if (!(ec.ec == std::errc{})) ++error;
      num = ec.ptr;
    }
    ++idx;

    /* get next line */
    fin.getline(line, MAX_VMF3_CHARS);
  }

  /* check for errors */
  if (error) {
    fprintf(stderr,
            "[ERROR] Failed collecting data from VMF3 grid file %s (traceback: "
            "%s)\n",
            fn, __func__);
    fprintf(stderr, "[ERROR] C'nued erronuous line is: %s (traceback: %s)\n",
            line, __func__);
    return 5;
  }

  /* check validity */
  if (idx == (std::size_t)grid->num_pts()) {
    return 0;
  } else {
    fprintf(stderr,
            "[ERROR] Failed collecting data from VMF3 grid file %s (traceback: "
            "%s)\n",
            fn, __func__);
    fprintf(stderr, "[ERROR] Unspecified error! (traceback: %s)\n", __func__);
    return 9;
  }
}