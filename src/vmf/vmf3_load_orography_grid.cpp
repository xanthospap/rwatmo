#include <charconv>
#include <cstdio>
#include <cstring>
#include <fstream>

#include "vmf3.hpp"

namespace {
constexpr const int MAX_VMF3_CHARS = 24;
const char *skipws(const char *line) noexcept {
  while (*line && *line == ' ') ++line;
  return line;
}
}  // namespace

int dso::Vmf3SiteHandler::load_sites_orography(
    const char *fn, std::size_t num_vals,
    dso::vmf3::GridVmf3Data *grid) noexcept {
  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(
        stderr,
        "[ERROR] Failed opening VMF3 orography data file %s (traceback: %s)\n",
        fn, __func__);
    return 1;
  }

  char line[MAX_VMF3_CHARS];
  int error = 0;

  std::size_t idx = 0;
  while (!(error) &&
         ((idx < num_vals) && (fin.getline(line, MAX_VMF3_CHARS)))) {
    auto ec = std::from_chars(skipws(line), line + std::strlen(line),
                              grid->data(idx)->oro_ell());
    if (!(ec.ec == std::errc{})) ++error;
  }

  return 0;
}
