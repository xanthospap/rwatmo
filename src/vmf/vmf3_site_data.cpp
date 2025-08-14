#include "vmf.hpp"
#include <charconv>
#include <cstring>

namespace {
const char *skipws(const char *str) noexcept {
  while (*str && *str == ' ')
    ++str;
  return str;
}

/** Resolve a VMF3 site-specific data line; note that the site name is omited.
 *
 * An example of such line is:
 * ADEA      58496.00  0.00120787  0.00065184  2.2275  0.0626   979.44
 * -3.07   2.67  -0.406   0.051   0.112   0.023
 *
 * For an explanation of the fields, see ege the yearly files in
 * https://vmf.geo.tuwien.ac.at/trop_products/DORIS/V3GR/V3GR_OP/
 *
 * Currently, only for DORIS sites, no other technique-specific data files
 * are checked.
 *
 * @param[in] line The data /record line to resolve.
 * @param[out] block A vmf3::SiteBlock holding the resolved data (apart from
 * date).
 * @param[out] t The epoch of the data line.
 * @return Anything other than 0 denotes an error.
 */
int resolve_vmf3gr_line(const char *line, dso::vmf3::SiteBlock &block,
                        dso::MjdEpoch &t) noexcept {
  /**/
  const int sz = std::strlen(line);

  /* first field is the site name */
  std::memcpy(block.site(), line, 4 * sizeof(char));

  double data[dso::vmf3::SiteBlock::NUM_DATA_FIELDS];
  int error = 0;

  /* resolve the date */
  {
    auto res = std::from_chars(skipws(line + 5), line + sz, data[0]);
    if (res.ec != std::errc{})
      ++error;
    double ipart;
    double fpart = std::modf(data[0], &ipart);
    t = dso::MjdEpoch((int)ipart, dso::FractionalSeconds(fpart * 86400e0));
    if (error) {
      fprintf(stderr,
              "[ERROR] Failed parsing date from VMF3 site-specific file "
              "(traceback: %s)\n",
              __func__);
      fprintf(stderr, "[ERROR] ct'nued line was: %s (traceback: %s)\n", line,
              __func__);
      return 1;
    }
  }

  /* resolve numeric values */
  const char *start = line + 19;
  for (int i = 0; i < dso::vmf3::SiteBlock::NUM_DATA_FIELDS; i++) {
    auto res = std::from_chars(skipws(start), line + sz, block.data()[i]);
    if (res.ec != std::errc{})
      ++error;
    start = res.ptr;
  }

  /* an error occured */
  if (error) {
    fprintf(stderr,
            "[ERROR] Failed parsing float field from VMF3 site-specific file "
            "(traceback: %s)\n",
            __func__);
    return error;
  }

  return 0;
}
} /* anonymous namespace */

int dso::Vmf3SiteStream::get_next_block(dso::vmf3::Block &block) noexcept {

  block.reset();
  dso::vmf3::SiteBlock bc;
  dso::MjdEpoch t, tc;
  int error = 0;

  /* already reached EOF */
  if ((mlast_line[0] == 'E' && mlast_line[1] == 'O') &&
      (mlast_line[2] == 'F')) {
    return -1;
  }

  /* first line: it is either buffered or empty (the latter means we are
   * probably reading the first block)
   */
  if (!mlast_line[0]) {
    mstream.getline(mlast_line, MAX_RECORD_CHARS);
  }
  error = resolve_vmf3gr_line(mlast_line, bc, t);
  block.append(bc);

  /* go on reading/parsing line untill we hit a different epoch */
  while (mstream.getline(mlast_line, MAX_RECORD_CHARS) && (!error)) {
    error = resolve_vmf3gr_line(mlast_line, bc, tc);

    if (!error) {
      /* if we encountered a new date, the block is over */
      if (tc != t) {
        error = -10;
      } else {
        /* else append data for new site */
        block.append(bc);
      }
    }
  }

  /* error encounterd! */
  if (error > 0) {
    fprintf(stderr,
            "[ERROR] Failed parsing data block for VMF3 (gr) file %s "
            "(traceback: %s)\n",
            mfilename.c_str(), __func__);
    return 1;
  }

  /* EOF encountered, set accordingly the last buffered line */
  if (mstream.eof()) {
    mlast_line[0] = 'E';
    mlast_line[1] = 'O';
    mlast_line[2] = 'F';
  }

  return 0;
}