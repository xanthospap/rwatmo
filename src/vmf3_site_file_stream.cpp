#include "vmf3.hpp"
#include <algorithm>
#include <charconv>
#include <cstring>
#include <stdexcept>
#include <vector>

namespace {
const char *skipws(const char *str) noexcept {
  while (*str && *str == ' ')
    ++str;
  return str;
}

/** Resolve a VMF3 site-specific data line; note that the site name is omited.
 *
 * An example of such line is:
 * ADEA      54832.00  0.00120058  0.00059661  2.2546  0.0323   991.72
 * -1.95   2.73 For an explanation of the fields, see ege the yearly files in
 * https://vmf.geo.tuwien.ac.at/trop_products/DORIS/VMF3/VMF3_OP/yearly/
 *
 * Currently, only for DORIS sites, no other technique-specific data files
 * are checked.
 *
 * @param[in] line The data /record line to resolve.
 * @param[out] out An Vmf3SiteData holding the resolved data, apart from site
 *                 name.
 * @return Anything other than 0 denotes an error.
 */
int resolve_vmf3_line(const char *line, dso::Vmf3SiteData &out) noexcept {
  const char *start = skipws(line + 10);
  int sz = std::strlen(line);
  double data[10];
  int error = 0;

  /* resolve numeric values */
  for (int i = 0; i < 8; i++) {
    auto res = std::from_chars(start, line + sz, data[i]);
    if (res.ec != std::errc{})
      ++error;
    start = skipws(res.ptr);
  }

  /* an error occured */
  if (error)
    return error;

  /* assign data to instance */
  out.ah() = data[1];
  out.aw() = data[2];
  out.zhd() = data[3];
  out.zwd() = data[4];
  out.pressure() = data[5];
  out.temperature() = data[6];
  out.water_vapour_pressure() = data[7];

  /* assign epoch to instance */
  {
    double ipart;
    double fpart = std::modf(data[0], &ipart);
    out.t() =
        dso::MjdEpoch((int)ipart, dso::FractionalSeconds(fpart * 86400e0));
  }

  return 0;
}

/** Resolve the date from a VMF3 site-specific data line.
 *
 * An example of such line is:
 * ADEA      54832.00  0.00120058  0.00059661  2.2546  0.0323   991.72
 * -1.95   2.73 For an explanation of the fields, see ege the yearly files in
 * https://vmf.geo.tuwien.ac.at/trop_products/DORIS/VMF3/VMF3_OP/yearly/
 *
 * Currently, only for DORIS sites, no other technique-specific data files
 * are checked.
 *
 * @param[in] line The data /record line to resolve.
 * @param[out] t   A MjdEpoch instance holding the resolved datetime
 * @return Anything other than 0 denotes an error.
 */
int resolve_date(const char *line, dso::MjdEpoch &t) noexcept {
  int sz = std::strlen(line);
  const char *start = skipws(line + 10);
  int error = 0;
  double data;
  /* resolve date (given in MJD) */
  auto res = std::from_chars(start, line + sz, data);
  if (res.ec != std::errc{})
    ++error;
  double ipart;
  double fpart = std::modf(data, &ipart);
  t = dso::MjdEpoch((int)ipart, dso::FractionalSeconds(fpart * 86400e0));
  return error;
}

int skip_comment_lines(std::ifstream &fin, char *line, int MAX_CHARS) noexcept {
  while ((*line && *line == '#') && fin.good())
    fin.getline(line, MAX_CHARS);
  return 0;
}

} /* unnamed namespace */

//dso::MjdEpoch dso::Vmf3SiteFileStream::buffered_epoch() const {
//  if (bline[0] == '#' && bline[1] == '\0')
//    return dso::MjdEpoch::max();
//  dso::MjdEpoch t;
//  if (resolve_date(bline, t)) {
//    throw std::runtime_error("[ERROR]\n");
//  }
//  return t;
//}

/* Skip the current block from the stream file.
 *
 * Starting from the already buffered line, keep on reading new lines from the
 * stream until we encounter a date that is later than the current one.
 *
 * At exit, the instance's bline will hold the last line read, i.e. the first
 * line (record) of the new block.
 *
 * @return Anything other than 0 denotes an error.
 */
int dso::Vmf3SiteFileStream::skip_block() noexcept {
  //printf("Skipping block ... %s\n", __func__);
  /* ti: epoch of current line in buffer 
   * t : epoch of block to be skipped
   */
  int error = 0;
  dso::MjdEpoch ti = mnext;
  dso::MjdEpoch t = ti;
  
  /* loop through stream lines until we meet a later date */
  while (mstream.getline(bline, LINE_SZ) && (!error) && (ti == t)) {
    if ((*bline && *bline != '#'))
      error += resolve_date(bline, ti);
  }

  /* something is wrong, report it */
  if (error || (!mstream.good())) {
    if (error) 
      fprintf(stderr, "[ERROR] Failed resolving date from line %s from VMF3 file %s (traceback: %s)\n", bline, fn(), __func__);
    else
      fprintf(stderr, "[ERROR] Failed reading from VMF3 file %s; last line was: %s (traceback: %s)\n", fn(), bline, __func__);
  }

  /* assign epochs */
  mcurrent = t;
  mnext = ti;

  return error;
}

/* Read through a given block (of same epoch) and store VMF3 data.
 *
 * Starting from the already buffered line, keep on reading and parsing
 * new lines from the stream untill we encounter a date that is later than
 * the given one.
 * If the line holds data for a site of interest, store VMF3 records in the
 * passed-in data vector.
 * At exit, the instance's bline will hold the last line read.
 *
 * @param[in] site A vector of sites (for DORIS 4-chars); each site name
 *            should be a null-terminated C-string
 * @param[out] data The vector where the data for the sites of interest will
 *            be stored. Data (i.e. VMF3 records) are copyied into the vector,
 *            NOT APPENDED, hence the vector should have the right size at
 *            input.
 *            If a given record line matches the site name at index j of the
 *            sites (input) vector, then the corresponding VMF3 (parsed) data
 *            will be stored at index j*recs_per_size+rec_offset.
 *            This "indexing" is used to allow multiple records of different
 *            datetimes to be stored in the data vector.
 * @param[in] recs_per_site Number of records per site in the data vector (see
 *            above).
 * @param[in] rec_offset Offset to be used for indexing within the data vector
 *            (see above).
 * @return Anything other than 0 denotes an error.
 */
int dso::Vmf3SiteFileStream::parse_block(const std::vector<const char *> &sites,
                                         std::vector<Vmf3SiteData> &data,
                                         int recs_per_site,
                                         int rec_offset) noexcept {
  // printf("Parsing block ... %s\n", __func__);
  /* get current time from buffered line (skip comments) */
  //skip_comment_lines(mstream, bline, LINE_SZ);
  //dso::MjdEpoch t;
  //if (resolve_date(bline, t)) {
  //  fprintf(stderr,
  //          "[ERROR] Failed resolving date from VMF3 line: %s; file was %s "
  //          "(traceback: %s)\n",
  //          bline, mfn, __func__);
  //  return 1;
  //}


  /* assign current epoch */
  mcurrent = mnext;

  /* loop through stream lines untill we meet a later date */
  int error = 0;
  dso::MjdEpoch ti = mcurrent;
  dso::Vmf3SiteData d;
  char site[5];
  site[4] = '\0';
  while ((!error) && (ti == mcurrent)) {
    if (*bline && *bline != '#') {
      /* copy site name so that we can compare it */
      std::memcpy(site, bline, 4);
      /* are we interested in the site ? */
      auto it = vmf3_details::find_if_sorted_string(site, sites);
      if (it != sites.end()) {
        /* if yes, resolve the line and place the data in the data vector */
        if (resolve_vmf3_line(bline, d)) {
          fprintf(stderr,
                  "[ERROR] Failed resolving VMF3 line: %s (traceback: %s)\n",
                  bline, __func__);
          ++error;
        }
        auto data_index =
            std::distance(sites.begin(), it) * recs_per_site + rec_offset;
        data[data_index] = d;
        ti = d.t();
      } else {
        error += resolve_date(bline, ti);
      }
    } 

    /* buffer next line to bline */
    mstream.getline(bline, LINE_SZ);
    error += (!mstream.good());
  }
  
  /* something is wrong, report it */
  if (error || (!mstream.good())) {
    if (error) 
      fprintf(stderr, "[ERROR] Failed resolving line % from VMF3 file %s (traceback: %s)\n", bline, fn(), __func__);
    else
      fprintf(stderr, "[ERROR] Failed reading from VMF3 file %s; lat line was: %s (traceback: %s)\n", fn(), bline, __func__);
  }

  /* assign next (buffered) epoch */
  mnext = ti;

  return error;
}
