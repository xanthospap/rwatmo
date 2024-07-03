#include "vmf3.hpp"
#include "datetime/calendar.hpp"
#include <algorithm>
#include <vector>

int dso::Vmf3SiteStream::collect_epoch(const dso::MjdEpoch &t) {
}

const char *skipws(const char *str) noexcept {
  while (*str && *str == ' ') ++str;
  return str;
}

/** Resolve a VMF3 site-specific data line; note that the site name is omited.
 *
 * An example of such line is:
 * ADEA      54832.00  0.00120058  0.00059661  2.2546  0.0323   991.72  -1.95   2.73
 * For an explanation of the fields, see ege the yearly files in 
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
  for (int i=0; i<8; i++) {
    auto res = std::from_chars(start, line+sz, data[i]);
    if (res.ec != std::errc{}) ++error;
    start = skipws(res.ptr);
  }

  /* an error occured */
  if (error) return error;

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
    out.t() = dso::MjdEpoch((int)ipart, dso::FractionalSeconds(fpart * 86400e0));
  }

  return 0;
}

/** Resolve the date from a VMF3 site-specific data line.
 *
 * An example of such line is:
 * ADEA      54832.00  0.00120058  0.00059661  2.2546  0.0323   991.72  -1.95   2.73
 * For an explanation of the fields, see ege the yearly files in 
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
  /* get current time from buffered line */
  dso::MjdEpoch t;
  if (resolve_date(bline, t)) {
    fprintf(stderr,
            "[ERROR] Failed resolving date from VMF3 line: %s; file was %s "
            "(traceback: %s)\n",
            bline, mfn, __func__);
    return 1;
  }

  /* loop through stream lines until we meet a later date */
  int error = 0;
  dso::MjdEpoch ti(t);
  while (mstream.getline(bline, Vmf3FileStream::LINE_SZ) && (!error) &&
         (ti == t)) {
    error += resolve_date(bline, ti);
  }

  return error;
}

auto find_if_sorted_string(const char *site,
                           const std::vector<const char *> &sites) noexcept {
  for (auto it = sites.begin(); it != sites.end(); i++) {
    int cmp = std::strcmp(*it, site);
    if (!cmp)
      return it;
    if (cmp > 0)
      return sites.end();
  }
  return sites.end();
}

/* 
 * Starting from the already buffered line, keep on reading and parsing 
 * new lines from the stream untill we encounter a date that is later than 
 * the given one.
 *
 * At exit, the instance's bline will hold the last line read.
 */
int dso::Vmf3SiteFileStream::parse_block(const std::vector<const char *> &sites) noexcept {
  /* get current time from buffered line */
  dso::MjdEpoch t;
  if (resolve_date(bline, t)) {
    fprintf(stderr,
            "[ERROR] Failed resolving date from VMF3 line: %s; file was %s "
            "(traceback: %s)\n",
            bline, mfn, __func__);
    return 1;
  }

  /* loop through stream lines untill we meet a later date */
  int error = 0;
  dso::MjdEpoch ti(t);
  dso::Vmf3Data d;
  char site[5];
  site[4] = '\0';
  while ((!error) && (ti == t)) {
    /* copy site name so that we can compare it */
    std::memcpy(site, bline, 4);
    /* are we interested in the site ? */
    auto it = std::find_if_sorted_string(site, msites);
    if (it != msites.end())
      
    mstream.getline(bline, Vmf3FileStream::LINE_SZ);
  }

  //int error = 0;
  //dso::Vmf3Data out;
  //dso::MjdEpoch ti;
  //while (!error) {
  //  error += resolve_vmf3_line(bline, out);
  //  error += (!mstream.is_good());
  //}

  return error;
}

int dso::Vmf3SiteStream::forward_search(const dso::MjdEpoch &t, char *bline) {
  /* first set second epoch to first epoch */
  swap_epochs();
  
  /* read blocks until we are placed in the suitable interval. Note that the 
   * first line to be read, is already buffered! 
   */
  char *line = mstream.bline;
  int error = 0;
  while (!error) {
    /* skip header/comments */
    while (line[0] == '#')
      mstream.getline(line, dso::Vmf3FileStream::LINE_SZ);
    dso::Vmf3Data d;
    if (resolve_vmf3_line(line, d)) {
      fprintf(stderr,
              "[ERROR] Failed resolving VMF3 data line: %s (traceback: %s)\n",
              line, __func__);
      return -1;
    }
    if (t >= mdata[0].t() && t < d.t()) {
      /* resolve block */
      return 0;
    } else if (t.diff<dso::DateTimeDifferenceType::FractionalDays>(d.t()) >
               1e0) {
      /* skip block */
    } else {
      /* resolve block */
      swap_epochs();
    }

    mstream.getline(line, dso::Vmf3FileStream::LINE_SZ);
    error += !(mstream.is_good());
  }
}

int dso::Vmf3SiteStream::set_at_epoch(const dso::MjdEpoch &t) {
  /* quick return, we already have the data for this epoch */
  if (mdata[0].t() >= t && mdata[1].t() < t) return 0;
  
  /* need to go further in file */
  if (t>= mdata[1].t()) {
    char line[SZ];
    mstream.getline(line, SZ);
    /* skip header/comments */
    while (line[0] == '#' ) mstream.getline(line, SZ);
    dso::Vmf3Data d;
    if (resolve_vmf3_line(line, d)) {
      fprintf(stderr, "[ERROR] Failed resolving VMF3 data line: %s (traceback: %s)\n", line, __func__);
      return -1;
    }
    while (d.t()) {
    }
  }

}
