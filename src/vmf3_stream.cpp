#include "vmf3.hpp"
#include <datetime/calendar.hpp>
#include <datetime/dtdatetime.hpp>
#include <datetime/dtfund.hpp>

int dso::Vmf3SiteStream::collect_epoch(const dso::MjdEpoch &t) {
}


const char *skipws(const char *str) noexcept {
  while (*str && *str == ' ') ++str;
  return str;
}

int resolve_vmf3_line(const char *line, dso::Vmf3Data &out) noexcept {
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

/* Starting from the already buffered line, keep on reading new lines from the 
 * stream untill we encounter a date that is later than the current one.
 * At exit, the instance's bline will hold the last line read.
 */
int dso::Vmf3FileStream::skip_block() noexcept {
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
  while (mstream.getline(bline, Vmf3FileStream::LINE_SZ) && (!error) &&
         (ti == t)) {
    error += resolve_date(bline, ti);
  }

  return error;
}

/* Starting from the already buffered line, keep on reading and parsing 
 * new lines from the stream untill we encounter a date that is later than the 
 * given one.
 * At exit, the instance's bline will hold the last line read.
 */
int dso::Vmf3SiteStream::parse_block() noexcept {
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
  char buf[5];
  buf[4] = '\0';
  while ((!error) && (ti == t)) {
    /* copy site name so that we can compare it */
    std::memcpy(buf, bline, 4);
    /* are we interested in the site ? */
    auto it = std::find_if(msites.begin(), msites.end(), [&](const char *s){return !std::strcmp(s,buf);});
    if (it != msites.end())
    /* TODO how are we going to add the new site? */
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
  
  /* read blocks untill we are placed in the suitable interval. Note that the 
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
