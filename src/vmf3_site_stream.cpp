#include "vmf3.hpp"
#include <algorithm>
#include <cstring>

dso::Vmf3SiteStream::Vmf3SiteStream(const char *fn,
                                    const std::vector<const char *> &sites)
    : mstream(fn) {
  initialize(sites);
}

void dso::Vmf3SiteStream::initialize(
    const std::vector<const char *> &sites) noexcept {
  /* delete memory */
  if (mmemsites)
    delete[] mmemsites;
  msites.clear();
  mdata.clear();

  /* allocate memory */
  const auto num_sites = sites.size();
  mmemsites = new char[num_sites * SITE_LEN];
  std::memset(mmemsites, '\0', num_sites * SITE_LEN);
  msites.reserve(num_sites);
  mdata.reserve(num_sites * RECS_PER_SITE);

  /* copy sites to new memory and set pointers */
  int k = 0;
  for (const auto it : sites) {
    std::strncpy(mmemsites + k * SITE_LEN, it, MAX_SITE_CHARS);
    msites.push_back(mmemsites + k * SITE_LEN);
    ++k;
  }

  /* order sites lexicographycally */
  std::sort(msites.begin(), msites.end(), [&](const char *a, const char *b) {
    return std::strncmp(a, b, MAX_SITE_CHARS) < 0;
  });

  /* initialize mdata with the correct number of records */
  for (int j = 0; j < sites.size(); j++) {
    mdata.push_back(dso::Vmf3SiteData{});
    for (int i = 1; i < RECS_PER_SITE; i++)
      mdata.push_back(dso::Vmf3SiteData{dso::MjdEpoch::min()});
  }

  /* rewind file stream */
  mstream.rewind();

  /* reset dates */
  t0 = dso::MjdEpoch::max();
  t1 = dso::MjdEpoch::min();

  return;
}

int dso::Vmf3SiteStream::forward_search(const dso::MjdEpoch &t) noexcept {
  //printf("Calling %s ...\n", __func__);
  /* first set second epoch to first epoch; t0 and t1 will also be swaped */
  swap_epochs();

  /* read blocks until we are placed in the suitable interval. Note that the
   * first line to be read, is already buffered!
   */
  int error = 0;
  while (!error) {
    /* epoch of buffered line */
    dso::MjdEpoch ct = mstream.next_epoch();
    if (t >= t0 && t < ct) {
      /* resolve block, place at index 1 */
      error = mstream.parse_block(msites, mdata, RECS_PER_SITE, 1);
      if (!error)
        t1 = mstream.current_epoch();
      return error;
    } else if (t.diff<dso::DateTimeDifferenceType::FractionalDays>(ct) > .74) {
      /* too far away, no need to parse the block */
      error += mstream.skip_block();
    } else {
      /* resolve block */
      error += mstream.parse_block(msites, mdata, RECS_PER_SITE, 1);
      if (!error)
        t1 = mstream.current_epoch();
      swap_epochs();
    }
    error += !(mstream.good());
  }

  return error;
}

int dso::Vmf3SiteStream::set_at_epoch(const dso::MjdEpoch &t) noexcept {
  //printf("Calling %s ...\n", __func__);
  /* quick return, we already have the data for this epoch */
  if (t >= interval_start() && t < interval_stop())
    return 0;

  /* maybe its the first call to the function */
  if ((interval_start() == MjdEpoch::max()) &&
      (interval_stop() == MjdEpoch::min())) {
    int error = 0;
    error += mstream.parse_block(msites, mdata, RECS_PER_SITE, 0);
    if (!error) t0 = mstream.current_epoch();
    error += mstream.parse_block(msites, mdata, RECS_PER_SITE, 1);
    if (!error) t1 = mstream.current_epoch();
    if (error) {
      fprintf(stderr,
              "[ERROR] Failed getting initial epochs from VMF3 stream %s "
              "(traceback: %s)\n",
              mstream.fn(), __func__);
      return error;
    }
    printf("Read first two epoch of stream, i.e. %.6f, %.6f\n", interval_start().imjd() + interval_start().fractional_days(), interval_stop().imjd() + interval_stop().fractional_days());
  }

  /* maybe the epoch is out of range in past */
  if (t < interval_start()) {
    fprintf(stderr,
            "[ERROR] Failed searching for VMF3 data for given epoch! Cannot go "
            "backwards in time! (traceback: %s)\n",
            __func__);
    return 1;
  }

  /* need to go further in file */
  return forward_search(t);
}

int dso::Vmf3SiteStream::operator()(const dso::MjdEpoch &t) noexcept {
  /* locate interval and buffer data */
  if (this->set_at_epoch(t)) {
    fprintf(stderr,
            "[ERROR] Failed locating suitable interval for epoch %.6f for VMF3 "
            "stream %s (traceback: %s)\n",
            t.imjd() + t.fractional_days(), mstream.fn(), __func__);
    return 1;
  }

  /**/
  //printf("Left buffer: %.6f right buffer: %.6f (value: %.6f)\n",
  //       interval_start().imjd() + interval_start().fractional_days(),
  //       interval_stop().imjd() + interval_stop().fractional_days(),
  //       t.imjd() + t.fractional_days());

  return 0;
}

int dso::Vmf3SiteStream::site_vmf3(const char *site, const dso::MjdEpoch &t,
                                   dso::Vmf3SiteData &vmf3) noexcept {
  /* find the site in the site list */
  const auto it = vmf3_details::find_if_sorted_string(site, msites);
  if (it == msites.end()) {
    fprintf(stderr,
            "[ERROR] Failed to find site %s in stream's site list! (traceback: "
            "%s)\n",
            site, __func__);
    return 1;
  }

  /* iterator to the data for this site */
  int index = std::distance(msites.cbegin(), it) * RECS_PER_SITE;

  /* validate that we have the correct epochs */
  if (!(t >= mdata[index].t() && t < mdata[index + 1].t())) {
    fprintf(
        stderr,
        "[ERROR] Buffered VMF3 data interval does not match! (traceback: %s)\n",
        __func__);
    return 100;
  }

  /* linear interpolation for all data */
  const auto dt10 = t1.diff<dso::DateTimeDifferenceType::FractionalDays>(t0);
  const auto dt0 = t.diff<dso::DateTimeDifferenceType::FractionalDays>(t0);
  const auto dt1 = t1.diff<dso::DateTimeDifferenceType::FractionalDays>(t);

  const double *y0 = mdata[index].data();
  const double *y1 = mdata[index + 1].data();
  double *__restrict__ y = vmf3.data();
  for (int i = 0; i < 7; i++) {
    // y[i] = (y0[i] * dt1 + y1[i] * dt0) / dt10;
    y[i] = y0[i] + dt0 * (y1[i]-y0[i]) / dt10;
  }

  /* assign epoch */
  vmf3.t() = t;

  return 0;
}
