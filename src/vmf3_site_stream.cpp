#include "vmf3.hpp"
#include <algorithm>
#include <cstdlib>
#include <cstring>

dso::Vmf3SiteStream::Vmf3SiteStream(const char *fn,
                                    const std::vector<const char *> &sites)
    : mstream(fn) {
  initialize(sites);
}

std::vector<const char *>::iterator
dso::Vmf3SiteStream::append_site(const char *site) noexcept {
  /* current size (in chars) */
  int psize = msites.size() * SITE_LEN;

  /* allocate (new)memory */
  auto ptr = new char[psize + SITE_LEN];

  /* copy sites to memmory (in lexicographical order) */
  int sitenr = 0;
  for (const auto &s : msites) {
    std::memcpy(ptr + sitenr * SITE_LEN, s, SITE_LEN);
    ++sitenr;
  }

  /* free memory */
  delete[] mmemsites;
  mmemsites = ptr;

  /* pointers are now invalidated; reset */
  for (int i = 0; i < msites.size(); i++)
    msites[i] = mmemsites + i * SITE_LEN;

  /* copy new site to internal memory */
  std::memcpy(mmemsites + psize, site, MAX_SITE_CHARS);
  mmemsites[psize + MAX_SITE_CHARS] = '\0';

  /* add new site */
  auto it = msites.insert(std::lower_bound(msites.begin(), msites.end(), site,
                                           [](const char *a, const char *b) {
                                             return std::strcmp(a, b) < 0;
                                           }),
                          mmemsites + psize);

  /* data vector right shift at index */
  int j = std::distance(msites.begin(), it) * RECS_PER_SITE;
  auto dit =
      mdata.insert(mdata.begin() + j, RECS_PER_SITE, dso::Vmf3SiteData{});

  /* set dates at newly added data elements */

  for (int i = 0; i < RECS_PER_SITE - 1; i++) {
    ++dit;
    dit->t() = dso::MjdEpoch::min();
  }

  return it;
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
    } else if (t.diff<dso::DateTimeDifferenceType::FractionalDays>(ct).days() >
               .74) {
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
  /* quick return, we already have the data for this epoch */
  if (t >= interval_start() && t < interval_stop())
    return 0;

  /* maybe its the first call to the function */
  if ((interval_start() == MjdEpoch::max()) &&
      (interval_stop() == MjdEpoch::min())) {
    int error = 0;
    error += mstream.parse_block(msites, mdata, RECS_PER_SITE, 0);
    if (!error)
      t0 = mstream.current_epoch();
    error += mstream.parse_block(msites, mdata, RECS_PER_SITE, 1);
    if (!error)
      t1 = mstream.current_epoch();
    if (error) {
      fprintf(stderr,
              "[ERROR] Failed getting initial epochs from VMF3 stream %s "
              "(traceback: %s)\n",
              mstream.fn(), __func__);
      return error;
    }
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

// int dso::Vmf3SiteStream::operator()(const dso::MjdEpoch &t) noexcept {
//   /* locate interval and buffer data */
//   if (this->set_at_epoch(t)) {
//     fprintf(stderr,
//             "[ERROR] Failed locating suitable interval for epoch %.6f for
//             VMF3 " "stream %s (traceback: %s)\n", t.imjd() +
//             t.fractional_days(), mstream.fn(), __func__);
//     return 1;
//   }
//
//   return 0;
// }

int dso::Vmf3SiteStream::site_vmf3(const char *site, const dso::MjdEpoch &t,
                                   dso::Vmf3SiteData &vmf3,
                                   int *site_index) noexcept {
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
  const auto dt10 =
      t1.diff<dso::DateTimeDifferenceType::FractionalDays>(t0).days();
  const auto dt0 =
      t.diff<dso::DateTimeDifferenceType::FractionalDays>(t0).days();
  const auto dt1 = t1.diff<dso::DateTimeDifferenceType::FractionalDays>(t);

  const double *y0 = mdata[index].data();
  const double *y1 = mdata[index + 1].data();
  double *__restrict__ y = vmf3.data();
  for (int i = 0; i < 7; i++) {
    // y[i] = (y0[i] * dt1 + y1[i] * dt0) / dt10;
    y[i] = y0[i] + dt0 * (y1[i] - y0[i]) / dt10;
  }

  /* assign epoch */
  vmf3.t() = t;

  /* assign site index */
  if (site_index)
    *site_index = std::distance(msites.cbegin(), it);

  return 0;
}

int dso::Vmf3SiteStream::site_vmf3(
    const char *site, const dso::MjdEpoch &t, double el,
    dso::Vmf3SiteStream::Vmf3Result &res) noexcept {
  /* interpolate VMF3 data */
  dso::Vmf3SiteData site_data;
  int site_index;
  if (site_vmf3(site, t, site_data, &site_index)) {
    fprintf(stderr,
            "[ERROR] Failed interpolating VMF3 data for site %s at %.6f "
            "(traceback: %s)\n",
            t.imjd() + t.fractional_days().days(), __func__);
    return 1;
  }

  /* the FullCoeffs should be at index site_index of the mvfc vector; compute
   * empirical coeffs for epoch
   */
  const auto ec = mvfc[site_index].computeCoeffs(t);

  /* VMF3 wet and dry maping functions for given elevation angle */
  double mfh, mfw;
  ec.mf(el, site_data.ah(), site_data.aw(), res.mfh, res.mfw);

  /* assign zenith path delays */
  res.zhd = site_data.zhd();
  res.zwd = site_data.zwd();

  /* all done */
  return 0;
}

std::vector<dso::vmf3_details::Vmf3FullCoeffs>::iterator
dso::Vmf3SiteStream::set_site_coordinates(
    const char *site, const dso::GeodeticCrdConstView &crd) noexcept {
  /* find site in list of sites */
  const auto it = vmf3_details::find_if_sorted_string(site, msites);
  if (it == msites.end()) {
    fprintf(stderr,
            "[ERROR] Failed to find site %s in stream's site list! (traceback: "
            "%s)\n",
            site, __func__);
    return mvfc.end();
  }

  /* compute vmf3FullCoeffSet for given site */
  int index = std::distance(msites.cbegin(), it);
  if (index >= mvfc.size()) {
    fprintf(stderr,
            "[ERROR] Incosistent size of Vmf3FullCoeffs in Vmf3SiteStream; "
            "site is %s (traceback: %s)\n",
            site, __func__);
    return mvfc.end();
  }

  mvfc[index] = dso::vmf3_details::vmf3FullCoeffSet(crd.lon(), crd.lat());
  return mvfc.begin() + index;
}
