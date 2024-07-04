#include "vmf3.hpp"
#include <algorithm>
#include <cstring>

dso::Vmf3SiteStream::Vmf3SiteStream(const char *fn,
                                    const std::vector<const char *> &sites)
    : mstream(fn) {
  set_local_vectors(sites);
}

void dso::Vmf3SiteStream::set_local_vectors(
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

  return;
}

int dso::Vmf3SiteStream::forward_search(const dso::MjdEpoch &t) noexcept {
  /* first set second epoch to first epoch */
  swap_epochs();

  /* read blocks until we are placed in the suitable interval. Note that the
   * first line to be read, is already buffered!
   */
  int error = 0;
  while (!error) {
    dso::MjdEpoch ct = mstream.buffered_epoch();
    if (t >= mdata[0].t() && t < ct) {
      /* resolve block, place at index 1 */
      return mstream.parse_block(msites, mdata, RECS_PER_SITE, 1);
    } else if (t.diff<dso::DateTimeDifferenceType::FractionalDays>(ct) >
               .74) {
      /* skip block alltogether, more than (3/4)day in past */
      mstream.skip_block();
    } else {
      /* resolve block */
      error += mstream.parse_block(msites, mdata, RECS_PER_SITE, 1);
      swap_epochs();
    }
    error += !(mstream.good());
  }

  return error;
}

int dso::Vmf3SiteStream::set_at_epoch(const dso::MjdEpoch &t) noexcept {
  /* quick return, we already have the data for this epoch */
  if (mdata[0].t() >= t && mdata[1].t() < t) return 0;

  /* maybe the epoch is out of range in past */
  if (t < mdata[0].t()) {
    fprintf(stderr,
            "[ERROR] Failed searching for VMF3 data for given epoch! Cannot go "
            "backwards in time! (traceback: %s)\n",
            __func__);
    return 1;
  }

  /* need to go further in file */
  return forward_search(t);
}
