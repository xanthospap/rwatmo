#include "vmf3.hpp"


std::vector<const char *>::const_iterator
dso::vmf3_details::find_if_sorted_string(const char *site,
                      const std::vector<const char *> &sites) noexcept {
  for (auto it = sites.begin(); it != sites.end(); ++it) {
    int cmp = std::strcmp(*it, site);
    if (!cmp)
      return it;
    if (cmp > 0)
      return sites.end();
  }
  return sites.end();
}
