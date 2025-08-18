#ifndef __DSO_VMF3_TROPO_MODELING_HPP__
#define __DSO_VMF3_TROPO_MODELING_HPP__

#include "vmf3_grid_stream.hpp"

namespace dso {
class Vmf3 {}; /* class Vmf3 */

class Vmf3SiteHandler {
  MjdEpoch mt0;
  MjdEpoch mt1;
  std::vector<vmf3::SiteBlock> msites;
  std::string mdata_dir;

  int load_sites_for_epoch(const MjdEpoch &t);

  auto vmf3(const MjdEpoch &t, double el);
}; /*class Vmf3SiteHandler */
} /* namespace dso  */

#endif