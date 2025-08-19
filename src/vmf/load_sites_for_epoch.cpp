#include <filesystem>

#include "geodesy/units.hpp"
#include "vmf.hpp"

int dso::Vmf3SiteHandler::load_sites_for_epoch(const ymd_date &ymd,
                                               int day_hours,
                                               int index) noexcept {
  /* construct the (possible) filename, i.e. V3GR_20230105.H06 */
  char filename[48];
  std::sprintf(filename, "V3GR_%d%0d%0d.H%0d", ymd.yr(), ymd.mn(), ymd.dm(),
               day_hours);
  /* concatenate path and filename */
  std::filesystem::path full = std::filesystem::path(mdata_dir) / filename;

  /* load the whole grid for the specified epoch */
  vmf3::GridVmf3Data grid;
  if (load_vmfgr3_grid_map(full.c_str(), &grid)) {
    fprintf(
        stderr,
        "[ERROR] Failed reading VMF3 GR grid from file %s (traceback: %s)\n",
        full.c_str(), __func__);
    return 1;
  }

  /* go on and extract data for the list of sites of interest */
  for (const auto &site : msites) {
    const double lat = dso::rad2deg(site.mcrd.lat());
    const double lon = dso::rad2deg(site.mcrd.lon());
    /* TODO guard for exceptions */
    if (grid.bilinear_interpolation(lat, lon, site.mdata_t1)) {
      return 1;
    }
  }

  return 0;
}
