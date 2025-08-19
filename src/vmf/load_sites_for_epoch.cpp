#include <filesystem>

#include "geodesy/units.hpp"
#include "vmf.hpp"

int dso::Vmf3SiteHandler::load_sites_for_epoch(const ymd_date &ymd,
                                               int day_hours) noexcept {
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

  /* grid limits (longitude is periodic!) */
  const double lat_min = std::min(grid.axis().tick_axis_outter().tick_start(),
                                  grid.axis().tick_axis_outter().tick_stop());
  const double lat_max = std::max(grid.axis().tick_axis_outter().tick_start(),
                                  grid.axis().tick_axis_outter().tick_stop());

  /* go on and extract data for the list of sites of interest */
  for (auto &site : msites) {
    const double lat = dso::rad2deg(site.mcrd.lat());
    /* longitudes in VMF3 grids are in the range [0-360] not [-180, 180] */
    const double lon = dso::norm_angle<dso::detail::AngleUnit::Degrees>(
        dso::rad2deg(site.mcrd.lon()));
    /* make sure point is within grid, so that we can perform bilinear in. */
    if (lat >= lat_min && lat < lat_max) {
      /* preform bilinear interpolation; store data in site's mdata_t1 */
      if (grid.bilinear_interpolation_nocheck(lat, lon, site.mdata_t1)) {
        return 1;
      }
    } else {
      fprintf(stderr,
              "[ERROR] Given site with lat=%.2f, lon=%.2f falls outside VMF3 "
              "grid; cannot use bilinear interpolation! (traceback: %s)\n",
              lat, lon, __func__);
      return 1;
    }
  }

  return 0;
}
