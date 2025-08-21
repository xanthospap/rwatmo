#include <filesystem>

#include "geodesy/units.hpp"
#include "vmf3.hpp"

int dso::Vmf3SiteHandler::load_sites_for_epoch(
    const ymd_date &ymd, int day_hours, bool load_orography_ell) noexcept {
  /* construct the (possible) filename, i.e. V3GR_20230105.H06 */
  char filename[48];
  std::sprintf(filename, "V3GR_%d%02d%02d.H%02d", ymd.yr().as_underlying_type(),
               ymd.mn().as_underlying_type(), ymd.dm().as_underlying_type(),
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

  /* if demanded, also load the orography grid file */
  if (load_orography_ell) {
    int resolution = grid.axis().tick_axis_inner().tick_step();
    std::sprintf(filename, "orography_ell_%dx%d", resolution, resolution);
    full = std::filesystem::path(mdata_dir) / filename;
    if (load_sites_orography(full.c_str(), grid.num_pts(), &grid)) {
      fprintf(stderr,
              "[ERROR] Failed reading the orography file %s (traceback: %s)\n",
              full.c_str(), __func__);
      return 2;
    }
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
      /* we will need the surrounding cell */
      const auto cell = grid.axis().cell(lat, lon);

      /* if we 'straddle the seam', we must rearrange lon values of the nodes
       * Example:
       * grid: 0, 2.5, …, 357.5 (periodic, 360°)
       * query: lon = 359°
       * segment picked from indices: (y0, y1) = (357.5, 0)
       * should become => y0 = 357.5, y1 = 360, v = 359
       */
      double lon_in = lon;
      auto [dummy_x1, y1] = grid.axis().val(cell.bl);
      auto [dummy_x2, y2] = grid.axis().val(cell.br);
      grid.axis().tick_axis_inner().unwrap_segment_around(lon_in, y1, y2);
      if (lon_in != lon) {
        fprintf(
            stderr,
            "[ERROR] Longitude for site %s seems to straddle the seam in "
            "longitude, but i cannot compute correct nodes! (traceback: %s)\n",
            site.msite.name(), __func__);
        fprintf(stderr,
                "[ERROR] Longitude that caused the problem was %.9f[deg] "
                "(traecback: %s)\n",
                lon, __func__);
        return 1;
      }

      /* load data for each cell, in the order: [bl, br, tl, tr] */
      vmf3::GridVmf3Data::Data *data = &(site.mdata_t1[0]);
      data->lat_deg() = grid.data(cell.bl)->lat_deg();
      data->lon_deg() = y1;  // grid.data(cell.bl)->lon_deg();
      data->oro_ell() = grid.data(cell.bl)->oro_ell();
      std::memcpy(data->data(), grid.data(cell.bl)->data(),
                  sizeof(double) * vmf3::GridVmf3Data::Data::NUM_DATA_ELEMENTS);
      /* bottom right */
      data = &(site.mdata_t1[1]);
      data->lat_deg() = grid.data(cell.br)->lat_deg();
      data->lon_deg() = y2;  // grid.data(cell.br)->lon_deg();
      data->oro_ell() = grid.data(cell.br)->oro_ell();
      std::memcpy(data->data(), grid.data(cell.br)->data(),
                  sizeof(double) * vmf3::GridVmf3Data::Data::NUM_DATA_ELEMENTS);
      /* top left */
      data = &(site.mdata_t1[2]);
      data->lat_deg() = grid.data(cell.tl)->lat_deg();
      data->lon_deg() = y1;  // grid.data(cell.tl)->lon_deg();
      data->oro_ell() = grid.data(cell.tl)->oro_ell();
      std::memcpy(data->data(), grid.data(cell.tl)->data(),
                  sizeof(double) * vmf3::GridVmf3Data::Data::NUM_DATA_ELEMENTS);
      /* top right */
      data = &(site.mdata_t1[3]);
      data->lat_deg() = grid.data(cell.tr)->lat_deg();
      data->lon_deg() = y2;  // grid.data(cell.tr)->lon_deg();
      data->oro_ell() = grid.data(cell.tr)->oro_ell();
      std::memcpy(data->data(), grid.data(cell.tr)->data(),
                  sizeof(double) * vmf3::GridVmf3Data::Data::NUM_DATA_ELEMENTS);

      /* just to be sure, verify site lays within the cell */
      if (!(grid.data(cell.bl)->lat_deg() <= lat &&
            grid.data(cell.bl)->lon_deg() <= lon &&
            grid.data(cell.br)->lat_deg() <= lat &&
            grid.data(cell.br)->lon_deg() > lon &&
            grid.data(cell.tl)->lat_deg() > lat &&
            grid.data(cell.tl)->lon_deg() <= lon &&
            grid.data(cell.tr)->lat_deg() > lat &&
            grid.data(cell.tr)->lon_deg() > lon)) {
        fprintf(stderr,
                "[ERROR] Unexpected error encountered while retrieving data "
                "from sourounding cell for site %s (tracebac: %s)\n",
                site.msite.name(), __func__);
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

  /* assign date */
  mt1 = MjdEpoch(ymd.yr(), ymd.mn(), ymd.dm(),
                 dso::FractionalSeconds(day_hours * 3600e0));

  return 0;
}
