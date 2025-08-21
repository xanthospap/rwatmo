#include "geodesy/units.hpp"
#include "vmf3_grid_data.hpp"

int dso::vmf3::GridVmf3Data::bilinear_interpolation(
    double lat_deg, double lon_deg,
    dso::vmf3::GridVmf3Data::Data &out) const noexcept {
  /* wrap longitude in range [0, 360] if needed */
  lon_deg = dso::norm_angle<dso::detail::AngleUnit::Degrees>(lon_deg);

  dso::TickAxis2D_NP::Cell cell;
  try {
    /* may throw */
    cell = axis_.cell(lat_deg, lon_deg);
  } catch (std::exception &e) {
    fprintf(stderr,
            "[ERROR] Exception encountered while specifiying cell for point "
            "lat=%.2f, lon=%.2f (traceback: %s)\n",
            lat_deg, lon_deg);
    fprintf(stderr, "[ERROR] error message is %s (traceback: %s)\n", e.what());
    return 1;
  }

  auto [x1, y1] = axis_.val(cell.bl);
  auto [x2, y2] = axis_.val(cell.tr);
  double x = lat_deg;
  double y = lon_deg;

  // Unwrap the inner (periodic) segment so y is numerically between y1 and y2
  axis_.tick_axis_inner().unwrap_segment_around(y, y1, y2);

  const double dx = (x2 - x1);
  const double dy = (y2 - y1);
  const double inv_dx = 1.0 / dx;
  const double inv_dy = 1.0 / dy;
  const double wx2 = (x2 - x) * inv_dx;  // weight for x1-side
  const double wx1 = (x - x1) * inv_dx;  // weight for x2-side
  const double wy2 = (y2 - y) * inv_dy;  // weight for y1-side
  const double wy1 = (y - y1) * inv_dy;  // weight for y2-side
  for (int k = 0; k < dso::vmf3::GridVmf3Data::Data::NUM_ELEMENTS; k++) {
    const double f11 = this->data(cell.bl)->data()[k];  // (x1,y1)
    const double f21 = this->data(cell.br)->data()[k];  // (x2,y1)
    const double f12 = this->data(cell.tl)->data()[k];  // (x1,y2)
    const double f22 = this->data(cell.tr)->data()[k];  // (x2,y2)
    const double fxy1 = wx2 * f11 + wx1 * f21;          // along x at y1
    const double fxy2 = wx2 * f12 + wx1 * f22;          // along x at y2
    out.data()[k] = wy2 * fxy1 + wy1 * fxy2;            // along y
  }

  return 0;
}
int dso::vmf3::GridVmf3Data::bilinear_interpolation_nocheck(
    double lat_deg, double lon_deg,
    dso::vmf3::GridVmf3Data::Data &out) const noexcept {
  /* Unchecked cell: may produce out-of-range indices if the caller didn’t
   * validate
   */
  lon_deg = dso::norm_angle<dso::detail::AngleUnit::Degrees>(lon_deg);
  dso::TickAxis2D_NP::Cell cell = axis_.cell_nocheck(lat_deg, lon_deg);

  auto [x1, y1] = axis_.val(cell.bl);
  auto [x2, y2] = axis_.val(cell.tr);
  double x = lat_deg, y = lon_deg;
  axis_.tick_axis_inner().unwrap_segment_around(y, y1, y2);

  const double dx = (x2 - x1), dy = (y2 - y1);
  const double inv_dx = 1.0 / dx, inv_dy = 1.0 / dy;
  const double wx2 = (x2 - x) * inv_dx, wx1 = (x - x1) * inv_dx;
  const double wy2 = (y2 - y) * inv_dy, wy1 = (y - y1) * inv_dy;

  for (int k = 0; k < dso::vmf3::GridVmf3Data::Data::NUM_ELEMENTS; ++k) {
    const double f11 = this->data(cell.bl)->data()[k];
    const double f21 = this->data(cell.br)->data()[k];
    const double f12 = this->data(cell.tl)->data()[k];
    const double f22 = this->data(cell.tr)->data()[k];
    const double fxy1 = wx2 * f11 + wx1 * f21;
    const double fxy2 = wx2 * f12 + wx1 * f22;
    out.data()[k] = wy2 * fxy1 + wy1 * fxy2;
  }
  return 0;
}