#include <cmath>
#include <stdexcept>

#include "vmf3_grid_stream.hpp"

dso::vmf3::TickAxis2D::index_type dso::vmf3::TickAxis2D::index(
    double val_out, double val_in) const noexcept {
  /* row ("bottom") */
  int i = outter_.index(val_out);
  /* col ("left") */
  int j = inner_.index(val_in);
  /* Linearize in signed 64-bit to avoid silent wrap on negatives. */
  return static_cast<std::int64_t>(i) *
             static_cast<std::int64_t>(inner_.num_pts()) +
         static_cast<std::int64_t>(j);
}

dso::vmf3::TickAxis2D::Cell dso::vmf3::TickAxis2D::cell(double val_out,
                                                        double val_in) const {
  const int nrows = outter_.num_pts();
  const int ncols = inner_.num_pts();
  const int i = outter_.index(val_out);
  const int j = inner_.index(val_in);

  /* Step directions (+1 for ascending, -1 for descending) */
  const int so = (outter_.tick_step() > 0) ? +1 : -1;
  const int si = (inner_.tick_step() > 0) ? +1 : -1;

  /* Validate they can form a cell (must allow i+1 and j+1) */
  const int io = i + so;
  const int ji = j + si;
  if (i < 0 || i >= nrows || io < 0 || io >= nrows || j < 0 || j >= ncols ||
      ji < 0 || ji >= ncols) {
    throw std::out_of_range(
        "[ERROR] no surrounding cell (indices out of range) for (" +
        std::to_string(val_out) + ", " + std::to_string(val_in) +
        ") (traceback: vmf3::TickAxis2D::cell)\n");
  }

  /* Validate the point actually lies between ticks i..i+1 and j..j+1 */
  if (!outter_.within_segment_inclusive(val_out, i, i + so) ||
      !inner_.within_segment_inclusive(val_in, j, j + si)) {
    throw std::out_of_range(
        "[ERROR] point not enclosed by adjacent ticks for (" +
        std::to_string(val_out) + ", " + std::to_string(val_in) +
        ") (traceback: vmf3::TickAxis2D::cell)\n");
  }

  /* Linearize (row-major; inner_ is the fast axis) */
  const int stride = ncols;
  const int bl = i * stride + j;    // bottom-left
  const int br = i * stride + ji;   // move along inner axis by si
  const int tl = io * stride + j;   // move along outer axis by so
  const int tr = io * stride + ji;  // both
  return Cell{bl, br, tl, tr};
}

dso::vmf3::TickAxis2D::Cell dso::vmf3::TickAxis2D::cell_nocheck(
    double val_out, double val_in) const noexcept {
  const int i = outter_.index(val_out);  // bottom on outer axis
  const int j = inner_.index(val_in);    // left   on inner axis

  const int so = (outter_.tick_step() > 0) ? +1 : -1;
  const int si = (inner_.tick_step() > 0) ? +1 : -1;

  const int stride = inner_.num_pts();
  const int bl = i * stride + j;
  const int br = i * stride + (j + si);
  const int tl = (i + so) * stride + j;
  const int tr = (i + so) * stride + (j + si);
  return Cell{bl, br, tl, tr};
}