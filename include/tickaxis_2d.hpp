#ifndef __DSO__GENERIC_GRID_2D_CLASS_HPP__
#define __DSO__GENERIC_GRID_2D_CLASS_HPP__

#include "tickaxis_1d.hpp"

namespace dso {
/**
 * @brief 2-D grid made of two TickAxis instances (outer = rows, inner = cols).
 *
 * @tparam POuter  IsPeriodic::NonPeriodic or IsPeriodic::Periodic for OUTER
 * axis
 * @tparam PInner  IsPeriodic::NonPeriodic or IsPeriodic::Periodic for INNER
 * axis
 *
 * Row-major linearization: idx = i * N_inner + j
 * - OUTER axis (“rows”) is the first template parameter.
 * - INNER axis (“cols”) is the second template parameter, contiguous in memory.
 */
template <IsPeriodic POuter, IsPeriodic PInner>
class TickAxis2D {
 public:
  using AxisOut = TickAxis<POuter>;
  using AxisIn = TickAxis<PInner>;
  using index_type = std::int64_t;

  struct Cell {
    index_type bl;  ///< bottom-left  (i,   j)
    index_type br;  ///< bottom-right (i,   j+step_dir_inner)
    index_type tl;  ///< top-left     (i+step_dir_outer, j)
    index_type tr;  ///< top-right    (i+step_dir_outer, j+step_dir_inner)
  };

 private:
  AxisOut outter_;  // rows
  AxisIn inner_;    // cols

 public:
  // ---- ctors ---------------------------------------------------------------

  /**
   * @brief Build from axis ranges.
   *
   * @param ostart, ostop, ostep   Outer axis range & step (inclusive).
   * @param istart, istop, istep   Inner axis range & step (inclusive).
   * @param eps_steps              Divisibility tolerance (fraction of |step|).
   * @param outer_period_hint      Period for OUTER (used only if periodic).
   * @param inner_period_hint      Period for INNER (used only if periodic).
   */
  TickAxis2D(double ostart, double ostop, double ostep, double istart,
             double istop, double istep, double eps_steps = 1e-12,
             double outer_period_hint = 0.0, double inner_period_hint = 0.0)
      : outter_(ostart, ostop, ostep, eps_steps, outer_period_hint),
        inner_(istart, istop, istep, eps_steps, inner_period_hint) {}

  TickAxis2D() = default;

  // ---- basic info ----------------------------------------------------------

  const AxisOut& tick_axis_outter() const noexcept { return outter_; }
  const AxisIn& tick_axis_inner() const noexcept { return inner_; }

  int rows() const noexcept { return outter_.num_pts(); }
  int cols() const noexcept { return inner_.num_pts(); }
  index_type num_pts() const noexcept {
    return static_cast<index_type>(rows()) * static_cast<index_type>(cols());
  }

  // ---- value access --------------------------------------------------------

  /// Map (i,j) -> (outer_value, inner_value). No bounds checks.
  std::pair<double, double> val(int i, int j) const noexcept {
    return {outter_.val(i), inner_.val(j)};
  }

  /// Map linear index -> (outer_value, inner_value). No bounds checks.
  std::pair<double, double> val(index_type idx) const noexcept {
    const index_type ncols = static_cast<index_type>(cols());
    const index_type i = idx / ncols;  // row
    const index_type j = idx % ncols;  // col (may be negative if idx<0)
    return {outter_.val(static_cast<int>(i)), inner_.val(static_cast<int>(j))};
  }

  // ---- indexing ------------------------------------------------------------

  /// Unchecked 2-D indices (outer, inner) for bottom/left rule.
  std::pair<int, int> index_ij_unchecked(double val_out,
                                         double val_in) const noexcept {
    return {outter_.index(val_out), inner_.index(val_in)};
  }

  /// Unchecked row-major linear index (may be negative or out-of-range).
  index_type index(double val_out, double val_in) const noexcept {
    const auto ij = index_ij_unchecked(val_out, val_in);
    return static_cast<index_type>(ij.first) * static_cast<index_type>(cols()) +
           static_cast<index_type>(ij.second);
  }

  // ---- cells ---------------------------------------------------------------

  /**
   * @brief Return surrounding cell (validated). Throws if none exists.
   *
   * Uses bottom/left indices i=outer.index(val_out), j=inner.index(val_in),
   * then neighbors in each axis’ forward direction (sign-aware). For periodic
   * axes, neighbors wrap; for non-periodic, both the base index and its forward
   * neighbor must be within [0..N-1].
   *
   * Also verifies geometric containment along each axis using
   * TickAxis::within_segment_inclusive.
   */
  Cell cell(double val_out, double val_in) const {
    const int nrows = rows();
    const int ncols = cols();
    if (nrows < 2 || ncols < 2)
      throw std::runtime_error(
          "TickAxis2D::cell: grid too small to form cells");

    const int i = outter_.index(val_out);
    const int j = inner_.index(val_in);
    const int io = outter_.neighbor_forward(i);  // ±1, wrapped iff periodic
    const int ji = inner_.neighbor_forward(j);   // ±1, wrapped iff periodic

    // Range validation (compile-time: only for non-periodic axes)
    if constexpr (!AxisOut::periodic()) {
      if (i < 0 || i >= nrows || io < 0 || io >= nrows)
        throw std::out_of_range("TickAxis2D::cell: outer index out of range");
    }
    if constexpr (!AxisIn::periodic()) {
      if (j < 0 || j >= ncols || ji < 0 || ji >= ncols)
        throw std::out_of_range("TickAxis2D::cell: inner index out of range");
    }

    // Geometric containment along each axis (handles unwrap internally if
    // periodic)
    if (!outter_.within_segment_inclusive(val_out, i))
      throw std::out_of_range(
          "TickAxis2D::cell: outer value not enclosed by adjacent ticks");
    if (!inner_.within_segment_inclusive(val_in, j))
      throw std::out_of_range(
          "TickAxis2D::cell: inner value not enclosed by adjacent ticks");

    const index_type stride = static_cast<index_type>(ncols);
    const index_type bl =
        static_cast<index_type>(i) * stride + static_cast<index_type>(j);
    const index_type br =
        static_cast<index_type>(i) * stride + static_cast<index_type>(ji);
    const index_type tl =
        static_cast<index_type>(io) * stride + static_cast<index_type>(j);
    const index_type tr =
        static_cast<index_type>(io) * stride + static_cast<index_type>(ji);
    return Cell{bl, br, tl, tr};
  }

  /// Unchecked variant: may return out-of-range indices; no validation.
  Cell cell_nocheck(double val_out, double val_in) const noexcept {
    const int i = outter_.index(val_out);
    const int j = inner_.index(val_in);
    const int io = outter_.neighbor_forward(i);
    const int ji = inner_.neighbor_forward(j);

    const index_type stride = static_cast<index_type>(cols());
    const index_type bl =
        static_cast<index_type>(i) * stride + static_cast<index_type>(j);
    const index_type br =
        static_cast<index_type>(i) * stride + static_cast<index_type>(ji);
    const index_type tl =
        static_cast<index_type>(io) * stride + static_cast<index_type>(j);
    const index_type tr =
        static_cast<index_type>(io) * stride + static_cast<index_type>(ji);
    return Cell{bl, br, tl, tr};
  }
};

/// Handy aliases
using TickAxis2D_NN =
    TickAxis2D<IsPeriodic::NonPeriodic, IsPeriodic::NonPeriodic>;
using TickAxis2D_NP = TickAxis2D<IsPeriodic::NonPeriodic, IsPeriodic::Periodic>;
using TickAxis2D_PN = TickAxis2D<IsPeriodic::Periodic, IsPeriodic::NonPeriodic>;
using TickAxis2D_PP = TickAxis2D<IsPeriodic::Periodic, IsPeriodic::Periodic>;
} /* namespace dso */

#endif