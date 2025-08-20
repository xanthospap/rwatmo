#ifndef __DSO__GENERIC_GRID_1D_CLASS_HPP__
#define __DSO__GENERIC_GRID_1D_CLASS_HPP__

#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace dso {
enum class IsPeriodic { NonPeriodic, Periodic };

/* Storage present only for the periodic specialization */
template <IsPeriodic P>
struct PeriodicStorage {};
template <>
struct PeriodicStorage<IsPeriodic::Periodic> {
  double period_{0.0};
};

/** @brief Discrete, inclusive 1-D axis with robust numerics.
 * @tparam P  IsPeriodic::Periodic for seam-aware axes (e.g., longitude),
 *            IsPeriodic::NonPeriodic otherwise.
 *
 * Inclusive means both start and stop are ticks. The sign of @p step
 * is normalized to match (stop - start). Divisibility is validated
 * within a tolerance.
 */
template <IsPeriodic P>
class TickAxis : private PeriodicStorage<P> {
  static constexpr bool kPeriodic = (P == IsPeriodic::Periodic);
  using PS = PeriodicStorage<P>;

  double start_{0.0};  ///< first tick (inclusive)
  double step_{1.0};   ///< signed spacing
  int num_pts_{0};     ///< inclusive number of ticks

 public:
  // accessor/setter that compile only the relevant branch
  constexpr double period() const noexcept {
    if constexpr (kPeriodic)
      return PS::period_;
    else
      return 0.0;  // or throw/unused
  }
  constexpr void set_period(double T) noexcept {
    if constexpr (kPeriodic)
      PS::period_ = T;
    else
      (void)T;  // no-op
  }

  /// Construct [start..stop] inclusive with spacing @p step.
  /// @param eps_steps tolerance as a fraction of |step| (default 1e-9)
  /// @param period_hint required period if periodic (e.g., 360.0); if 0, infer
  TickAxis(double start, double stop, double step, double eps_steps = 1e-12,
           double period_hint = 0.0) {
    if (step == 0.0) throw std::runtime_error("TickAxis: step cannot be zero");
    if (eps_steps < 0.0) throw std::runtime_error("TickAxis: eps_steps < 0");

    start_ = start;
    const double dir = (stop >= start) ? +1.0 : -1.0;
    step_ = dir * std::abs(step);

    // robust count (inclusive)
    const double raw = (stop - start_) / step_;
    const long n = static_cast<long>(std::llround(raw));
    const long pts_ll = n + 1;
    if (pts_ll <= 0 ||
        pts_ll > static_cast<long>(std::numeric_limits<int>::max()))
      throw std::runtime_error(
          "TickAxis: invalid number of points after rounding");
    num_pts_ = static_cast<int>(pts_ll);

    // divisibility check
    const double snapped_stop = start_ + static_cast<double>(n) * step_;
    const double tol = std::abs(step_) * eps_steps;
    if (std::abs(snapped_stop - stop) > tol)
      throw std::runtime_error(
          "TickAxis: (stop-start) not a multiple of step within tolerance");

    if constexpr (kPeriodic) {
      if (period_hint > 0.0) {
        set_period(period_hint);
        const double inferred = std::abs(step_) * static_cast<double>(num_pts_);
        if (std::abs(inferred - period()) >
            std::abs(step_) * eps_steps * num_pts_)
          throw std::runtime_error(
              "TickAxis: period_hint inconsistent with step*num_pts");
      } else {
        set_period(std::abs(step_) *
                   static_cast<double>(num_pts_));  // e.g., 1°*360 = 360°
      }
    }
  }

  /// Default: empty axis.
  TickAxis() = default;

  // ---------- trivial getters ----------
  double tick_start() const noexcept { return start_; }
  double tick_step() const noexcept { return step_; }
  double tick_stop() const noexcept { return start_ + step_ * (num_pts_ - 1); }
  int num_pts() const noexcept { return num_pts_; }

  // ---------- periodic helpers ----------
  static constexpr bool periodic() noexcept { return kPeriodic; }

  /// Wrap index into [0, num_pts()-1] when periodic; otherwise return as-is.
  int wrap_index(int idx) const noexcept {
    if constexpr (!kPeriodic) return idx;
    const int n = num_pts_;
    const int m = idx % n;
    return (m < 0) ? (m + n) : m;
  }

  /// Forward neighbor (±1 in axis direction). Wrapped if periodic.
  int neighbor_forward(int idx) const noexcept {
    const int s = (step_ > 0) ? +1 : -1;
    if constexpr (kPeriodic)
      return wrap_index(idx + s);
    else
      return idx + s;
  }

  /// Wrap a value into the canonical interval when periodic; otherwise return
  /// it.
  double wrap_value(double v) const noexcept {
    if constexpr (!kPeriodic) return v;
    const double N = static_cast<double>(num_pts_);
    double u = (v - start_) / step_;  // normalized in ticks
    u -= N * std::floor(u / N);       // reduce to [0,N)
    return start_ + u * step_;        // back to axis units
  }

  // ---------- main mapping ----------
  /// Map a coordinate to the bottom/left tick index (no bounds checks).
  int index(double v) const noexcept {
    if constexpr (kPeriodic) v = wrap_value(v);
    const double u = (v - start_) / step_;
    constexpr double eps = 1e-12;
    int k = (step_ > 0) ? static_cast<int>(std::floor(u - eps))
                        : static_cast<int>(std::ceil(u + eps));
    if constexpr (kPeriodic)
      return wrap_index(k);
    else
      return k;
  }

  /// Map integer index to tick value (no bounds checks).
  double val(int idx) const noexcept { return start_ + step_ * idx; }

  /// Non-periodic: check against span; Periodic: always true.
  bool in_range(double v) const noexcept {
    if constexpr (kPeriodic) return true;
    const double lo = std::min(start_, tick_stop());
    const double hi = std::max(start_, tick_stop());
    return (v >= lo) && (v <= hi);
  }

  /// True if v lies between idx0 and its forward neighbor. Periodic: unwrapped
  /// around v.
  bool within_segment_inclusive(double v, int idx0) const noexcept {
    const int idx1 = neighbor_forward(idx0);
    double a0 = val(idx0), a1 = val(idx1);
    if constexpr (kPeriodic) unwrap_pair_around(v, a0, a1);
    const double lo = std::min(a0, a1);
    const double hi = std::max(a0, a1);
    const double tol = std::abs(step_) * 1e-12;
    return (v >= lo - tol) && (v <= hi + tol);
  }

  // Does the same as unwrap_pair_around but also modifies the value if needed.
  // A no-op if the Axis is not periodic.
  // usefull e.g. when interpolating near bounds (for the periodic case)
  // Example:
  //  Periodic longitude axis: period = 360°, ticks: 0, 2.5, 5, …, 357.5 (step
  //  +2.5). Site longitude given as v = -1° (same as 359° but in a different
  //  wrap). The two neighbor grid longitudes for the cell you’re interpolating
  //  across are a0 = 357.5° and a1 = 0.0° (this segment crosses the seam).
  //  double v  = -1.0;     // site lon
  //  double a0 = 357.5;  // left grid lon
  //  double a1 = 0.0;    // right grid lon
  //  axis.unwrap_segment_around(v, a0, a1);
  //  Result : (v, a0, a1) = (-1.0, -2.5, 0.0) — a clean, monotonic segment with
  //  v inside it.
  void unwrap_segment_around(double& v, double& a0, double& a1) const noexcept {
    if constexpr (kPeriodic) {
      // reuse the private logic:
      // shift (a0,a1) near v and ensure segment direction matches step sign,
      // also nudge v into that unwrapped interval.
      const double per = period();
      const double k = std::floor((v - a0) / per + 0.5);
      a0 += k * per;
      a1 += k * per;
      if (step_ > 0) {
        if (a1 < a0) a1 += per;
        if (v < a0) v += per;
      } else {
        if (a1 > a0) a1 -= per;
        if (v > a0) v -= per;
      }
    } else {
      // non-periodic: nothing to do
      (void)a0;
      (void)a1;
      (void)v;
    }
  }

 private:
  // For periodic axes: shift (a0,a1) by integer multiples of period so that v
  // lies between them.
  void unwrap_pair_around(double v, double& a0, double& a1) const noexcept {
    if constexpr (!kPeriodic) return;
    const double per = period();
    const double k = std::floor((v - a0) / per + 0.5);
    a0 += k * per;
    a1 += k * per;
    if (step_ > 0) {
      if (a1 < a0) a1 += per;
      if (v < a0) v += per;  // local
    } else {
      if (a1 > a0) a1 -= per;
      if (v > a0) v -= per;  // local
    }
  }
};

/// Convenient aliases
using TickAxisNonPeriodic = TickAxis<IsPeriodic::NonPeriodic>;
using TickAxisPeriodic = TickAxis<IsPeriodic::Periodic>;

} /* namespace dso */

#endif
