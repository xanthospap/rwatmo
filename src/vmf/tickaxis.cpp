#include <cmath>
#include <stdexcept>

#include "vmf3_grid_stream.hpp"

dso::vmf3::TickAxis::TickAxis(double start, double stop, double step,
                              double eps_steps) {
  if (!std::isfinite(start) || !std::isfinite(stop) || !std::isfinite(step))
    throw std::runtime_error(
        "[ERROR] non-finite start/stop/step (traceback: "
        "vmf3::TickAxis::TickAxis)\n");

  if (step == 0.0)
    throw std::runtime_error(
        "[ERROR] step cannot be zero (traceback: vmf3::TickAxis::TickAxis)\n");

  if (eps_steps < 0.0)
    throw std::runtime_error(
        "[ERROR] eps_steps must be non-negative (traceback: "
        "vmf3::TickAxis::TickAxis)\n");

  start_ = start;

  /* Ensure step direction matches the axis direction. */
  const double dir = (stop >= start) ? +1.0 : -1.0;
  step_ = dir * std::abs(step);

  /* Theoretical number of steps */
  const double delta = stop - start_;
  const double raw = delta / step_;

  /* Nearest integer number of steps */
  const long long n = static_cast<long long>(std::llround(raw));

  /* Reconstruct endpoint and check divisibility within tolerance */
  const double snapped_stop = start_ + static_cast<double>(n) * step_;
  const double tol =
      std::abs(step_) * eps_steps;  // absolute tolerance in axis units
  if (std::abs(snapped_stop - stop) > tol) {
    throw std::runtime_error(
        "[ERROR] (stop - start) is not an integer multiple of step within "
        "tolerance (traceback: vmf3::TickAxis::TickAxis)\n");
  }

  /* Inclusive count */
  const long long pts_ll = n + 1;
  if (pts_ll <= 0 ||
      pts_ll > static_cast<long long>(std::numeric_limits<int>::max()))
    throw std::runtime_error(
        "[ERROR] invalid number of points after rounding (traceback: "
        "vmf3::TickAxis::TickAxis)\n");

  num_pts_ = static_cast<int>(pts_ll);
}

int dso::vmf3::TickAxis::index(double val) const noexcept {
  const double u = (val - start_) / step_;
  /* Nudge exact-grid hits to the lower cell to avoid ambiguity */
  constexpr double eps = 1e-12;
  if (step_ > 0) {
    /* coordinate increases with index; "bottom/left" => floor */
    return static_cast<int>(std::floor(u - eps));
  } else {
    /* coordinate decreases with index; "bottom/left" => ceil */
    return static_cast<int>(std::ceil(u + eps));
  }
}