#ifndef __DSO_RWATMO_CORE_MATH_UTILS_HPP__
#define __DSO_RWATMO_CORE_MATH_UTILS_HPP__

namespace dso {
  double splini(const double *__restrict__ xa, const double *__restrict__ ya,
                const double *__restrict__ y2a, int n, double x) noexcept;
  double splint(const double *__restrict__ xa, const double *__restrict__ ya,
                const double *__restrict__ y2a, int n, double x) noexcept;
  void spline(const double *__restrict__ x, const double *__restrict__ y, int n,
              double yp1, double ypn, double *__restrict__ y2,
              double *__restrict__ u) noexcept;
} /* namespace dso */

#endif
