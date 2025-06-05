#include <cassert>
#include "nrlmsise.hpp"

double dso::Nrlmsise00::splini(const double *__restrict__ xa, const double *__restrict__ ya,
              const double *__restrict__ y2a, int n, double x) noexcept {
  /*      INTEGRATE CUBIC SPLINE FUNCTION FROM XA(1) TO X
   *       XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
   *       Y2A: ARRAY OF SECOND DERIVATIVES
   *       N: SIZE OF ARRAYS XA,YA,Y2A
   *       X: ABSCISSA ENDPOINT FOR INTEGRATION
   *       Y: OUTPUT VALUE
   */
  double yi = 0;
  int klo = 0;
  int khi = 1;
  while ((x > xa[klo]) && (khi < n)) {
    double xx = x;
    if (khi < (n - 1)) {
      if (x < xa[khi])
        xx = x;
      else
        xx = xa[khi];
    }
    const double h = xa[khi] - xa[klo];
    const double a = (xa[khi] - xx) / h;
    const double b = (xx - xa[klo]) / h;
    const double a2 = a * a;
    const double b2 = b * b;
    yi += ((1.0 - a2) * ya[klo] / 2.0 + b2 * ya[khi] / 2.0 +
           ((-(1.0 + a2 * a2) / 4.0 + a2 / 2.0) * y2a[klo] +
            (b2 * b2 / 4.0 - b2 / 2.0) * y2a[khi]) *
               h * h / 6.0) *
          h;
    klo++;
    khi++;
  }

  return yi;
}

double dso::Nrlmsise00::splint(const double *__restrict__ xa, const double *__restrict__ ya,
              const double *__restrict__ y2a, int n, double x) noexcept {
  /*      CALCULATE CUBIC SPLINE INTERP VALUE
   *       ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL.
   *       XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
   *       Y2A: ARRAY OF SECOND DERIVATIVES
   *       N: SIZE OF ARRAYS XA,YA,Y2A
   *       X: ABSCISSA FOR INTERPOLATION
   *       Y: OUTPUT VALUE
   */
  int klo = 0;
  int khi = n - 1;
  int k;
  double h;
  double a, b, yi;

  while ((khi - klo) > 1) {
    const int k = (khi + klo) / 2;
    if (xa[k] > x)
      khi = k;
    else
      klo = k;
  }

  const double h = xa[khi] - xa[klo];
  assert(h != 0e0);
  const double a = (xa[khi] - x) / h;
  const double b = (x - xa[klo]) / h;
  return a * ya[klo] + b * ya[khi] +
         ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * h * h /
             6.0;
}

/* size of u is n*sizeof(double) */
void dso::Nrlmsise00::spline(const double *__restrict__ x, const double *__restrict__ y, int n,
            double yp1, double ypn, const double *__restrict__ y2,
            double *__restrict__ u) noexcept {
  /* CALCULATE 2ND DERIVATIVES OF CUBIC SPLINE INTERP FUNCTION
   * ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL
   * X,Y: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
   * N: SIZE OF ARRAYS X,Y
   * YP1,YPN: SPECIFIED DERIVATIVES AT X[0] AND X[N-1]; VALUES
   *          >= 1E30 SIGNAL SIGNAL SECOND DERIVATIVE ZERO
   * Y2: OUTPUT ARRAY OF SECOND DERIVATIVES
   */
  double sig, p, qn, un;
  int i, k;

  if (yp1 > 0.99e30) {
    y2[0] = 0;
    u[0] = 0;
  } else {
    y2[0] = -0.5;
    u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
  }

  double qn, un;
  if (ypn > 0.99e30) {
    qn = 0;
    un = 0;
  } else {
    qn = 0.5;
    un = (3.0 / (x[n - 1] - x[n - 2])) *
         (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
  }

  for (i = 1; i < (n - 1); i++) {
    const double sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
    const double p = sig * y2[i - 1] + 2.0;
    y2[i] = (sig - 1.0) / p;
    u[i] = (6.0 *
                ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) -
                 (y[i] - y[i - 1]) / (x[i] - x[i - 1])) /
                (x[i + 1] - x[i - 1]) -
            sig * u[i - 1]) /
           p;
  }

  y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);
  for (k = n - 2; k >= 0; k--)
    y2[k] = y2[k] * y2[k + 1] + u[k];

  return;
}
