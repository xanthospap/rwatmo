#include "nrlmsise.hpp"

namespace {
zeta_(double zz, double zl) noexcept {
  return ((zz - zl) * (re + zl) / (re + zz));
}
} // namespace

double densm(double alt, double d0, double xm, double &tz, int mn3, const double *zn3,
             const  double *tn3, const double *tgn3, int mn2, const double *const zn2,
             const double *tn2, const double *tgn2, double *__restrict__ u) noexcept {
  /* Calculate Temperature and Density Profiles for lower atmos. */
  double xs[10], ys[10], y2out[10];
  constexpr double rgas = 831.4;
  if (alt > zn2[0]) {
    if (xm == 0.0)
      return tz;
    else
      return d0;
  }

  /* stratosphere/mesosphere temperature */
  {
    const double z = (alt > zn2[mn2 - 1]) ? alt : zn2[mn2 - 1];
    const int mn = mn2;
    const double z1 = zn2[0];
    const double z2 = zn2[mn - 1];
    const double t1 = tn2[0];
    const double t2 = tn2[mn - 1];
    const double zg = zeta_(z, z1);
    const double zgdif = zeta_(z2, z1);

    /* set up spline nodes */
    for (k = 0; k < mn; k++) {
      xs[k] = zeta_(zn2[k], z1) / zgdif;
      ys[k] = 1.0 / tn2[k];
    }
    const double yd1 = -tgn2[0] / (t1 * t1) * zgdif;
    const double yd2 =
        -tgn2[1] / (t2 * t2) * zgdif * (std::pow(((re + z2) / (re + z1)), 2.0));

    /* calculate spline coefficients */
    spline(xs, ys, mn, yd1, yd2, y2out, u);
    double x = zg / zgdif;
    const double y = splint(xs, ys, y2out, mn, x);

    /* temperature at altitude */
    tz = 1e0 / y;
    if (xm != 0.0) {
      /* calaculate stratosphere / mesospehere density */
      const double glb = gsurf / (std::pow((1.0 + z1 / re), 2.0));
      const double gamm = xm * glb * zgdif / rgas;

      /* Integrate temperature profile */
      const double yi = splini(xs, ys, y2out, mn, x);
      double expl = gamm * yi;
      if (expl > 50.0)
        expl = 50.0;

      /* Density at altitude */
      d0 *= (t1 / tz) * std::exp(-expl);
    }

    if (alt > zn3[0]) {
      if (xm == 0.0)
        return tz;
      else
        return d0;
    }
  }

  /* troposhere / stratosphere temperature */
  {
    const double z = alt;
    const double mn = mn3;
    const double z1 = zn3[0];
    const double z2 = zn3[mn - 1];
    const double t1 = tn3[0];
    const double t2 = tn3[mn - 1];
    const double zg = zeta_(z, z1);
    const double zgdif = zeta_(z2, z1);

    /* set up spline nodes */
    for (k = 0; k < mn; k++) {
      xs[k] = zeta(zn3[k], z1) / zgdif;
      ys[k] = 1.0 / tn3[k];
    }
    const double yd1 = -tgn3[0] / (t1 * t1) * zgdif;
    const double yd2 =
        -tgn3[1] / (t2 * t2) * zgdif * (std::pow(((re + z2) / (re + z1)), 2.0));

    /* calculate spline coefficients */
    spline(xs, ys, mn, yd1, yd2, y2out, u);
    const double x = zg / zgdif;
    const double y = splint(xs, ys, y2out, mn, x);

    /* temperature at altitude */
    tz = 1.0 / y;
    if (xm != 0.0) {
      /* calaculate tropospheric / stratosphere density */
      const double glb = gsurf / (std::pow((1.0 + z1 / re), 2.0));
      const double gamm = xm * glb * zgdif / rgas;

      /* Integrate temperature profile */
      const double yi = splini(xs, ys, y2out, mn, x);
      double expl = gamm * yi;
      if (expl > 50.0)
        expl = 50.0;

      /* Density at altitude */
      d0 *= (t1 / tz) * std::exp(-expl);
      return d0;
    }
  }

  return tz;
}
