#ifndef __DSO_RADIO_WAVE_TROPOSPHERE_GEN_HPP__
#define __DSO_RADIO_WAVE_TROPOSPHERE_GEN_HPP__

#include <cmath>

namespace dso {
/** Zenith hydrostatic delay according to Saastamoinen and Davis.
 *
 * Saastamoinen, J., Atmospheric correction for the troposphere and 
 * stratosphere in radio ranging of satellites. The use of artificial 
 * satellites for geodesy, Geophys. Monogr. Ser. 15, Amer. Geophys. Union, 
 * pp. 274-251, 1972.
 *
 * Davis, J.L, T.A. Herring, I.I. Shapiro, A.E.E. Rogers, and G. Elgered, 
 * Geodesy by Radio Interferometry: Effects of Atmospheric Modeling Errors 
 * on Estimates of Baseline Length, Radio Science, Vol. 20, No. 6, 
 * pp. 1593-1607, 1985
 *
 * @param[in] p    Pressure at site [hPa]
 * @param[in] lat  Latitude of site [rad]
 * @param[in] hell Ellipsoidal height at site [m]
 * @return Zenith hydrostatic delay in [m]
 */
inline double zhd(double p, double lat, double hell) noexcept {
  return (0.0022768*p) / (1e0-0.00266 * std::cos(2e0*lat) - 0.28e-6*hell);
}

/** Zenith wet delay based on Aske and Nordius (1987)
 *
 * Askne and Nordius, Estimation of tropospheric delay for microwaves from 
 * surface weather data, Radio Science, Vol 22(3): 379-386, 1987.
 *
 * param[in] e   Water vapor pressure [hPa]
 * param[in] tm  Mean temperature [Kelvin]
 * param[in] lambda Water vapor lapse rate (see definition in Askne and  
 *                 Nordius 1987)
 * @return Zenith wet delay in [m]
 */
inline double zwd(double e, double tm, double lambda) noexcept  {
  constexpr const double k1 = 77.604e0;
  constexpr const double k2 = 64.790e0;
  constexpr const double k3 = 3.776e5;
  /* molar mass of water */
  constexpr const double Mw = 18.0152e0; 
  /* molar mass of dry air */
  constexpr const double Md = 28.9644e0; 
  /* k2' */
  constexpr const double k2d = k2 - k1 * Mw / Md;
  /* molar gas constant */
  constexpr const double R = 8.314e0; 
  /* specific gas constant for dry constituents */
  constexpr const double Rd = R / Md;
  /* mean gravity */
  constexpr const double gm = 9.80665e0;
  /* Eq. (18) from Askne and Nordius */
  return 1e-3 * (k2d+k3/tm) * (Rd*e) / ((lambda+1e0)*gm);
}
} /* namespace dso */

#endif
