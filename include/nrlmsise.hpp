#ifndef __DSO_NRLMSISE00_UPPERATMO_HPP__
#define __DSO_NRLMSISE00_UPPERATMO_HPP__

namespace dso {
class Nrlmsise00 {

  /* POWER7 */
  static constexpr double pt[150];
  static constexpr double pd[9][150];
  static constexpr double ps[150];
  static constexpr double pdl[2][25];
  static constexpr double ptl[4][100];
  static constexpr double pma[10][100];
  static constexpr double sam[100];
  /* LOWER7 */
  static constexpr double ptm[10];
  static constexpr double pdm[8][10];
  static constexpr double pavgm[10];


  /* TODO why are these here? */
  double gsurf;
  double re;
  /* GTS3C */
  double dd;
  /* DMIX */
  double dm04, dm16, dm28, dm32, dm40, dm01, dm14;
  /* MESO7 */
  double meso_tn1[5];
  double meso_tn2[4];
  double meso_tn3[5];
  double meso_tgn1[2];
  double meso_tgn2[2];
  double meso_tgn3[2];
  /* LPOLY */
  double dfa;
  double plg[4][9];
  double ctloc, stloc;
  double c2tloc, s2tloc;
  double s3tloc, c3tloc;
  double apdf, apt[4];
};/* Nrlmsise00 */
} /* namespace dso */

#endif
