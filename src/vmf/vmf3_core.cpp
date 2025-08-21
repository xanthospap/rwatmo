#include "iers/gravity.hpp"
#include "vmf3.hpp"

int dso::Vmf3::vmf3_spatial_coeffs_impl(
    const Eigen::Vector3d &rsta, dso::vmf3::Vmf3FullCoeffs &coeffs) noexcept {
  dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> V(
      MAX_DEGREE + 1, MAX_ORDER + 1);
  dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> W(
      MAX_DEGREE + 1, MAX_ORDER + 1);

  /* calculate Legendre polynomials */
  if (dso::gravity::sh_basis_cs_exterior(rsta.normalized(), MAX_DEGREE,
                                         MAX_ORDER, V, W)) {
    fprintf(stderr,
            "[ERROR] Failed computing spherical harmonics basis "
            "functions! (traceback: %s)\n",
            __func__);
    return 2;
  }

  /* reconstruct to nullify */
  coeffs = dso::vmf3::Vmf3FullCoeffs();

  /* determine the coefficients bh, bw, ch and cw */
  int i = 0;
  for (int n = 0; n <= MAX_DEGREE; n++) {
    for (int m = 0; m <= n; m++) {
      coeffs.bh_A0 += (anm_bh_A0()[i] * V(n, m) + bnm_bh_A0()[i] * W(n, m));
      coeffs.bh_A1 += (anm_bh_A1()[i] * V(n, m) + bnm_bh_A1()[i] * W(n, m));
      coeffs.bh_B1 += (anm_bh_B1()[i] * V(n, m) + bnm_bh_B1()[i] * W(n, m));
      coeffs.bh_A2 += (anm_bh_A2()[i] * V(n, m) + bnm_bh_A2()[i] * W(n, m));
      coeffs.bh_B2 += (anm_bh_B2()[i] * V(n, m) + bnm_bh_B2()[i] * W(n, m));

      coeffs.bw_A0 += (anm_bw_A0()[i] * V(n, m) + bnm_bw_A0()[i] * W(n, m));
      coeffs.bw_A1 += (anm_bw_A1()[i] * V(n, m) + bnm_bw_A1()[i] * W(n, m));
      coeffs.bw_B1 += (anm_bw_B1()[i] * V(n, m) + bnm_bw_B1()[i] * W(n, m));
      coeffs.bw_A2 += (anm_bw_A2()[i] * V(n, m) + bnm_bw_A2()[i] * W(n, m));
      coeffs.bw_B2 += (anm_bw_B2()[i] * V(n, m) + bnm_bw_B2()[i] * W(n, m));

      coeffs.ch_A0 += (anm_ch_A0()[i] * V(n, m) + bnm_ch_A0()[i] * W(n, m));
      coeffs.ch_A1 += (anm_ch_A1()[i] * V(n, m) + bnm_ch_A1()[i] * W(n, m));
      coeffs.ch_B1 += (anm_ch_B1()[i] * V(n, m) + bnm_ch_B1()[i] * W(n, m));
      coeffs.ch_A2 += (anm_ch_A2()[i] * V(n, m) + bnm_ch_A2()[i] * W(n, m));
      coeffs.ch_B2 += (anm_ch_B2()[i] * V(n, m) + bnm_ch_B2()[i] * W(n, m));

      coeffs.cw_A0 += (anm_cw_A0()[i] * V(n, m) + bnm_cw_A0()[i] * W(n, m));
      coeffs.cw_A1 += (anm_cw_A1()[i] * V(n, m) + bnm_cw_A1()[i] * W(n, m));
      coeffs.cw_B1 += (anm_cw_B1()[i] * V(n, m) + bnm_cw_B1()[i] * W(n, m));
      coeffs.cw_A2 += (anm_cw_A2()[i] * V(n, m) + bnm_cw_A2()[i] * W(n, m));
      coeffs.cw_B2 += (anm_cw_B2()[i] * V(n, m) + bnm_cw_B2()[i] * W(n, m));

      ++i;
    }
  }

  return 0;
}

void dso::Vmf3::mf(const dso::MjdEpoch &t, double el, double ah, double aw,
                   const dso::vmf3::Vmf3FullCoeffs &bc, double &mfh,
                   double &mfw) noexcept {
  /* recompute temporal dependence if needed (for b's and c's) */
  if (std::abs(t.diff<DateTimeDifferenceType::FractionalSeconds>(last_epoch_)
                   .seconds()) > MAX_SEC_DIFF) {
    this->date_trigs(t);
  }

  /* add the seasonal amplitudes for the specified doy to the mean values */
  const double bh = bc.bh_A0 + bc.bh_A1 * cd2pi_ + bc.bh_B1 * sd2pi_ +
                    bc.bh_A2 * cd4pi_ + bc.bh_B2 * sd4pi_;
  const double bw = bc.bw_A0 + bc.bw_A1 * cd2pi_ + bc.bw_B1 * sd2pi_ +
                    bc.bw_A2 * cd4pi_ + bc.bw_B2 * sd4pi_;
  const double ch = bc.ch_A0 + bc.ch_A1 * cd2pi_ + bc.ch_B1 * sd2pi_ +
                    bc.ch_A2 * cd4pi_ + bc.ch_B2 * sd4pi_;
  const double cw = bc.cw_A0 + bc.cw_A1 * cd2pi_ + bc.cw_B1 * sd2pi_ +
                    bc.cw_A2 * cd4pi_ + bc.cw_B2 * sd4pi_;

  /* calculating the hydrostatic and wet mapping factors */
  const double sel = std::sin(el);
  mfh =
      (1 + (ah / (1 + bh / (1 + ch)))) / (sel + (ah / (sel + bh / (sel + ch))));
  mfw =
      (1 + (aw / (1 + bw / (1 + cw)))) / (sel + (aw / (sel + bw / (sel + cw))));
}