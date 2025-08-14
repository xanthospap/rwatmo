#ifndef __DSO_TUWIEN_VMF_TROPO_DELAY_HPP__
#define __DSO_TUWIEN_VMF_TROPO_DELAY_HPP__

#include "datetime/calendar.hpp"
#include "geodesy/transformations.hpp"
#include "iers/gravity.hpp"
#include <fstream>
#include <string>

namespace dso {

namespace vmf3 {

/** @brief VMF data file hold site names with 4 characters */
struct SiteName {
  char data[5] = {'\0'};
  bool operator==(const SiteName &other) const noexcept {
    return (!std::strcmp(data, other.data));
  }
  bool operator!=(const SiteName &other) const noexcept {
    return !(*this == other);
  }
};

struct SiteBlock {
  static constexpr const int NUM_DATA_FIELDS = 11;
  SiteName msite;
  double mdata[NUM_DATA_FIELDS];
  // hydrostatic "a" coefficient
  // wet "a" coefficient
  // zenith hydrostatic delay [m]
  // zenith wet delay [m]
  // pressure at the site [hPa]
  // temperature at the site [°C]
  // water vapor pressure at the site [hPa]
  // hydrostatic north gradient Gn_h [mm]
  // hydrostatic east gradient Ge_h [mm]
  // wet north gradient Gn_w [mm]
  // wet east gradient Ge_w [mm]

  const double *data() const noexcept { return mdata; }
  double *data() noexcept { return mdata; }
  const char *const site() const noexcept { return msite.data; }
  char *site() noexcept { return msite.data; }
}; /* struct SiteBlock */

/** @brief A whole block of site-specific VMF? parameters for a given epoch.
 */
struct Block {
  /* reference epoch */
  MjdEpoch mt;
  /* site-specific records */
  std::vector<SiteBlock> msite_data;

  void append(const SiteBlock &b) noexcept { msite_data.emplace_back(b); }

  void reset() {
    msite_data.clear();
    mt = dso::MjdEpoch::min();
  }
}; /* struct Block */

/** Full set of coefficients for computing b and c 'empirical' VMF3
 * coefficients. These coefficients have only a spatial dependence, hence can
 * be computed once for every site of interest.
 *
 * To fully compute the coefficients, temporal dependence must be considered,
 * which is done in the function Vmf3::compute().
 */
struct Vmf3FullCoeffs {
  double bh_A0{0e0};
  double bh_A1{0e0};
  double bh_B1{0e0};
  double bh_A2{0e0};
  double bh_B2{0e0};
  double bw_A0{0e0};
  double bw_A1{0e0};
  double bw_B1{0e0};
  double bw_A2{0e0};
  double bw_B2{0e0};
  double ch_A0{0e0};
  double ch_A1{0e0};
  double ch_B1{0e0};
  double ch_A2{0e0};
  double ch_B2{0e0};
  double cw_A0{0e0};
  double cw_A1{0e0};
  double cw_B1{0e0};
  double cw_A2{0e0};
  double cw_B2{0e0};
}; /* Vmf3FullCoeffs*/

} /* namespace vmf3 */

class Vmf3SiteStream {
  /* max characters in a VMF3 data file line */
  static constexpr int MAX_RECORD_CHARS = 256;
  /* the input stream */
  std::ifstream mstream;
  /* last line read, first to be resolved for new block */
  char mlast_line[MAX_RECORD_CHARS] = {'\0'};
  /* last block read, NOT current */
  vmf3::Block mlast_block;
  /* iterator, holds the current block */
  iterator mit;
  /* filename */
  std::string mfilename;

  int get_next_block(vmf3::Block &b) noexcept;

  Vmf3SiteStream(const char *fn)
      : mstream(fn), mlast_block(), mit(), mfilename(fn) {
    if (!mstream.is_open()) {
      throw std::runtime_error(
          "[ERROR] Failed opening VMF3 site-specific file " + mfilename + "\n");
    }
  }

  class iterator {
  public:
    using value_type = vmf3::Block;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::input_iterator_tag;
    using pointer = const vmf3::Block *;
    using reference = const vmf3::Block &;

    iterator() noexcept = default; /* end() */
    reference operator*() const noexcept { return current_; }
    pointer operator->() const noexcept { return &current_; }

    /* pre-increment: advance to next block */
    iterator &operator++() {
      /* already end() */
      if (!owner_)
        return *this;

      int rc = owner_->get_next_block(current_);
      if (rc == 0) {
        return *this; /* ok, now positioned at next block */
      }
      /* reached EOF, set iterator to end() */
      if (rc < 0) {
        owner_ = nullptr;
        return *this;
      }
      /* rc > 0 → error */
      throw std::runtime_error(
          "[ERROR] Failed getting next data block from MyDataFile");
    }

    /* post-increment: returns the old value (input iterator semantics) */
    iterator operator++(int) {
      iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    friend bool operator==(const iterator &a, const iterator &b) noexcept {
      return a.owner_ == b.owner_;
    }
    friend bool operator!=(const iterator &a, const iterator &b) noexcept {
      return !(a == b);
    }

  private:
    explicit iterator(MyDataFile *owner) : owner_(owner) {
      // prime first block
      int rc = owner_->get_next_block(current_);
      if (rc == 0)
        return;
      if (rc < 0) {
        owner_ = nullptr;
        return;
      } // empty file → end()
      throw std::runtime_error(
          "[ERROR] Failed getting first data block from MyDataFile");
    }

    /* non-owning; end() when nullptr */
    Vmf3SiteStream *owner_ = nullptr;
    /* current block (all sites) */
    vmf3::Block current_;
  }; /* class Vmf3SiteStream::iterator */

}; /* class Vmf3SiteStream */

class Vmf3 {
  /* coefficients for SH expansion (spatial dependence) */
  static constexpr const double anm_bh[5][91];
  static constexpr const double anm_bw[5][91];
  static constexpr const double anm_ch[5][91];
  static constexpr const double anm_cw[5][91];
  static constexpr const double bnm_bh[5][91];
  static constexpr const double bnm_bw[5][91];
  static constexpr const double bnm_ch[5][91];
  static constexpr const double bnm_cw[5][91];

  /* degree and order of SH expansion (spatial dependence) */
  static constexpr const int MAX_DEGREE = 12;
  static constexpr const int MAX_ORDER = 12;

  /* last epoch used for computing temporal dependence, i.e. trigs based on
   * date of interest
   */
  static MjdEpoch last_epoch = MjdEpoch::min();
  /* if two epochs are less than MAX_SEC_DIFF away, do not recompute temporal
   * dependence. use the trigs as stored in [cs]d[24]pi
   */
  static constexpr const double MAX_SEC_DIFF = 60;
  /* coefficients (i.e. trigs) for temporal dependence. computed at last_epoch
   */
  static double sd2pi;
  static double cd2pi;
  static double sd4pi;
  static double cd4pi;

  /* recompute temporal dependence, i.e. trigs [cs]d[24]pi */
  static void date_trigs(const MjdEpoch &t) noexcept {
    const dso::ydoy_date yd(t.to_ydoy());
    const double doy =
        yd.dy().as_underlying_type() + t.fractional_days().days();
    sd2pi = std::sin(doy / 365.25 * 2 * dso::DPI);
    cd2pi = std::cos(doy / 365.25 * 2 * dso::DPI);
    sd4pi = std::sin(doy / 365.25 * 4 * dso::DPI);
    cd4pi = std::cos(doy / 365.25 * 4 * dso::DPI);
    last_epoch = t;
    return;
  }

  /* easy access for coefficient arrays */
  static const double *const anm_bh_A0() noexcept { return anm_bh[1 - 1]; }
  static const double *const anm_bh_A1() noexcept { return anm_bh[2 - 1]; }
  static const double *const anm_bh_B1() noexcept { return anm_bh[3 - 1]; }
  static const double *const anm_bh_A2() noexcept { return anm_bh[4 - 1]; }
  static const double *const anm_bh_B2() noexcept { return anm_bh[5 - 1]; }
  static const double *const anm_bw_A0() noexcept { return anm_bw[1 - 1]; }
  static const double *const anm_bw_A1() noexcept { return anm_bw[2 - 1]; }
  static const double *const anm_bw_B1() noexcept { return anm_bw[3 - 1]; }
  static const double *const anm_bw_A2() noexcept { return anm_bw[4 - 1]; }
  static const double *const anm_bw_B2() noexcept { return anm_bw[5 - 1]; }
  static const double *const anm_ch_A0() noexcept { return anm_ch[1 - 1]; }
  static const double *const anm_ch_A1() noexcept { return anm_ch[2 - 1]; }
  static const double *const anm_ch_B1() noexcept { return anm_ch[3 - 1]; }
  static const double *const anm_ch_A2() noexcept { return anm_ch[4 - 1]; }
  static const double *const anm_ch_B2() noexcept { return anm_ch[5 - 1]; }
  static const double *const anm_cw_A0() noexcept { return anm_cw[1 - 1]; }
  static const double *const anm_cw_A1() noexcept { return anm_cw[2 - 1]; }
  static const double *const anm_cw_B1() noexcept { return anm_cw[3 - 1]; }
  static const double *const anm_cw_A2() noexcept { return anm_cw[4 - 1]; }
  static const double *const anm_cw_B2() noexcept { return anm_cw[5 - 1]; }
  static const double *const bnm_bh_A0() noexcept { return bnm_bh[1 - 1]; }
  static const double *const bnm_bh_A1() noexcept { return bnm_bh[2 - 1]; }
  static const double *const bnm_bh_B1() noexcept { return bnm_bh[3 - 1]; }
  static const double *const bnm_bh_A2() noexcept { return bnm_bh[4 - 1]; }
  static const double *const bnm_bh_B2() noexcept { return bnm_bh[5 - 1]; }
  static const double *const bnm_bw_A0() noexcept { return bnm_bw[1 - 1]; }
  static const double *const bnm_bw_A1() noexcept { return bnm_bw[2 - 1]; }
  static const double *const bnm_bw_B1() noexcept { return bnm_bw[3 - 1]; }
  static const double *const bnm_bw_A2() noexcept { return bnm_bw[4 - 1]; }
  static const double *const bnm_bw_B2() noexcept { return bnm_bw[5 - 1]; }
  static const double *const bnm_ch_A0() noexcept { return bnm_ch[1 - 1]; }
  static const double *const bnm_ch_A1() noexcept { return bnm_ch[2 - 1]; }
  static const double *const bnm_ch_B1() noexcept { return bnm_ch[3 - 1]; }
  static const double *const bnm_ch_A2() noexcept { return bnm_ch[4 - 1]; }
  static const double *const bnm_ch_B2() noexcept { return bnm_ch[5 - 1]; }
  static const double *const bnm_cw_A0() noexcept { return bnm_cw[1 - 1]; }
  static const double *const bnm_cw_A1() noexcept { return bnm_cw[2 - 1]; }
  static const double *const bnm_cw_B1() noexcept { return bnm_cw[3 - 1]; }
  static const double *const bnm_cw_A2() noexcept { return bnm_cw[4 - 1]; }
  static const double *const bnm_cw_B2() noexcept { return bnm_cw[5 - 1]; }

public:
  /** @brief This function will compute the spatial part of the  b and c
   * 'empirical' VMF3 coefficients
   *
   * These coefficients depend only on the coordinates of the site, and do
   * not involve any temporal changes. Hence, they can be computed once for
   * every site of interest.
   *
   * @param[in] rsta Cartesian ECEF coordinates of site [m]
   * @param[out] coeffs An instance of vmf3::Vmf3FullCoeffs which can be later
   * used to compute the VMF3 mapping functions via mf()
   */
  static void vmf3_spatial_coeffs(dso::CartesianCrdConstView &rsta,
                                  vmf3::Vmf3FullCoeffs &coeffs) noexcept {
    dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> V(
        MAX_DEGREE, MAX_ORDER);
    dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> W(
        MAX_DEGREE, MAX_ORDER);

    const Eigen::Vector3d r = rsta.mv;

    /* calculate Legendre polynomials */
    if (dso::gravity::sh_basis_cs_exterior(r / r.normalized(), MAX_DEGREE,
                                           MAX_ORDER, V, W)) {
      fprintf(stderr,
              "[ERROR] Failed computing spherical harmonics basis "
              "functions! (traceback: %s)\n",
              __func__);
      return 2;
    }

    /* reconstruct to nullify */
    coeffs = vmf3::Vmf3FullCoeffs();

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

    return;
  }

  /** @brief Compute the VMF3 mapping functions, mfh and mfw.
   *
   * @param[in] t The epoch of interest
   * @param[in] el Elevation angle in [rad]
   * @param[in] bc Spatial b and c empirical coefficients as
   * vmf3::Vmf3FullCoeffs instance for the site of interest. These can be
   * computed via a call to vmf3_spatial_coeffs
   * @param[out] mfh Value of the VMF3 hydrostatic mapping function [m]
   * @param[out] mfw Value of the VMF3 wet mapping function [m]
   */
  void mf(MjdEpoch &t, double el, const vmf3::Vmf3FullCoeffs &bc, double &mfh,
          double &mfw) noexcept {

    if (std::abs(t.diff<DateTimeDifferenceType::FractionalSeconds>(last_epoch)
                     .seconds()) > MAX_SEC_DIFF) {
      date_trigs(t);
    }

    /* adding the seasonal amplitudes for the specified doy to the mean values
     */
    const double bh = bc.bh_A0 + bc.bh_A1 * cd2pi + bc.bh_B1 * sd2pi +
                      bc.bh_A2 * cd4pi + bc.bh_B2 * sd4pi;
    const double bw = bc.bw_A0 + bc.bw_A1 * cd2pi + bc.bw_B1 * sd2pi +
                      bc.bw_A2 * cd4pi + bc.bw_B2 * sd4pi;
    const double ch = bc.ch_A0 + bc.ch_A1 * cd2pi + bc.ch_B1 * sd2pi +
                      bc.ch_A2 * cd4pi + bc.ch_B2 * sd4pi;
    const double cw = bc.cw_A0 + bc.cw_A1 * cd2pi + bc.cw_B1 * sd2pi +
                      bc.cw_A2 * cd4pi + bc.cw_B2 * sd4pi;

    /* calculating the hydrostatic and wet mapping factors */
    const double sel = std::sin(el);
    mfh = (1 + (ah / (1 + bh / (1 + ch)))) /
          (sel + (ah / (sel + bh / (sel + ch))));
    mfw = (1 + (aw / (1 + bw / (1 + cw)))) /
          (sel + (aw / (sel + bw / (sel + cw))));
  }

}; /* struct Vmf3 */

} /* namespace dso */

#endif
