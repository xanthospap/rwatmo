#include <datetime/core/fundamental_types_generic_utilities.hpp>
#include <geodesy/core/crdtype_warppers.hpp>
#include <limits>

#include "geodesy/transformations.hpp"
#include "vmf3_grid_data.hpp"

namespace dso {

namespace vmf3 {

/** Full set of coefficients for computing b and c 'empirical' VMF3
 * coefficients. These coefficients have only a spatial dependence, hence can
 * be computed once for every site of interest.
 *
 * To fully compute the coefficients, temporal dependence must be considered.
 * This can be done via the Vmf3::vmf3_spatial_coeffs function.
 */
struct Vmf3FullCoeffs {
  static constexpr const double MISSING = std::numeric_limits<double>::min();
  double bh_A0{MISSING};
  double bh_A1{MISSING};
  double bh_B1{MISSING};
  double bh_A2{MISSING};
  double bh_B2{MISSING};
  double bw_A0{MISSING};
  double bw_A1{MISSING};
  double bw_B1{MISSING};
  double bw_A2{MISSING};
  double bw_B2{MISSING};
  double ch_A0{MISSING};
  double ch_A1{MISSING};
  double ch_B1{MISSING};
  double ch_A2{MISSING};
  double ch_B2{MISSING};
  double cw_A0{MISSING};
  double cw_A1{MISSING};
  double cw_B1{MISSING};
  double cw_A2{MISSING};
  double cw_B2{MISSING};
  bool is_missing() const noexcept { return bh_A0 == MISSING; }
  void set_missing() noexcept { bh_A0 = MISSING; }
}; /* Vmf3FullCoeffs*/
} /* namespace vmf3 */

class Vmf3 {
  /* coefficients for SH expansion (spatial dependence) */
  static const double anm_bh[5][91];
  static const double anm_bw[5][91];
  static const double anm_ch[5][91];
  static const double anm_cw[5][91];
  static const double bnm_bh[5][91];
  static const double bnm_bw[5][91];
  static const double bnm_ch[5][91];
  static const double bnm_cw[5][91];

  /* degree and order of SH expansion (spatial dependence) */
  static constexpr const int MAX_DEGREE = 12;
  static constexpr const int MAX_ORDER = 12;

  /* if two epochs are less than MAX_SEC_DIFF away, do not recompute temporal
   * dependence. use the trigs as stored in [cs]d[24]pi
   */
  static constexpr const double MAX_SEC_DIFF = 60;

  /* last epoch used for computing temporal dependence, i.e. trigs based on
   * date of interest
   */
  MjdEpoch last_epoch_ = MjdEpoch::min();

  /* coefficients (i.e. trigs) for temporal dependence. computed at last_epoch
   */
  double sd2pi_;
  double cd2pi_;
  double sd4pi_;
  double cd4pi_;

  /* recompute temporal dependence, i.e. trigs [cs]d[24]pi */
  void date_trigs(const MjdEpoch &t) noexcept {
    const dso::ydoy_date yd(t.to_ydoy());
    const double doy =
        yd.dy().as_underlying_type() + t.fractional_days().days();
    sd2pi_ = std::sin(doy / 365.25 * 2 * dso::DPI);
    cd2pi_ = std::cos(doy / 365.25 * 2 * dso::DPI);
    sd4pi_ = std::sin(doy / 365.25 * 4 * dso::DPI);
    cd4pi_ = std::cos(doy / 365.25 * 4 * dso::DPI);
    last_epoch_ = t;
    return;
  }

  /* easy access for coefficient arrays */
  static const double *anm_bh_A0() noexcept { return anm_bh[1 - 1]; }
  static const double *anm_bh_A1() noexcept { return anm_bh[2 - 1]; }
  static const double *anm_bh_B1() noexcept { return anm_bh[3 - 1]; }
  static const double *anm_bh_A2() noexcept { return anm_bh[4 - 1]; }
  static const double *anm_bh_B2() noexcept { return anm_bh[5 - 1]; }
  static const double *anm_bw_A0() noexcept { return anm_bw[1 - 1]; }
  static const double *anm_bw_A1() noexcept { return anm_bw[2 - 1]; }
  static const double *anm_bw_B1() noexcept { return anm_bw[3 - 1]; }
  static const double *anm_bw_A2() noexcept { return anm_bw[4 - 1]; }
  static const double *anm_bw_B2() noexcept { return anm_bw[5 - 1]; }
  static const double *anm_ch_A0() noexcept { return anm_ch[1 - 1]; }
  static const double *anm_ch_A1() noexcept { return anm_ch[2 - 1]; }
  static const double *anm_ch_B1() noexcept { return anm_ch[3 - 1]; }
  static const double *anm_ch_A2() noexcept { return anm_ch[4 - 1]; }
  static const double *anm_ch_B2() noexcept { return anm_ch[5 - 1]; }
  static const double *anm_cw_A0() noexcept { return anm_cw[1 - 1]; }
  static const double *anm_cw_A1() noexcept { return anm_cw[2 - 1]; }
  static const double *anm_cw_B1() noexcept { return anm_cw[3 - 1]; }
  static const double *anm_cw_A2() noexcept { return anm_cw[4 - 1]; }
  static const double *anm_cw_B2() noexcept { return anm_cw[5 - 1]; }
  static const double *bnm_bh_A0() noexcept { return bnm_bh[1 - 1]; }
  static const double *bnm_bh_A1() noexcept { return bnm_bh[2 - 1]; }
  static const double *bnm_bh_B1() noexcept { return bnm_bh[3 - 1]; }
  static const double *bnm_bh_A2() noexcept { return bnm_bh[4 - 1]; }
  static const double *bnm_bh_B2() noexcept { return bnm_bh[5 - 1]; }
  static const double *bnm_bw_A0() noexcept { return bnm_bw[1 - 1]; }
  static const double *bnm_bw_A1() noexcept { return bnm_bw[2 - 1]; }
  static const double *bnm_bw_B1() noexcept { return bnm_bw[3 - 1]; }
  static const double *bnm_bw_A2() noexcept { return bnm_bw[4 - 1]; }
  static const double *bnm_bw_B2() noexcept { return bnm_bw[5 - 1]; }
  static const double *bnm_ch_A0() noexcept { return bnm_ch[1 - 1]; }
  static const double *bnm_ch_A1() noexcept { return bnm_ch[2 - 1]; }
  static const double *bnm_ch_B1() noexcept { return bnm_ch[3 - 1]; }
  static const double *bnm_ch_A2() noexcept { return bnm_ch[4 - 1]; }
  static const double *bnm_ch_B2() noexcept { return bnm_ch[5 - 1]; }
  static const double *bnm_cw_A0() noexcept { return bnm_cw[1 - 1]; }
  static const double *bnm_cw_A1() noexcept { return bnm_cw[2 - 1]; }
  static const double *bnm_cw_B1() noexcept { return bnm_cw[3 - 1]; }
  static const double *bnm_cw_A2() noexcept { return bnm_cw[4 - 1]; }
  static const double *bnm_cw_B2() noexcept { return bnm_cw[5 - 1]; }

 public:
  /** @brief This function will compute the spatial part of the b and c
   * 'empirical' VMF3 coefficients
   *
   * These coefficients depend only on the coordinates of the site, and do
   * not involve any temporal dependence. Hence, they can be computed once for
   * every site of interest.
   *
   * @param[in] rsta Cartesian ECEF coordinates of site [m]
   * @param[out] coeffs An instance of vmf3::Vmf3FullCoeffs which can be later
   * used to compute the VMF3 mapping functions via mf()
   */
  static int vmf3_spatial_coeffs(dso::CartesianCrdConstView &rsta,
                                 vmf3::Vmf3FullCoeffs &coeffs) noexcept;

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
  void mf(const MjdEpoch &t, double el, double ah, double aw,
          const vmf3::Vmf3FullCoeffs &bc, double &mfh, double &mfw) noexcept;

}; /* struct Vmf3 */

namespace vmf3 {
/* Contains a Site and grid data from interpolation (for some epoch). */
struct SiteBlock {
  /* site name +  domes */
  vmf3::Site msite;
  /* ellispoidal coordinates */
  dso::GeodeticCrd mcrd;
  /* grid data for surrounding nodes in the order: [bl, br, tl, tr] */
  vmf3::GridVmf3Data::Data mdata_t0[4];
  /* grid data for surrounding nodes in the order: [bl, br, tl, tr] */
  vmf3::GridVmf3Data::Data mdata_t1[4];
  /* coeffs for each surrounding node, [bl, br, tl, tr] */
  vmf3::Vmf3FullCoeffs msitebc[4];

  int compute_spatial_vmf3_coeffs() noexcept {
    dso::CartesianCrd rsta = dso::geodetic2cartesian(mcrd);
  Vmf3::vmf3_spatial_coeffs(CartesianCrdConstView(rsta), 
                                 vmf3::Vmf3FullCoeffs &coeffs) noexcept;
  }

  void cp_t0tot1() noexcept {
    for (int i = 0; i < 4; i++) {
      std::memcpy(mdata_t1[i].raw_ptr(), mdata_t0[i].raw_ptr(),
                  sizeof(double) * GridVmf3Data::Data::NUM_ENTRIES);
    }
  }
  void cp_t1tot0() noexcept {
    for (int i = 0; i < 4; i++) {
      std::memcpy(mdata_t0[i].raw_ptr(), mdata_t1[i].raw_ptr(),
                  sizeof(double) * GridVmf3Data::Data::NUM_ENTRIES);
    }
  }

  template <ellipsoid E, typename C = CartesianCrd>
  SiteBlock(const char *site, const C &crd) noexcept {
    static_assert(dso::CoordinateTypeTraits<C>::isCartesian);
    std::strcpy(msite.data, site);
    mcrd = dso::cartesian2geodetic<E>(crd);
    for (int i = 0; i < 4; i++) msitebc[i].set_missing();
  }

}; /* struct SiteBlock */

struct Vmf3Result {
  double zhd_, zwd_, mfh_, mfw_;
  double zhd() const noexcept { return zhd_; }
  double zwd() const noexcept { return zwd_; }
  double mfh() const noexcept { return mfh_; }
  double mfw() const noexcept { return mfw_; }
  double &zhd() noexcept { return zhd_; }
  double &zwd() noexcept { return zwd_; }
  double &mfh() noexcept { return mfh_; }
  double &mfw() noexcept { return mfw_; }
};
} /* namespace vmf3 */

class Vmf3SiteHandler {
  /* one file per 6 hours */
  static constexpr const int FILE_TIME_INTERVAL = 6;
  /* hour designator of last file for a given day (i.e. H18) */
  static constexpr const int MAX_FILE_HOUR = 24 - FILE_TIME_INTERVAL;

  /* t0 epoch, corresponding to mdata_t0 */
  MjdEpoch mt0 = MjdEpoch::max();
  /* t1 epoch, corresponding to mdata_t1 */
  MjdEpoch mt1 = MjdEpoch::min();
  /* site-specific data for each site of interest, for t0 and t1 */
  std::vector<vmf3::SiteBlock> msites{};
  /* data dir for locating grid files */
  std::string mdata_dir;
  /* a Vmf3 instance for computations */
  Vmf3 mvmf3;

  int compute_spatial_vmf3_coeffs() noexcept {
    for (auto &site : msites) {
    }
  }

  /* replace t0 with t1; t1 is left as is */
  void left_shift() noexcept {
    mt0 = mt1;
    for (auto &site : msites) {
      site.cp_t1tot0();
    }
  }
  void right_shift() noexcept {
    mt1 = mt0;
    for (auto &site : msites) {
      site.cp_t0tot1();
    }
  }

 public:
  Vmf3SiteHandler()
      : mt0(MjdEpoch::max()),
        mt1(MjdEpoch::min()),
        msites(),
        mdata_dir{},
        mvmf3() {};

  template <typename T>
  void append_site(const char *site, const T &crd) noexcept {}

  int initialize(const MjdEpoch &t) noexcept {
    /* find surounding epochs for t */
    const double hours = t.seconds().seconds() / 3600e0;
    const int intrv = (int)(hours / (double)FILE_TIME_INTERVAL);
    const int hp = intrv * FILE_TIME_INTERVAL;
    const int hn =
        ((hp + FILE_TIME_INTERVAL) >= 24) ? 0 : (hp + FILE_TIME_INTERVAL);
    if (hp < 0 || hn > 18) {
      fprintf(stderr,
              "[ERROR] Failed locating a time interval for given date "
              "(traceback: %s)\n",
              __func__);
      return 1;
    }

    /* load the file prior to t */
    if (load_sites_for_epoch(t.to_ymd(), hp)) {
      fprintf(
          stderr,
          "[ERROR] Failed loading VMF3(GR) grid data file (traceback: %s)\n",
          __func__);
      return 1;
    }
    /* good, but t0 is now loaded as t1; left shift the data */
    this->left_shift();

    /* load the file next to t */
    auto ymdn = t.to_ymd();
    if (hn == 0) ymdn = t.add_seconds(FractionalSeconds(3599e0)).to_ymd();
    if (load_sites_for_epoch(ymdn, hn)) {
      fprintf(
          stderr,
          "[ERROR] Failed loading VMF3(GR) grid data file (traceback: %s)\n",
          __func__);
      return 1;
    }
    /* good, t1 now holds the data immidiate after t */
    this->left_shift();

    // TODO load orography
    return 0;
  }

  int load_correct_interval(const MjdEpoch &t) noexcept {
    constexpr const double INTRV = FILE_TIME_INTERVAL * 3600e0;
    if (t >= mt0 && t < mt1) return 0;

    /* the most probable case would be that we overrun the interval [t0, t1]
     * and we now need [t1, t1+1]
     */
    if (mt1 <= t &&
        t.diff<dso::DateTimeDifferenceType::FractionalSeconds>(mt1).seconds() <
            INTRV) {
      /* first move t1 to t0 */
      this->left_shift();
      /* find date and hours of day for next file */
      auto tn = mt0.add_seconds(FractionalSeconds(INTRV));
      const double hours = tn.seconds().seconds() / 3600e0;
      const int intrv = (int)(hours / (double)FILE_TIME_INTERVAL);
      const int hn = intrv * FILE_TIME_INTERVAL;
      /* load the file prior to t */
      if (load_sites_for_epoch(tn.to_ymd(), hn)) {
        fprintf(
            stderr,
            "[ERROR] Failed loading VMF3(GR) grid data file (traceback: %s)\n",
            __func__);
        return 1;
      }
    }

    /* if this is not the case fuck it, re-initialize */
    if (initialize(t)) {
      fprintf(
          stderr,
          "[ERROR] Failed loading VMF3(GR) grid data file(s) (traceback: %s)\n",
          __func__);
      return 1;
    }

    // TODO load orography
    return 0;
  }

 private:
  /** @brief Load a VMF3GR data grid and populate instance's  msites.
   *
   * The function will try to find a VMF3GR grid file located within the
   * m_data_dir directory, using the date passed in. If it finds a
   * corresponding data file, it will open it and load the grid. Then, it
   * will loop through the individual sites stored in msites, and for each
   * site it will compute the surounding 4 nodes (i.e. the cell) and load
   * the data for each of these 4 nodes. The data are going to be stored
   * in the calling instance's msites vector. There are two places where
   * the parse data can be stored: either mdata_t0 or mdata_t1 (for each
   * site). The function will store the data in the latter, i.e. mdata_t1.
   *
   * Note that mdata_t1 also hold an orography height, which will be NOT
   * set by the function. You will need a different call to extract this,
   * using a corresponding orography file.
   *
   * On success, mt1 will be assigned the passed-in date (i.e. from ymd
   * and day_hours).
   *
   * @param[in] ymd The date of request
   * @param[in] day_hours The hour of day (of request). This (and the
   * above) parameters are used to form a possible filename for a VMF3GR
   * filename. Thus, the day_hours parameter should be a multiple of 3
   * with a 24-hour interval (i.e. 0, 3, 6, ...)
   * @return Anything other than zero denotes an error.
   */
  int load_sites_for_epoch(const ymd_date &ymd, int day_hours) noexcept;

  int load_sites_orography(const char *fn, std::size_t num_vals,
                           vmf3::GridVmf3Data *grid) noexcept;

  int vmf3(const char *site, const MjdEpoch &t, double el,
           vmf3::Vmf3Result &result) noexcept;
}; /*class Vmf3SiteHandler */

} /* namespace dso */
