#ifndef __DSO_VMF3_GRID_DATA_STREAM_HPP__
#define __DSO_VMF3_GRID_DATA_STREAM_HPP__

#include <cstdint>
#include <cstdlib>
#include <limits>

#include "datetime/calendar.hpp"
#include "geodesy/transformations.hpp"

namespace dso {
namespace vmf3 {
/** @brief 4-character id (optionally including domes) */
struct Site {
  static constexpr const int name_start_at = 0;
  static constexpr const int domes_start_at = 5;

  char data[5 + 10] = {'\0'};
  const char *name() const { return data + name_start_at; }
  const char *dones() const { return data + domes_start_at; }
};

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

/** @brief Discrete, inclusive 1-D axis with robust floating-point handling.
 *
 * The TickAxis class models a one-dimensional, regularly spaced grid of
 * tick marks. The axis is **inclusive**: both @c start and @c stop are
 * valid ticks, so the number of points is @f$N = \mathrm{round}((stop -
 * start)/step) + 1@f$, validated within a user-specified tolerance.
 *
 * The constructor normalizes the sign of @c step to match the geometric
 * direction from @c start to @c stop (ascending or descending). A divisibility
 * check ensures the span is (within tolerance) an integer multiple of the step.
 *
 * Designed for geographic grids (e.g., latitude [89.5,-89.5], step=-1),
 * but general enough for any numeric axis.
 */
class TickAxis {
  double start_; /**< @brief First tick value on the axis (inclusive). */
  double step_;  /**< @brief Signed spacing between consecutive ticks. */
  int num_pts_;  /**< @brief Number of ticks on the axis (inclusive count). */

 public:
  /** @brief Construct an inclusive axis from @p start to @p stop with
   * spacing @p step.
   *
   * The constructor:
   *   1. Validates inputs are finite and @p step ≠ 0.
   *   2. Sets the internal step sign to match the axis direction:
   *      - If @p stop >= @p start, @c tick_step() > 0 (ascending).
   *      - Else, @c tick_step() < 0 (descending).
   *   3. Computes the theoretical number of steps
   *      @f$ raw = (stop - start)/step @f$, rounds to nearest integer,
   *      and validates that the span is divisible within a tolerance
   *      @f$ tol = |step|\cdot eps\_steps @f$.
   *   4. Sets the inclusive point count to @f$ n+1 @f$.
   *
   * @param start      First tick value (inclusive).
   * @param stop       Last tick value (inclusive).
   * @param step       Nominal step size (magnitude used; sign normalized
   * internally).
   * @param eps_steps  Tolerance as a fraction of @c |step| used to validate
   *                   divisibility (default: 1e-9).
   *
   * @throws std::runtime_error
   *   - If any of @p start, @p stop, @p step are non-finite.
   *   - If @p step == 0.
   *   - If @p eps_steps < 0.
   *   - If (stop - start) is not (within tolerance) an integer multiple of @p
   * step.
   *   - If the resulting number of points is non-positive or exceeds @c
   * INT_MAX.
   *
   * @par Examples
   * @code
   * // Ascending:
   * TickAxis x(0.0, 10.0, 1.0);   // step becomes +1.0, num_pts = 11
   *
   * // Descending (latitude-like):
   * TickAxis lat(89.5, -89.5, 1.0); // step becomes -1.0, num_pts = 180
   * @endcode
   */
  TickAxis(double start, double stop, double step, double eps_steps = 1e-9);

  /** @brief Default constructor: empty axis.
   *
   * Initializes an empty axis with @c start_=0, @c step_=1, @c num_pts_=0.
   * @note An empty axis has no valid ticks.
   */
  TickAxis() : start_(0), step_(1), num_pts_(0) {};
  /**
   * @brief Get the first tick value (inclusive).
   * @return The start tick value.
   */
  double tick_start() const noexcept { return start_; }

  /** @brief Get the signed spacing between ticks.
   * @return The step size; positive for ascending axes, negative for
   * descending.
   */
  double tick_step() const noexcept { return step_; }

  /** @brief Get the last tick value (inclusive).
   * @return @c tick_start() + tick_step() * (num_pts() - 1).
   */
  double tick_stop() const noexcept { return start_ + step_ * (num_pts_ - 1); }

  /** @brief Get the number of ticks on the axis (inclusive count).
   * @return The number of points @f$N \ge 0@f$.
   */
  int num_pts() const noexcept { return num_pts_; }

  /** @brief Map a coordinate value to the index of the **bottom/left** tick
   * below it.
   *
   * Computes the normalized coordinate @f$ u = (val - start)/step @f$ and
   * returns:
   *   - @c floor(u - eps) for ascending axes (step > 0),
   *   - @c ceil(u + eps)  for descending axes (step < 0),
   * where @c eps = 1e-12 biases exact-on-grid values to the lower cell to
   * avoid ambiguity at tick boundaries.
   *
   * This is useful when you need the bottom-left vertex of the cell
   * surrounding a point, i.e., the largest tick not greater than @p val
   * in geometric terms.
   *
   * @param val  Coordinate to index.
   * @return Integer index corresponding to the bottom/left tick.
   *
   * @note No bounds checking is performed; callers should clamp to
   *       [0, num_pts()-2] when they require a valid cell (for i+1 access).
   * @see val(int)
   */
  int index(double val) const noexcept;

  /** @brief Map an integer index to the corresponding tick value.
   * @param index  Zero-based tick index.
   * @return @c tick_start() + index * tick_step().
   *
   * @note No bounds checking is performed; valid indices are in [0,
   * num_pts()-1].
   */
  double val(int index) const noexcept { return start_ + step_ * index; }

  bool in_range(double val) const noexcept {
    const double lo = std::min(tick_start(), tick_stop());
    const double hi = std::max(tick_start(), tick_stop());
    return val >= lo && val < hi;
  }

  /* Helper: does v lie between nodes k and k+1 (lower inclusive, upper
   * inclusive within tol)?
   */
  bool within_segment_inclusive(double val, int idx0, int idx1) const noexcept {
    const double a0 = this->val(idx0), a1 = this->val(idx1);
    const double lo = std::min(a0, a1), hi = std::max(a0, a1);
    const double tol = std::abs(tick_step()) * 1e-12;
    return (val >= lo - tol) && (val <= hi + tol);
  }

}; /* class TickAxis */

class TickAxis2D {
  TickAxis outter_;
  TickAxis inner_;

  /* helper: get (i,j) without checks. */
  inline std::pair<int, int> index_ij_unchecked(double val_out,
                                                double val_in) const noexcept {
    return {outter_.index(val_out), inner_.index(val_in)};
  }

 public:
  using index_type = int64_t;
  
  struct Cell {
    int bl, br, tl, tr;
  };

  index_type num_pts() const noexcept {
    return outter_.num_pts() * inner_.num_pts();
  }
  const TickAxis &tick_axis_outter() const noexcept { return outter_; }
  const TickAxis &tick_axis_inner() const noexcept { return inner_; }

  TickAxis2D(double ostart, double ostop, double ostep, double istart,
             double istop, double istep) noexcept
      : outter_(ostart, ostop, ostep), inner_(istart, istop, istep) {};

  TickAxis2D() noexcept : outter_(), inner_() {};

  index_type index(double val_out, double val_in) const noexcept;

  Cell cell(double val_out, double val_in) const;
  Cell cell_nocheck(double val_out, double val_in) const noexcept;

   /** @brief Map a row-major linear index to its tick values on the outer and inner axes.
   *
   * Uses the same linearization as TickAxis2D::index():
   * \f$ \text{idx} = i \cdot N_{\text{inner}} + j \f$ where @c N_inner = inner_.num_pts().
   *
   * @param idx Linear index (no bounds checking; may be negative).
   * @return std::pair<double,double> { outter_.val(i), inner_.val(j) }.
   *
   * @warning No bounds checking is performed. If @c inner_.num_pts() == 0 the division
   *          is undefined. Ensure the grid is initialized. Negative @p idx will
   *          produce negative @c i/@c j by truncating division toward zero.
   */
  std::pair<double,double> val(index_type idx) const noexcept {
    const index_type ncols = static_cast<index_type>(inner_.num_pts());
    const index_type i = idx / ncols; // row
    const index_type j = idx % ncols; // col
    return { outter_.val(static_cast<int>(i)),
             inner_.val (static_cast<int>(j)) };
  }
}; /* class TickAxis2D */

class GridVmf3Data {
 public:
  using index_type = TickAxis2D::index_type;

  struct Data {
//(1) 	latitude [°]
//(2) 	longitude [°]
//(3) 	hydrostatic "a" coefficient
//(4) 	wet "a" coefficient
//(5) 	zenith hydrostatic delay [m]
//(6) 	zenith wet delay [m]
//(7) 	hydrostatic north gradient [mm]
//(8) 	hydrostatic east gradient [mm]
//(9) 	wet north gradient [mm]
//(10) 	wet east gradient [mm]
#ifdef DEBUG
    static constexpr const int data_start_at_ = 2;
    static constexpr const int NUM_ELEMENTS = 10;
    double lat() const noexcept { return data_[0]; }
    double lon() const noexcept { return data_[1]; }
    double &lat() noexcept { return data_[0]; }
    double &lon() noexcept { return data_[1]; }
#else
    static constexpr const int data_start_at_ = 0;
    static constexpr const int NUM_ELEMENTS = 8;
#endif
    double data_[NUM_ELEMENTS];
    double *data() noexcept { return data_ + data_start_at_; }
    const double *data() const noexcept { return data_ + data_start_at_; }
  };

 private:
  TickAxis2D axis_;
  Data *data_ = nullptr;
  std::size_t capacity_ = 0;

 public:
  GridVmf3Data(double ostart, double ostop, double ostep, double istart,
               double istop, double istep)
      : axis_(ostart, ostop, ostep, istart, istop, istep),
        data_((Data *)std::malloc(axis_.num_pts() * sizeof(Data))),
        capacity_(axis_.num_pts()) {};
  GridVmf3Data() : axis_(), data_(nullptr), capacity_(0) {};
  ~GridVmf3Data() noexcept {
    std::free(data_);
    capacity_ = 0;
  }

  const TickAxis &tick_axis_outter() const noexcept {
    return axis_.tick_axis_outter();
  }
  const TickAxis &tick_axis_inner() const noexcept {
    return axis_.tick_axis_inner();
  }

  /* no copy allowed */
  GridVmf3Data(const GridVmf3Data &) = delete;

  /* move c'tor */
  GridVmf3Data(GridVmf3Data &&other) noexcept
      : axis_(other.axis_), data_(other.data_), capacity_(other.capacity_) {
    other.capacity_ = 0;
    other.data_ = nullptr;
  }

  /* no assignment operator */
  GridVmf3Data &operator=(const GridVmf3Data &) = delete;

  /* move assignment operator */
  GridVmf3Data &operator=(GridVmf3Data &&other) noexcept {
    axis_ = other.axis_;
    data_ = other.data_;
    capacity_ = other.capacity_;
    other.data_ = nullptr;
    other.capacity_ = 0;
    return *this;
  }

  void set_range(double ostart, double ostop, double ostep, double istart,
                 double istop, double istep) noexcept {
    axis_ = TickAxis2D(ostart, ostop, ostep, istart, istop, istep);
    std::size_t npts = axis_.num_pts();
    if (npts > capacity_) {
      if (data_) std::free(data_);
      data_ = (Data *)std::malloc(axis_.num_pts() * sizeof(Data));
      capacity_ = axis_.num_pts();
    }
  }

  index_type num_pts() const noexcept { return axis_.num_pts(); }

  TickAxis2D::Cell cell(double lat, double lon) const {
    return axis_.cell(lat, lon);
  }

  TickAxis2D::Cell cell_nocheck(double lat, double lon) const noexcept {
    return axis_.cell_nocheck(lat, lon);
  }
  
  std::pair<double,double> val(TickAxis2D::index_type idx) const noexcept {
    return axis_.val(idx);
  }
  
  Data *data(std::size_t idx) noexcept { return data_ + idx; }
  const Data *data(std::size_t idx) const noexcept { return data_ + idx; }

  int bilinear_interpolation(double lat_deg, double lon_deg, Data &out) const noexcept {
    try {
      const auto cell = this->cell(lat_deg, lon_deg);
      auto [x1, y1] = this->val(cell.bl);
      auto [x2, y2] = this->val(cell.tr);
      const double  x = lat_deg;
      const double  y = lon_deg;
      for (int i=0; i<Data::NUM_ELEMENTS; i++) {
        const double f11 = this->data(cell.bl)->data()[i];
        const double f21 = this->data(cell.br)->data()[i];
        const double f12 = this->data(cell.tl)->data()[i];
        const double f22 = this->data(cell.tr)->data()[i];
        const double fxy1 = (x2-x)/(x2-x1)*f11 + (x-x1)/(x2-x1)*f21;
        const double fxy2 = (x2-x)/(x2-x1)*f12 + (x-x1)/(x2-x1)*f22;
        out.data()[i] = (y2-y)/(y2-y1)*fxy1 + (y-y1)/(y2-y1)*fxy2;
      }
    } catch (std::exception &e) {
      return 1;
    }

    return 0;
  }

  index_type bl(double lat, double lon) const noexcept {
    return axis_.index(lat, lon);
  }
}; /* class GridVmf3Data */

/* Contains a Site and grid data from interpolation (for some epoch). */
struct SiteBlock {
  Site msite;
  dso::GeodeticCrd mcrd;
  GridVmf3Data::Data mdata_t0;
  GridVmf3Data::Data mdata_t1;
  Vmf3FullCoeffs msitebc;

  void linear_interpolation(const MjdEpoch &t) const noexcept;
}; /* struct SiteBlock */

} /* namespace vmf3 */

/** @brief Load a VMF3 grid data file into a  vmf3::GridVmf3Data.
 *
 * The grid file can be of any resolution (i.e. 1x1 or 5x5). The function
 * will first read the header and set the domensions of grid accordingly.
 * Hence, the grid instance can be anything at input, it is irrelevant.
 *
 * The grid data file should be of the 'V3GR' series, i.e. it is expected to
 * contain 11 columns of data (see https://vmf.geo.tuwien.ac.at/products.html)
 *
 * @param[in] fn    The filename of the grid file
 * @param[out] grid Instance holding the grid file at output
 * @return Anythin other than zero denotes an error.
 *
 */
int load_vmfgr3_grid_map(const char *fn, vmf3::GridVmf3Data *grid) noexcept;
} /* namespace dso */

#endif
