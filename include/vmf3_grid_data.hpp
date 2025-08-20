#ifndef __DSO_VMF3_GRID_DATA_STREAM_HPP__
#define __DSO_VMF3_GRID_DATA_STREAM_HPP__

#include <cstdint>
#include <cstdlib>
#include <limits>

#include "datetime/calendar.hpp"
#include "geodesy/transformations.hpp"
#include "tickaxis_2d.hpp"

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

class GridVmf3Data {
 public:
  struct Data {
    //(1) 	latitude [°]
    //(2) 	longitude [°]
    //(3)  height from orography (optional)
    //(4) 	hydrostatic "a" coefficient
    //(5) 	wet "a" coefficient
    //(6) 	zenith hydrostatic delay [m]
    //(7) 	zenith wet delay [m]
    //(8) 	hydrostatic north gradient [mm]
    //(9) 	hydrostatic east gradient [mm]
    //(10) 	wet north gradient [mm]
    //(11) 	wet east gradient [mm]
    static constexpr const int data_start_at_ = 3;
    static constexpr const int NUM_DATA_ELEMENTS = 8;
    static constexpr const int NUM_ENTRIES = NUM_DATA_ELEMENTS + data_start_at_;

    double data_[NUM_ENTRIES];
    double *data() noexcept { return data_ + data_start_at_; }
    const double *data() const noexcept { return data_ + data_start_at_; }

    double lat_deg() const noexcept { return data_[0]; }
    double lon_deg() const noexcept { return data_[1]; }
    double oro_ell() const noexcept { return data_[2]; }
    double &lat_deg() noexcept { return data_[0]; }
    double &lon_deg() noexcept { return data_[1]; }
    double &oro_ell() noexcept { return data_[2]; }

    double ah() const noexcept { return data_[data_start_at_]; }
    double aw() const noexcept { return data_[data_start_at_ + 1]; }
    double zhd() const noexcept { return data_[data_start_at_ + 2]; }
    double zwd() const noexcept { return data_[data_start_at_ + 3]; }
  };

 private:
  TickAxis2D_NP axis_;
  Data *data_ = nullptr;
  std::size_t capacity_ = 0;

 public:
  using index_type = TickAxis2D_NP::index_type;

  const TickAxis2D_NP &axis() const noexcept { return axis_; }

  GridVmf3Data(double ostart, double ostop, double ostep, double istart,
               double istop, double istep)
      : axis_(ostart, ostop, ostep, istart, istop, istep, 1e-9, 0e0, 360e0),
        data_((Data *)std::malloc(axis_.num_pts() * sizeof(Data))),
        capacity_(axis_.num_pts()) {};
  GridVmf3Data() : axis_(), data_(nullptr), capacity_(0) {};
  ~GridVmf3Data() noexcept {
    std::free(data_);
    capacity_ = 0;
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
                 double istop, double istep, double eps_steps = 1e-12,
                 double inner_period_hint = 360e0) noexcept {
    axis_ = TickAxis2D_NP(ostart, ostop, ostep, istart, istop, istep, eps_steps,
                          0e0, inner_period_hint);
    index_type npts = axis_.num_pts();
    if ((std::size_t)npts > capacity_) {
      if (data_) std::free(data_);
      data_ = (Data *)std::malloc(axis_.num_pts() * sizeof(Data));
      capacity_ = axis_.num_pts();
    }
  }

  index_type num_pts() const noexcept { return axis_.num_pts(); }

  Data *data(std::size_t idx) noexcept { return data_ + idx; }
  const Data *data(std::size_t idx) const noexcept { return data_ + idx; }

  //[[no_discard]] int bilinear_interpolation(double lat_deg, double lon_deg,
  //                                          Data &out) const noexcept;
  // int bilinear_interpolation_nocheck(double lat_deg, double lon_deg,
  //                                   Data &out) const noexcept;

}; /* class GridVmf3Data */

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
