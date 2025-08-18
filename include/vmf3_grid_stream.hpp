#ifndef __DSO_VMF3_GRID_DATA_STREAM_HPP__
#define __DSO_VMF3_GRID_DATA_STREAM_HPP__

#include "datetime/calendar.hpp"
#include <cstdlib>

namespace dso
{
  namespace vmf3
  {
    /** @brief 4-character id (optionally including domes) */
    struct Site
    {
      static constexpr const int name_start_at = 0;
      static constexpr const int domes_start_at = 5;

      char data[5 + 10] = {'\0'};
      const char *name() const { return data + name_start_at; }
      const char *dones() const { return data + domes_start_at; }
    };

    /* Contains a Site and grid data from interpolation (for some epoch). */
    struct SiteBlock
    {
      static constexpr const int NUM_DATA_FIELDS = 11;
      Site msite;
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
    }; /* struct SiteBlock */

    /** @brief A whole block of site-specific VMF3 parameters for a given epoch.
     */
    struct Block
    {
      /* reference epoch */
      MjdEpoch t_;
      /* site-specific records */
      std::vector<SiteBlock> site_data_;

      void append(const SiteBlock &b) noexcept { site_data_.emplace_back(b); }

      void reset()
      {
        site_data_.clear();
        t_ = dso::MjdEpoch::min();
      }
    }; /* struct Block */

    /** Full set of coefficients for computing b and c 'empirical' VMF3
     * coefficients. These coefficients have only a spatial dependence, hence can
     * be computed once for every site of interest.
     *
     * To fully compute the coefficients, temporal dependence must be considered,
     * which is done in the function Vmf3::compute().
     */
    struct Vmf3FullCoeffs
    {
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

    class TickAxis
    {
      double start_;
      double step_;
      int num_pts_;

    public:
      /** C'tor; note that limits are inclusive, i.e. both start and stop are valid ticks on the axis, which contains the range [start, stop] */
      TickAxis(double start, double stop, double step) : start_(start), step_(step), num_pts_((stop - start) / step + 1)
      {
        if (num_pts_ < 0)
          throw std::runtime_error("[ERROR] Failed creating TickAxis from specified marks!\n");
      };
      TickAxis() : start_(0), step_(1), num_pts_(0) {};
      double tick_start() const noexcept { return start_; }
      double tick_step() const noexcept { return step_; }
      double tick_stop() const noexcept { return start_ + step_ * num_pts_; }

      int num_pts() const noexcept { return num_pts_; }
      int index(double val) const noexcept { return (val - start_) / step_; }
      double val(int index) const noexcept { return start_ + step_ * index; }
    }; /* class TickAxis */

    class TickAxis2D
    {
      TickAxis outter_;
      TickAxis inner_;

    public:
      struct Cell
      {
        int bl, br, ul, ur;
      };

      std::size_t num_pts() const noexcept { return outter_.num_pts() * inner_.num_pts(); }
      const TickAxis &tick_axis_outter() const noexcept { return outter_; }
      const TickAxis &tick_axis_inner() const noexcept { return inner_; }

      TickAxis2D(double ostart, double ostop, double ostep, double istart, double istop, double istep) noexcept : outter_(ostart, ostop, ostep), inner_(istart, istop, istep) {};
      TickAxis2D() noexcept : outter_(), inner_() {};
      int index(double val_out, double val_in) const noexcept
      {
        return outter_.index(val_out) * inner_.num_pts() + inner_.index(val_in);
      }
      Cell cell(double val_out, double val_in) const noexcept
      {
        int bl = index(val_out, val_in);
        return Cell{bl, bl + inner_.num_pts(), bl + 1, bl};
      }
    }; /* class TickAxis2D */

    class GridVmf3Data
    {
      struct Data
      {
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
        static constexpr const int NUM_ELEMENTS = 10;
        double data_[NUM_ELEMENTS];
        double *data() noexcept { return data_; }
      };

      TickAxis2D axis_;
      Data *data_ = nullptr;
      std::size_t capacity_ = 0;

    public:
      GridVmf3Data(double ostart, double ostop, double ostep, double istart, double istop, double istep) : axis_(ostart, ostop, ostep, istart, istop, istep),
                                                                                                           data_((Data *)std::malloc(axis_.num_pts() * sizeof(Data))),
                                                                                                           capacity_(axis_.num_pts()) {};
      GridVmf3Data() : axis_(),
                       data_(nullptr),
                       capacity_(0) {};
      ~GridVmf3Data() noexcept
      {
        std::free(data_);
        capacity_ = 0;
      }

      const TickAxis &tick_axis_outter() const noexcept
      {
        return axis_.tick_axis_outter();
      }
      const TickAxis &tick_axis_inner() const noexcept
      {
        return axis_.tick_axis_inner();
      }

      /* no copy allowed */
      GridVmf3Data(const GridVmf3Data &) = delete;

      /* move c'tor */
      GridVmf3Data(GridVmf3Data &&other) noexcept : axis_(other.axis_), data_(other.data_), capacity_(other.capacity_)
      {
        other.capacity_ = 0;
        other.data_ = nullptr;
      }

      /* no assignment operator */
      GridVmf3Data &operator=(const GridVmf3Data &) = delete;

      /* move assignment operator */
      GridVmf3Data &operator=(GridVmf3Data &&other) noexcept
      {
        axis_ = other.axis_;
        data_ = other.data_;
        capacity_ = other.capacity_;
        other.data_ = nullptr;
        other.capacity_ = 0;
        return *this;
      }

      void set_range(double ostart, double ostop, double ostep, double istart, double istop, double istep) noexcept
      {
        axis_ = TickAxis2D(ostart, ostop, ostep, istart, istop, istep);
        std::size_t npts = axis_.num_pts();
        if (npts > capacity_)
        {
          if (data_)
            std::free(data_);
          data_ = (Data *)std::malloc(axis_.num_pts() * sizeof(Data));
          capacity_ = axis_.num_pts();
        }
      }

      std::size_t num_pts() const noexcept { return axis_.num_pts(); }

      TickAxis2D::Cell cell(double lat, double lon) const noexcept
      {
        return axis_.cell(lat, lon);
      }

      Data *data(std::size_t idx) noexcept { return data_ + idx; }
    }; /* class GridData */

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