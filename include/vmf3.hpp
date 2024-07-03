#ifndef __DSO_TUWIEN_VMF_TROPO_PRODUCTS_HPP__
#define __DSO_TUWIEN_VMF_TROPO_PRODUCTS_HPP__

#include <vector>
#include <fstream>
#include <cstring>
#include "datetime/calendar.hpp"

namespace dso {

/** @class
 * A simple class to hold all element records of a VMF3 record (based on 
 * VMF3 products published by TU WIEN).
 */
class Vmf3Data {
  double data[7];
  dso::MjdEpoch mt;

public:
  /* get datetime */
  const auto &t() const noexcept {return mt;}
  auto &t() noexcept {return mt;}
/* (3) hydrostatic mf coefficient a_h */
  const double &ah() const noexcept {return data[0];}
  double &ah() noexcept {return data[0];}
/* (4) wet mf coefficient a_w */
  const double &aw() const noexcept {return data[1];}
  double &aw() noexcept {return data[1];}
/* (5) zenith hydrostatic delay (m) */
  const double &zhd() const noexcept {return data[2];}
  double &zhd() noexcept {return data[2];}
/* (6) zenith wet delay (m) */
  const double &zwd() const noexcept {return data[3];}
  double &zwd() noexcept {return data[3];}
/* (7) pressure at the site (hPa) */
  const double &pressure() const noexcept {return data[4];}
  double &pressure() noexcept {return data[4];}
/* (8) temperature at the site (C) */
  const double &temperature() const noexcept {return data[5];}
  double &temperature() noexcept {return data[5];}
/* (9) water vapour pressure at the site (hPa) */
  const double &water_vapour_pressure() const noexcept {return data[6];}
  double &water_vapour_pressure() noexcept {return data[6];}
}; /* class Vmf3Data */

class Vmf3FileStream {
  static constexpr const int LINE_SZ = 124;
  char mfn[256];
  std::ifstream mstream;
  /* end of header position within the file */
  std::ifstream::pos_type meoh{-1};
  /* last line buffered */
  char bline[LINE_SZ] = "\0";

public:
  Vmf3FileStream(const char *fn) : mfn(std::strcpy(mfn, fn)), mstream(mfn) {}
  /* no copy */
  Vmf3FileStream(const Vmf3FileStream &) = delete;
  /* no assign operator */
  Vmf3FileStream& operator=(const Vmf3FileStream &) = delete;

}; /* class Vmf3FileStream */

class Vmf3SiteStream {
  std::vector<const char *> msites;
  std::vector<Vmf3Data> mdata;
  Vmf3FileStream mstream;
  
  void swap_epochs() noexcept {
    for (auto it = mdata.begin(); it != mdata.end(); it+=2)
      std::swap(*it, *(it+1));
    return;
  }

public:

}; /* class Vmf3SiteStream */

} /* namespace dso */

#endif
