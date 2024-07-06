#ifndef __DSO_TUWIEN_VMF_TROPO_PRODUCTS_HPP__
#define __DSO_TUWIEN_VMF_TROPO_PRODUCTS_HPP__

#include <vector>
#include <fstream>
#include <cstring>
#include "datetime/calendar.hpp"

namespace dso {

namespace vmf3_details {
std::vector<const char *>::const_iterator
find_if_sorted_string(const char *site,
                      const std::vector<const char *> &sites) noexcept;
}/* namespace vmf3_details */

/** @class
 * A simple class to hold all element records of a VMF3 record (based on 
 * VMF3 products published by TU WIEN).
 * The data correspond to a given site, at a given epoch.
 * The data can be found at:
 * https://vmf.geo.tuwien.ac.at/trop_products/DORIS/VMF3/VMF3_OP/
 */
class Vmf3SiteData {
  MjdEpoch mt;
  double mdata[7];

public:
  /* get datetime */
  const auto &t() const noexcept {return mt;}
  auto &t() noexcept {return mt;}
/* (3) hydrostatic mf coefficient a_h */
  const double &ah() const noexcept {return mdata[0];}
  double &ah() noexcept {return mdata[0];}
/* (4) wet mf coefficient a_w */
  const double &aw() const noexcept {return mdata[1];}
  double &aw() noexcept {return mdata[1];}
/* (5) zenith hydrostatic delay (m) */
  const double &zhd() const noexcept {return mdata[2];}
  double &zhd() noexcept {return mdata[2];}
/* (6) zenith wet delay (m) */
  const double &zwd() const noexcept {return mdata[3];}
  double &zwd() noexcept {return mdata[3];}
/* (7) pressure at the site (hPa) */
  const double &pressure() const noexcept {return mdata[4];}
  double &pressure() noexcept {return mdata[4];}
/* (8) temperature at the site (C) */
  const double &temperature() const noexcept {return mdata[5];}
  double &temperature() noexcept {return mdata[5];}
/* (9) water vapour pressure at the site (hPa) */
  const double &water_vapour_pressure() const noexcept {return mdata[6];}
  double &water_vapour_pressure() noexcept {return mdata[6];}

  double* data() noexcept {return mdata;}

  Vmf3SiteData(MjdEpoch t = MjdEpoch::max(), double *data = nullptr)
      : mt(t) {
    if (data)
      std::memcpy(mdata, data, sizeof(double) * 7);
    else
      std::memset((void *)mdata, 0, sizeof(double) * 7);
  }
}; /* class Vmf3SiteData */

class Vmf3SiteFileStream {
  static constexpr const int LINE_SZ = 124;
  std::ifstream mstream;
  /* end of header position within the file */
  // std::ifstream::pos_type meoh{-1};
  /* last line buffered */
  char bline[LINE_SZ] = "#\0";
  char mfn[256];
  MjdEpoch mcurrent{MjdEpoch::max()};
  MjdEpoch mnext{MjdEpoch::min()};

public:
  const char *fn() const noexcept {return mfn;}

  Vmf3SiteFileStream(const char *fn) noexcept : mstream(fn) {
    std::strcpy(mfn, fn);
  }
  /* no copy */
  Vmf3SiteFileStream(const Vmf3SiteFileStream &) = delete;
  /* no assign operator */
  Vmf3SiteFileStream& operator=(const Vmf3SiteFileStream &) = delete;

  int skip_block() noexcept;
  int parse_block(const std::vector<const char *> &sites,
                  std::vector<Vmf3SiteData> &data, int recs_per_site = 1,
                  int rec_offset = 0) noexcept;
  /* return the epoch in the buffered line */
  // MjdEpoch buffered_epoch() const ;
  /* stream state according to std::basic_ios<CharT,Traits>::good */
  auto good() const noexcept {return mstream.good();}
  /* set stream to top of file */
  void rewind() noexcept {mstream.seekg(0);}
  /* Current (just skipped/parsed) epoch */
  const MjdEpoch &current_epoch() const noexcept { return mcurrent; }
  /*  Next epoch (to be skipped/parsed); first line already buffered. */
  const MjdEpoch &next_epoch() const noexcept { return mnext; }

}; /* class Vmf3SiteFileStream */

class Vmf3SiteStream {
  static constexpr const int RECS_PER_SITE = 2;
  static constexpr const int MAX_SITE_CHARS = 4;
  static constexpr const int SITE_LEN = MAX_SITE_CHARS + 1;

  std::vector<const char *> msites;
  std::vector<Vmf3SiteData> mdata;
  /* current interval */
  MjdEpoch t0{MjdEpoch::max()};
  MjdEpoch t1{MjdEpoch::min()};
  Vmf3SiteFileStream mstream;
  char *mmemsites{nullptr};

  /* @brief Initialize the msites and mdata vectors given a list of sites.
   *
   * The given sites (within the sites vector) need not be null-terminated 
   * strings, but the strings should start with the site name and hold the 
   * name with  MAX_SITE_CHARS characters.
   * E.g. the following is correct:
   * sites[1] = "DIOA" or "DIOA  ", or "DIOA\0" or "DIOAAAAAAA" ...
   * but the following is incorrect:
   * sites[1] = " DIOA"
   *
   * The sites will be copyied to the instances' msites vector in 
   * lexicographical order. Any previous data within the msites vector will 
   * be erased/lost.
   * 
   * The mdata vector will be initialized with the correct number of elements 
   * (corresponding to the number of siztes given multiplied by RECS_PER_SITE).
   *
   * Note that the function will also rewind the stream to the top-of-the-file
   * position (aka we will be reading from the top of the data stream).
   *
   * This function should be used to 'initialize' the instance.
   *
   * @param[in] sites List of sites (see description).
   */
  void initialize(const std::vector<const char *> &sites) noexcept;
  
  int forward_search(const MjdEpoch &t) noexcept;

public:
  Vmf3SiteStream(const char *fn, const std::vector<const char *> &sites);
  ~Vmf3SiteStream() noexcept {if (mmemsites) delete[] mmemsites;}

  int set_at_epoch(const MjdEpoch &t) noexcept;
  int operator()(const MjdEpoch &t) noexcept;

  const std::vector<const char *> &sites() const noexcept {return msites;}

  void swap_epochs() noexcept {
    for (auto it = mdata.begin(); it != mdata.end(); it+=2)
      std::swap(*it, *(it+1));
    std::swap(t0,t1);
    return;
  }

  const MjdEpoch interval_start() const noexcept {return t0;}
  const MjdEpoch interval_stop() const noexcept {return t1;}

int site_vmf3(const char *site, const dso::MjdEpoch &t,
                                   dso::Vmf3SiteData &vmf3) noexcept;


}; /* class Vmf3SiteStream */

} /* namespace dso */

#endif
