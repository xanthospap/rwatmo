#ifndef __DSO_TUWIEN_VMF_TROPO_PRODUCTS_HPP__
#define __DSO_TUWIEN_VMF_TROPO_PRODUCTS_HPP__

#include "datetime/calendar.hpp"
#include "geodesy/geodesy.hpp"
#include <cstring>
#include <fstream>
#include <vector>

namespace dso {

namespace vmf3_details {
/** Iterator to the first element matching the given string, assuming a sorted
 * vector.
 *
 * @param[in] site  Null-terminated string to match (normally a 4-char id)
 * @param[in] sites A (lexicogrtaphically) sorted vector of null-terminated
 *                  site names to match against.
 * @return Iterator to the first vector element matching the given string, or
 *         sites.end() if no match.
 */
std::vector<const char *>::const_iterator
find_if_sorted_string(const char *site,
                      const std::vector<const char *> &sites) noexcept;

/** b and c 'empirical' VMF3 coefficients (wet and dry component) */
struct Vmf3EmpiricalCoeffs {
  double bh, bw, ch, cw;

  /** VMF3 wet and dry maping functions given the 'a' coefficients
   *
   * The instance already holds the 'empirical' b and c coefficients. Given
   * the 'a' coefficients (e.g. from site-specific or grid files) and the
   * elevationa angle, the function will compute the wet and dry maping
   * functions.
   *
   * @param[in] el Elevation [rad]
   * @param[in] ah Hydrostatic 'a' coefficient
   * @param[in] aw Wet 'a' coefficient
   * @param[out] mfh Mapping function, hydrostatic component [m]
   * @param[out] mfw Mapping function, wet component [m]
   * @return Always 0
   */
  int mf(double el, double ah, double aw, double &mfh,
         double &mfw) const noexcept {
    const double sel = std::sin(el);
    /* compute the hydrostatic and wet mapping factors */
    mfh = (1 + (ah / (1 + bh / (1 + ch)))) /
          (sel + (ah / (sel + bh / (sel + ch))));
    mfw = (1 + (aw / (1 + bw / (1 + cw)))) /
          (sel + (aw / (sel + bw / (sel + cw))));
    return 0;
  }
}; /* Vmf3EmpiricalCoeffs */

/** Full set of coefficients for computing b and c 'empirical' VMF3
 * coefficients, given an epoch.
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

  /** Given an epoch, compute the b and c 'empirical' VMF3 coefficients */
  Vmf3EmpiricalCoeffs computeCoeffs(const dso::MjdEpoch &t) const noexcept;
}; /* Vmf3FullCoeffs */

/** Given longitude and latitude of a site, compute the Vmf3FullCoeffs */
Vmf3FullCoeffs vmf3FullCoeffSet(double lon, double lat) noexcept;
} /* namespace vmf3_details */

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
  const auto &t() const noexcept { return mt; }
  auto &t() noexcept { return mt; }
  /* (3) hydrostatic mf ('a') coefficient a_h */
  const double &ah() const noexcept { return mdata[0]; }
  double &ah() noexcept { return mdata[0]; }
  /* (4) wet mf ('a') coefficient a_w */
  const double &aw() const noexcept { return mdata[1]; }
  double &aw() noexcept { return mdata[1]; }
  /* (5) zenith hydrostatic delay [m] */
  const double &zhd() const noexcept { return mdata[2]; }
  double &zhd() noexcept { return mdata[2]; }
  /* (6) zenith wet delay [m] */
  const double &zwd() const noexcept { return mdata[3]; }
  double &zwd() noexcept { return mdata[3]; }
  /* (7) pressure at the site [hPa] */
  const double &pressure() const noexcept { return mdata[4]; }
  double &pressure() noexcept { return mdata[4]; }
  /* (8) temperature at the site [C] */
  const double &temperature() const noexcept { return mdata[5]; }
  double &temperature() noexcept { return mdata[5]; }
  /* (9) water vapour pressure at the site [hPa] */
  const double &water_vapour_pressure() const noexcept { return mdata[6]; }
  double &water_vapour_pressure() noexcept { return mdata[6]; }
  /* return pointer to internal storage */
  double *data() noexcept { return mdata; }

  Vmf3SiteData(MjdEpoch t = MjdEpoch::max(), double *data = nullptr) : mt(t) {
    if (data)
      std::memcpy(mdata, data, sizeof(double) * 7);
    else
      std::memset((void *)mdata, 0, sizeof(double) * 7);
  }
}; /* class Vmf3SiteData */

/** A site-specific VMF3 data file stream. Note that this instance is meant to
 * assist parsing a data file in a forward sense, i.e. from one epoch to the
 * next going forward into the file.
 */
class Vmf3SiteFileStream {
  /* max number of characters in a line */
  static constexpr const int LINE_SZ = 124;
  /* the file sttream */
  std::ifstream mstream;
  /* last line buffered */
  char bline[LINE_SZ] = "#\0";
  /* the filename */
  char mfn[256];
  /* last epoch read */
  MjdEpoch mcurrent{MjdEpoch::max()};
  /* epoch of new block; first line already stored in bline */
  MjdEpoch mnext{MjdEpoch::min()};

public:
  /* return the filename */
  const char *fn() const noexcept { return mfn; }

  /* constructor; only the filename is set, the stream is NOT opened */
  Vmf3SiteFileStream(const char *fn) noexcept : mstream(fn) {
    std::strcpy(mfn, fn);
  }

  /* no copy */
  Vmf3SiteFileStream(const Vmf3SiteFileStream &) = delete;

  /* no assign operator */
  Vmf3SiteFileStream &operator=(const Vmf3SiteFileStream &) = delete;

  /* Skip the next block.
   *
   * A 'block' is a series of data records of a common epoch. The first line
   * of this block is already stored in bline. Keep reading and skipping lines
   * untill we meet a data record of a later epoch.
   * When this happens, the last line (i.e. first line of next block) is
   * stored in bline and the function returns, having set the current and next
   * epochs.
   *
   * @return Anything other than 0 denotes an error.
   */
  int skip_block() noexcept;

  /* Parse the next block.
   *
   * A 'block' is a series of data records of a common epoch. The first line
   * of this block is already stored in bline. Keep reading and parsing lines
   * untill we meet a data record of a later epoch.
   * When this happens, the last line (i.e. first line of next block) is
   * stored in bline and the function returns, having set the current and next
   * epochs.
   *
   * Data lines that hold records for sites of interest are stored in the
   * data vector.
   *
   * @param[in] site A vector of sites (for DORIS 4-chars); each site name
   *            should be a null-terminated C-string
   * @param[out] data The vector where the data for the sites of interest will
   *            be stored. Data (i.e. VMF3 records) are copyied into the vector,
   *            NOT APPENDED, hence the vector should have the right size at
   *            input.
   *            If a given record line matches the site name at index j of the
   *            sites (input) vector, then the corresponding VMF3 (parsed) data
   *            will be stored at index j*recs_per_size+rec_offset.
   *            This "indexing" is used to allow multiple records of different
   *            datetimes to be stored in the data vector.
   * @param[in] recs_per_site Number of records per site in the data vector (see
   *            above).
   * @param[in] rec_offset Offset to be used for indexing within the data vector
   *            (see above).
   * @return Anything other than 0 denotes an error.
   */
  int parse_block(const std::vector<const char *> &sites,
                  std::vector<Vmf3SiteData> &data, int recs_per_site = 1,
                  int rec_offset = 0) noexcept;

  /* stream state according to std::basic_ios<CharT,Traits>::good */
  auto good() const noexcept { return mstream.good(); }

  /* set stream to top of file */
  void rewind() noexcept { mstream.seekg(0); }

  /* Current (just skipped/parsed) epoch */
  const MjdEpoch &current_epoch() const noexcept { return mcurrent; }

  /*  Next epoch (to be skipped/parsed); first line already buffered. */
  const MjdEpoch &next_epoch() const noexcept { return mnext; }

}; /* class Vmf3SiteFileStream */

/** A VMF3 data site stream.
 *
 * This class is meant to be used for sequential computation of site-specific
 * VMF3 refraction. Internally it uses a list of specified sites and a
 * corresponding VMF3 site data file, and it parses through the file finding
 * VMF3 data for the sites of interest in requested, sequential eopchs.
 * Example usage is when performing POD, where VMF3 computations are needed
 * for a list of sites, going forward in time.
 *
 * VMF3 coefficients (and meto data included in the data file) are computed
 * using linear interpolation between adjacent epochs.
 */
class Vmf3SiteStream {
  /* for each site we store data for two adjacent epochs */
  static constexpr const int RECS_PER_SITE = 2;
  /* each site is specified by a 4char id */
  static constexpr const int MAX_SITE_CHARS = 4;
  /* site length (of string) is MAX_SITE_CHARS + '\0' */
  static constexpr const int SITE_LEN = MAX_SITE_CHARS + 1;

  /* list of sites of interest */
  std::vector<const char *> msites;
  /* data for each site; size is len(msites) * RECS_PER_SITE. data for each
   * site are stored sequentially (i.e. mdata[0] and mdata[1] hold data
   * for epochs t_i and t_{i+1} of the first site in the msites vector
   */
  std::vector<Vmf3SiteData> mdata;
  /* VMF3 full coefficient set for each site (1-1 correspondance with
   * msites). These only depend on site coordinates hence can be computed only
   * once (at initialization).
   */
  std::vector<vmf3_details::Vmf3FullCoeffs> mvfc;
  /* starting epoch of current interval */
  MjdEpoch t0{MjdEpoch::max()};
  /* ending epoch of current interval */
  MjdEpoch t1{MjdEpoch::min()};
  /* the data file stream */
  Vmf3SiteFileStream mstream;
  /* internal memmory to hold site strings (i.e. pointers in msites are
   * directed in here)
   */
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

  /** Find and collect epochs such that t0 <= t < t1
   *
   * Go forward in the data stream untill we meet blocks such that 
   * t0 <= t < t1. Collect the blocks (for the sites of interest) for t0 and 
   * t1 and store them in the  data vector.
   *
   * Note that we can only go forward in the data file.
   */
  int forward_search(const MjdEpoch &t) noexcept;

  std::vector<vmf3_details::Vmf3FullCoeffs>::iterator
  set_site_coordinates(const char *site,
                       const GeodeticCrdConstView &crd) noexcept;

public:
  /* Constructor using the VMF3 (site) data file and a list of sites */
  Vmf3SiteStream(const char *fn, const std::vector<const char *> &sites);

  /* Destructor */
  ~Vmf3SiteStream() noexcept {
    if (mmemsites)
      delete[] mmemsites;
  }

  /** Find and collect epochs such that t0 <= t < t1
   *
   * Go forward in the data stream untill we meet blocks such that 
   * t0 <= t < t1. Collect the blocks (for the sites of interest) for t0 and 
   * t1 and store them in the  data vector.
   *
   * Note that we can only go forward in the data file.
   */
  int set_at_epoch(const MjdEpoch &t) noexcept;

  /* Return the site vector */
  const std::vector<const char *> &sites() const noexcept { return msites; }

  /* Append a new site */
  std::vector<const char *>::iterator append_site(const char *site) noexcept;

  /* Assuming that the data vector holds two epochs (entries) for each site, 
   * one for t0 and one for t1, swap the two elements
   */
  void swap_epochs() noexcept {
    static_assert(RECS_PER_SITE == 2);
    for (auto it = mdata.begin(); it != mdata.end(); it += 2)
      std::swap(*it, *(it + 1));
    std::swap(t0, t1);
    return;
  }

  /* get the current starting epoch of the interval, i.e. t0 */
  const MjdEpoch interval_start() const noexcept { return t0; }
  
  /* get the current ending epoch of the interval, i.e. t1 */
  const MjdEpoch interval_stop() const noexcept { return t1; }

  /* Interpolate the VMF3 data for the requested site at the requested epoch 
   * It is auumed that t0 <= t < t1. Hence, you should first have called the 
   * set_at_epoch function for this epoch.
   */
  int site_vmf3(const char *site, const dso::MjdEpoch &t,
                dso::Vmf3SiteData &vmf3) noexcept;

}; /* class Vmf3SiteStream */

} /* namespace dso */

#endif
