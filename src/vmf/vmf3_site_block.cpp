#include "geodesy/units.hpp"
#include "vmf3.hpp"

bool dso::vmf3::SiteBlock::grid_matches_between_epochs() const noexcept {
  int status = 0;
  for (int k = 0; k < 4; k++) {
    status += (mdata_t0[k].lat_deg() == mdata_t1[k].lat_deg());
    status += (mdata_t0[k].lon_deg() == mdata_t1[k].lon_deg());
  }
  return (status == 8);
}

int dso::vmf3::SiteBlock::compute_spatial_vmf3_coeffs() noexcept {
  /* The point of computation is the four nodes of a cell (around the point
   * mcrd). The lat/lon coordinates of these points should be stored in
   * mdata_t0 (or mdata_t1).
   * These points are on the ellipsoid. Computed coeffs are stored in
   * msitebc
   * @return Anything other than zero denotes an error.
   *
   * @warning Assumes mdata_t0 is already populated!
   */
  int error = 0;
  for (int k = 0; k < 4 && (!error); k++) {
    const double lat = dso::deg2rad(mdata_t0[k].lat_deg());
    const double lon = dso::deg2rad(mdata_t0[k].lon_deg());
    /* unit vector (doesn't have to be a unit vectoe, but nevertheless ...
     * )*/
    Eigen::Vector3d rsta;
    rsta.x() = std::sin(dso::DPI / 2. - lat) * std::cos(lon);
    rsta.y() = std::sin(dso::DPI / 2. - lat) * std::sin(lon);
    rsta.z() = std::cos(dso::DPI / 2. - lat);
    error += Vmf3::vmf3_spatial_coeffs(CartesianCrdConstView(rsta), msitebc[k]);
  }
  return error;
}

/** Copy data that refer to epoch t0 to t1 (i.e. mdata_t1 <- mdata_t0)
 * Now mdata_t1 will hold an exact copy of the values in mdata_t0.
 */
void dso::vmf3::SiteBlock::cp_t0tot1() noexcept {
  for (int i = 0; i < 4; i++) {
    std::memcpy(mdata_t1[i].raw_ptr(), mdata_t0[i].raw_ptr(),
                sizeof(double) * GridVmf3Data::Data::NUM_ENTRIES);
  }
}

/** Copy data that refer to epoch t1 to t0 (i.e. mdata_t0 <- mdata_t1)
 * Now mdata_t0 will hold an exact copy of the values in mdata_t1.
 */
void dso::vmf3::SiteBlock::cp_t1tot0() noexcept {
  for (int i = 0; i < 4; i++) {
    std::memcpy(mdata_t0[i].raw_ptr(), mdata_t1[i].raw_ptr(),
                sizeof(double) * GridVmf3Data::Data::NUM_ENTRIES);
  }
}