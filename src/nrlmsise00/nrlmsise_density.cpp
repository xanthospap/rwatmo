#include "geodesy/transformations.hpp"
#include "nrlmsise.hpp"

double dso::Nrlmsise00::density(Eigen::Vector3d rsat, const dso::MjdEpoch &tt,
                                const dso::Msise00Data &data)
{

  /* cartesian to geodetic satellite position */
  const auto llh = dso::cartesian2geodetic<dso::ellipsoid::grs80>(
      dso::CartesianCrdConstView(rsat));

  /* temperatures and densities */
  double densities[9];
  double temperatures[2];

  /* */
  return gts7(tt, llh, data.f107A, data.f107, data.ap, densities, temperatures);
}
