#include "datetime/datetime_write.hpp"
#include "vmf3.hpp"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

using namespace dso;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Error. Usage: %s [VMF3GR DATA_DIR]\n", argv[0]);
    return 1;
  }

  /* setup some test sites */
  const char *sites[] = {"PDOC", "SVBC", "BEMB", "FUUC",
                         "COBB", "OWGC", "SCSC", "GAVC"};
  std::vector<Eigen::Vector3d> crds;
  Eigen::Vector3d tmp;
  tmp << 4551593.731, -2186898.736, 3883409.506;
  crds.push_back(tmp);
  tmp << 1201300.042, 251874.432, 6238000.308;
  crds.push_back(tmp);
  tmp << 1106046.627, -763739.010, -6214243.195;
  crds.push_back(tmp);
  tmp << -6178323.859, -202689.926, -1566023.236;
  crds.push_back(tmp);
  tmp << -3484298.102, -1084773.365, 5213542.261;
  crds.push_back(tmp);
  tmp << -4584391.539, -290926.021, -4410054.112;
  crds.push_back(tmp);
  tmp << -33879.940, -6377521.954, -82085.574;
  crds.push_back(tmp);
  tmp << 4783636.154, 2140691.252, 3623260.379;
  crds.push_back(tmp);

  const int num_sites = (int)crds.size();

  /* instance */
  Vmf3SiteHandler tropo(argv[1]);

  /* add sites */
  for (int i = 0; i < num_sites; i++) {
    assert(!tropo.append_site(sites[i], dso::CartesianCrdConstView(crds[i])));
  }

  /* some starting date */
  MjdEpoch t(60310, FractionalSeconds(0));

  vmf3::Vmf3Result res;
  char buf[64];
  /* results here */
  while (t < MjdEpoch(60312, FractionalSeconds(0))) {
    for (int i = 0; i < num_sites; i++) {
      if (tropo.vmf3(sites[i], t, DPI / 2., res)) {
        fprintf(stderr, "Error. Failed computing vmf3 for site %s at t=%.9f\n",
                sites[i], t.as_mjd());
        return 9;
      }
      printf("%s %s %.3f %.3f\n",
             to_char<YMDFormat::YYYYMMDD, HMSFormat::HHMMSSF>(t, buf), sites[i],
             res.zhd(), res.zwd());
    }
    t.add_seconds_inplace(FractionalSeconds(30));
  }

  return 0;
}
