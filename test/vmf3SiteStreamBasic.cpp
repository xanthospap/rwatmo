#include "vmf3.hpp"
#include <cstdio>
#include <vector>

using namespace dso;

const char *doris_sites[] = {"STKB", "PDOC", "SOCA", "GALA",
                             "FUUB", "WEUC", "ADEA"};
const int NUM_SITES = 7;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s <VMF3 DORIS FILE>\n", argv[0]);
    return 1;
  }

  std::vector<const char *> sites;
  for (int i = 0; i < NUM_SITES; i++)
    sites.push_back(doris_sites[i]);

  Vmf3SiteStream vstream(argv[1], sites);

  for (const auto ptr: vstream.sites()) {
    printf("Site: [%s]\n", ptr);
  }

  if (vstream.set_at_epoch(MjdEpoch(59792, FractionalSeconds(3600e0)))) {
    fprintf(stderr, "Failed getting to epoch\n");
    return 2;
  }

  Vmf3SiteData d;
  for (int sec=0e0; sec <= 86400; sec+=1) {
    const auto t = MjdEpoch(59792, FractionalSeconds((double)sec));
    /* first buffer the interval requested */
    if (vstream.set_at_epoch(t)) return 5;
    if (vstream.site_vmf3("SOCA", t, d)) return 3;
    printf("%.9f %.9f %.9f %.4f %.4f %.2f %.2f %.2f\n",
           t.imjd() + t.fractional_days().days(), d.ah(), d.aw(), d.zhd(), d.zwd(),
           d.pressure(), d.temperature(), d.water_vapour_pressure());
  }

  vstream.append_site("RAQB");
  vstream.append_site("DIOB");
  for (const auto ptr: vstream.sites()) {
    printf("Site: [%s]\n", ptr);
  }

  return 0;
}
