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

  Vmf3SiteStream(argv[1], sites);

  return 0;
}
