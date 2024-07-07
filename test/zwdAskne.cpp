#include "tropo.hpp"
#include <cassert>

using namespace dso;

int main() {
  /* example copyied from https://vmf.geo.tuwien.ac.at/codes/asknewet.f */
  assert(std::abs(zwd(10.9621e0, 273.8720e0, 2.8071e0) - 0.1176e0) < 1e-4);

  return 0;
}
