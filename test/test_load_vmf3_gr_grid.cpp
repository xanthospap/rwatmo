#include <random>

#include "vmf3_grid_stream.hpp"

/* IMPORTANT
 * This only works in debug mode
 */

using namespace dso;
constexpr const int num_tests = 10000;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Error. Usage: %s [VMF3GR GRID FILE]\n", argv[0]);
    return 1;
  }

  vmf3::GridVmf3Data grid;
  if (load_vmfgr3_grid_map(argv[1], &grid)) {
    fprintf(stderr, "Error. Failed loading grid file %s\n", argv[1]);
    return 1;
  }

  std::mt19937 rng(std::random_device{}());  // seed once
  std::uniform_real_distribution<double> rlon(
      std::min(grid.tick_axis_inner().tick_start(),
               grid.tick_axis_inner().tick_stop()),
      std::max(grid.tick_axis_inner().tick_start(),
               grid.tick_axis_inner().tick_stop()));
  std::uniform_real_distribution<double> rlat(
      std::min(grid.tick_axis_outter().tick_start(),
               grid.tick_axis_outter().tick_stop()),
      std::max(grid.tick_axis_outter().tick_start(),
               grid.tick_axis_outter().tick_stop()));

#ifdef DEBUG
  for (int i = 0; i < num_tests; i++) {
    double lat = rlat(rng);
    double lon = rlon(rng);
    try {
      const auto c = grid.cell(lat, lon);
      // printf(
      //     "(%+.2f, %+.2f)                                           (%+.2f, "
      //     "%+.2f)\n",
      //     grid.data(c.tl)->lat(), grid.data(c.tl)->lon(),
      //     grid.data(c.tr)->lat(), grid.data(c.tr)->lon());
      // printf("                                (%.2f, %.2f)\n", lat, lon);
      // printf(
      //     "(%+.2f, %+.2f)                                           (%+.2f, "
      //     "%+.2f)\n",
      //     grid.data(c.bl)->lat(), grid.data(c.bl)->lon(),
      //     grid.data(c.br)->lat(), grid.data(c.br)->lon());

      assert(grid.data(c.bl)->lat() <= lat);
      assert(grid.data(c.bl)->lon() <= lon);
      assert(grid.data(c.br)->lat() <= lat);
      assert(grid.data(c.br)->lon() > lon);
      assert(grid.data(c.tl)->lat() > lat);
      assert(grid.data(c.tl)->lon() <= lon);
      assert(grid.data(c.tr)->lat() > lat);
      assert(grid.data(c.tr)->lon() > lon);

    } catch (std::exception &e) {
      fprintf(stderr, "Failed locating cell for point (%.3f, %.3f)\n", lat,
              lon);
      fprintf(stderr, "Grid range: %+.3f < lat < %+.3f\n",
              grid.tick_axis_outter().tick_start(),
              grid.tick_axis_outter().tick_stop());
      fprintf(stderr, "          : %+.3f < lon < %+.3f\n",
              grid.tick_axis_inner().tick_start(),
              grid.tick_axis_inner().tick_stop());
      fprintf(stderr, "%s\n", e.what());
    }
  }

  for (int i = 0; i < num_tests; i++) {
    double lat = rlat(rng);
    double lon = rlon(rng);
    const auto c = grid.cell_nocheck(lat, lon);
    assert(grid.data(c.bl)->lat() <= lat);
    assert(grid.data(c.bl)->lon() <= lon);
    assert(grid.data(c.br)->lat() <= lat);
    assert(grid.data(c.br)->lon() > lon);
    assert(grid.data(c.tl)->lat() > lat);
    assert(grid.data(c.tl)->lon() <= lon);
    assert(grid.data(c.tr)->lat() > lat);
    assert(grid.data(c.tr)->lon() > lon);

    const auto c2 = grid.cell(lat, lon);
    assert(c2.bl == c.bl);
    assert(c2.br == c.br);
    assert(c2.tl == c.tl);
    assert(c2.tr == c.tr);
  }
#endif

  return 0;
}