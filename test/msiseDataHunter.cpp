#include "nrlmsise.hpp"
#include <cstdio>
#include <datetime/fractional_types.hpp>
#include <datetime/tpdate.hpp>
#include "space_weather.hpp"

int main(int argc, char *argv[]) {
  if (argc !=2) {
    fprintf(stderr, "Usage: %s [CELESTRAK CSV FILE]\n", argv[0]);
    return 1;
  }

  /* collect space weather data fro CSV */
  const auto swdv = dso::load_celestrak_sw(argv[1], dso::MjdEpoch(60310, dso::FractionalSeconds(0)), dso::MjdEpoch(60315, dso::FractionalSeconds(0)));
  if (swdv.size() != 4) {
    fprintf(stderr, "Data collection not successeful! expected records 4, got %d\n", (int)swdv.size());
    return 1;
  }

  dso::NrlmsiseDataHunter hunter;
  hunter.set_data_ptr(swdv);

  dso::TwoPartDateUTC utc(60312, dso::FractionalSeconds(60*60*3+1));
  hunter.get_data(utc.utc2tt());
  
  const dso::Msise00Data *ptr = &(hunter.msise_data());
  printf("f107 : %.2f\n", ptr->f107);
  printf("f107A: %.2f\n", ptr->f107A);
  for (int i=0; i<7; i++) printf("ap[%d]: %.2f\n", i, ptr->ap[i]);

  return 0;
}
