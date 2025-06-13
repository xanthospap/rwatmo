#include "nrlmsise.hpp"
#include "space_weather.hpp"
#include <cstdio>
#include <datetime/fractional_types.hpp>
#include <datetime/tpdate.hpp>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s [CELESTRAK CSV FILE]\n", argv[0]);
    return 1;
  }

  dso::TwoPartDateUTC d1(60309, dso::FractionalSeconds(0)); // 2023/12/31
  dso::TwoPartDateUTC d2(60315, dso::FractionalSeconds(0)); // 2024/01/06

  /* collect space weather data fro CSV, in range [d1,d2) */
  const auto swdv = dso::load_celestrak_sw(argv[1], d1.utc2tt(), d2.utc2tt());
  for (const auto &it : swdv) {
    printf("Date: %.1f [TT] ", it.tt().as_mjd());
    for (int i = 0; i < 22; i++) {
      printf("%d ", it.int_array()[i]);
    }
    printf("\n");
  }
  if (swdv.size() != 6) {
    fprintf(stderr,
            "Data collection not successeful! expected records 6, got %d\n",
            (int)swdv.size());
    return 1;
  }

  dso::NrlmsiseDataHunter hunter;
  hunter.set_data_ptr(swdv);

  /* 2024/01/03 is MJD 60312 */
  dso::TwoPartDateUTC utc(60312, dso::FractionalSeconds(60 * 60 * 3 - 10.));
  hunter.get_data(utc.utc2tt());

  const dso::Msise00Data *ptr = &(hunter.msise_data());
  printf("f107A: %.2f\n", ptr->f107A);
  printf("f107 : %.2f\n", ptr->f107);
  for (int i = 0; i < 7; i++)
    printf("ap[%d]: %.2f\n", i, ptr->ap[i]);

  utc = dso::TwoPartDateUTC(60312, dso::FractionalSeconds(60 * 60 * 3 + 10.));
  hunter.get_data(utc.utc2tt());

  printf("f107A: %.2f\n", ptr->f107A);
  printf("f107 : %.2f\n", ptr->f107);
  for (int i = 0; i < 7; i++)
    printf("ap[%d]: %.2f\n", i, ptr->ap[i]);

  return 0;
}
