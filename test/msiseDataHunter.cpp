#include "nrlmsise.hpp"
#include "space_weather.hpp"
#include <cstdio>
#include "datetime/calendar.hpp"
#include "datetime/datetime_write.hpp"
#include <vector>

std::vector<dso::TwoPartDateUTC> indates;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s [CELESTRAK CSV FILE]\n", argv[0]);
    return 1;
  }

  dso::TwoPartDateUTC d1(59580, dso::FractionalSeconds(0)); // 2022/01/01
  dso::TwoPartDateUTC d2(60795, dso::FractionalSeconds(0)); // 2025/04/30

  /* collect space weather data fro CSV, in range [d1,d2) */
  const auto swdv = dso::load_celestrak_sw(argv[1], d1.utc2tt(), d2.utc2tt());
  //for (const auto &it : swdv) {
  //  printf("Date: %.1f [TT] ", it.tt().as_mjd());
  //  for (int i = 0; i < 22; i++) {
  //    printf("%d ", it.int_array()[i]);
  //  }
  //  printf("\n");
  //}
  
  dso::NrlmsiseDataHunter hunter;
  hunter.set_data_ptr(swdv);

  /* 2024/01/03 is MJD 60312 */
  dso::TwoPartDateUTC utc(60312, dso::FractionalSeconds(60 * 60 * 3 - 10.));
  hunter.get_data(utc.utc2tt());

  const dso::Msise00Data *ptr = &(hunter.msise_data());
  //printf("f107A: %.2f\n", ptr->f107A);
  //printf("f107 : %.2f\n", ptr->f107);
  //for (int i = 0; i < 7; i++)
  //  printf("ap[%d]: %.2f\n", i, ptr->ap[i]);

  //utc = dso::TwoPartDateUTC(60312, dso::FractionalSeconds(60 * 60 * 3 + 10.));
  //hunter.get_data(utc.utc2tt());

  //printf("f107A: %.2f\n", ptr->f107A);
  //printf("f107 : %.2f\n", ptr->f107);
  //for (int i = 0; i < 7; i++)
  //  printf("ap[%d]: %.2f\n", i, ptr->ap[i]);

  /* input, test dates */
  indates.emplace_back(dso::TwoPartDateUTC(60693, dso::FractionalSeconds(60*60*10+1)));
  indates.emplace_back(dso::TwoPartDateUTC(60786, dso::FractionalSeconds(60*60*23+1)));
  indates.emplace_back(dso::TwoPartDateUTC(60618, dso::FractionalSeconds(60*60* 6+1)));
  indates.emplace_back(dso::TwoPartDateUTC(60532, dso::FractionalSeconds(60*60*15+1)));
  indates.emplace_back(dso::TwoPartDateUTC(60439, dso::FractionalSeconds(60*60* 1+1)));
  indates.emplace_back(dso::TwoPartDateUTC(60309, dso::FractionalSeconds(60*60* 0+1)));
  indates.emplace_back(dso::TwoPartDateUTC(60197, dso::FractionalSeconds(60*60*12+1)));
  indates.emplace_back(dso::TwoPartDateUTC(60075, dso::FractionalSeconds(60*60*20+1)));
  indates.emplace_back(dso::TwoPartDateUTC(59963, dso::FractionalSeconds(60*60*19+1)));
  indates.emplace_back(dso::TwoPartDateUTC(59614, dso::FractionalSeconds(60*60* 2+1)));

  /* process input dates */
  for (const auto &d : indates) {
    try {
    hunter.get_data(d.utc2tt());
    printf("%.2f,", ptr->f107A);
    printf("%.2f,", ptr->f107);
    for (int i = 0; i < 7; i++)
      printf("%.2f,", ptr->ap[i]);
    } catch (std::exception &e) {
      fprintf(stderr, "Failed getting space weather data for date; what: %s\n", e.what());
      char buf[64];
      fprintf(stderr, "Error triggered for date %s\n", dso::to_char<dso::YMDFormat::YYYYMMDD, dso::HMSFormat::HHMMSS>(d, buf));
    }
    printf("\n");
  }

  return 0;
}
