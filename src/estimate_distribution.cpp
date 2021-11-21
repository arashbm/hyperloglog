#include <cmath>
#include <iostream>
#include <fstream>
#include <numeric>

#include "../include/hyperloglog.hpp"
template <unsigned short precision>
void record_estimates(size_t trials, unsigned int points, size_t max) {

  std::vector<hll::hyperloglog<precision, 18>> hlls;

  for (size_t j = 0; j < trials; j++)
    hlls.emplace_back(true);

  size_t inc = std::max(1ul, max/points);

  for(unsigned long i = 1; i <= max; i++) {
    for (size_t j = 0; j < trials; j++)
      hlls[j].insert(std::to_string(i) + ";" + std::to_string(j));

    if (i%inc == 0) {
      double bias = 0;
      double se = 0;
      for (const auto& h: hlls) {
        double diff = h.estimate() - (double)i;
        bias += diff;
        se += std::pow(diff, 2);
      }
      std::cout << i << " " << trials
        << " " << bias/(double)trials/(double)i << " "
        << std::sqrt(se/(double)trials)/std::sqrt(trials)/(double)i << "\n";
    }
  }
}

int main() {
  size_t sample_size = 5000;
  size_t max = 320'000'000;
  unsigned int points = 50000;

  record_estimates<10>(sample_size, 10000, 50000);
  record_estimates<10>(sample_size, points, max);
}
