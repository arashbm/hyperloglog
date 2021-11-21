#include <cmath>
#include <iostream>
#include <fstream>
#include <numeric>

#include "../include/hll/hyperloglog.hpp"
template <uint8_t precision>
void record_estimates(
    std::size_t trials,
    std::size_t points,
    std::size_t max) {
  std::vector<hll::hyperloglog<precision, 18>> hlls;

  for (size_t j = 0; j < trials; j++)
    hlls.emplace_back(true);

  std::size_t inc = std::max(1ul, max/points);

  for (std::size_t i = 1; i <= max; i++) {
    for (std::size_t j = 0; j < trials; j++)
      hlls[j].insert(std::to_string(i) + ";" + std::to_string(j));

    if (i%inc == 0) {
      double bias = 0;
      double se = 0;
      for (const auto& h: hlls) {
        double diff = h.estimate() - static_cast<double>(i);
        bias += diff;
        se += std::pow(diff, 2);
      }
      std::cout << i << " " << trials
        << " "
        << bias/static_cast<double>(trials)/static_cast<double>(i) << " "
        << std::sqrt(se/static_cast<double>(trials))/
            std::sqrt(trials)/static_cast<double>(i)<< "\n";
    }
  }
}

int main() {
  std::size_t sample_size = 5000;
  std::size_t max = 320'000'000;
  std::size_t points = 50'000;

  record_estimates<10>(sample_size, 10000, 50000);
  record_estimates<10>(sample_size, points, max);
}
