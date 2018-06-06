#include "hyperloglog.hpp"


#include <cmath>
#include <iostream>
#include <numeric>

int main() {

  constexpr int p = 17, sample_size = 1000;

  std::vector<hll::HyperLogLog<p, p+1>> hlls(sample_size);

  unsigned long max = (1u << p)*5+1, inc = std::max(1ul, max/400);


  for(unsigned long i = 1; i <= max; i++) {
    for (int j = 0; j < sample_size; j++)
      hlls[j].insert(std::to_string(i) + "-" + std::to_string(j));

    if (i%inc == 0) {
      long double est = 0;
      for (const auto& h: hlls)
        est += h.estimate();//measure_error(i);
      long double avg_estimate = est/sample_size;
/*
      double est_dists = 0;
      for (const auto& h: hlls)
        est_dists += std::pow(h.estimate() - avg_estimate, 2);
      double sigma = std::sqrt(est_dists/sample_size);
      */

      std::cout << i
        << ", " << avg_estimate << std::endl;
        //<< ", " << sigma << std::endl;
    }
  }

  return 0;
}
