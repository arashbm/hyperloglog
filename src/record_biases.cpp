#include <cmath>
#include <iostream>
#include <fstream>
#include <numeric>
#include <limits>
#include <iomanip>

#include "../include/hyperloglog.hpp"

template <unsigned short precision>
void record_bias(size_t trials, unsigned int points,
    std::string out_dir) {
  std::ofstream out;
  out.open(out_dir + std::to_string(precision),
      std::ios::trunc | std::ios::out);

  std::vector<hll::hyperloglog<precision, (unsigned short int)(precision+1)>>
    hlls(trials, true);

  unsigned long max = (1u << precision)*6, inc = std::max(1ul, max/points);


  out << "{";
  out << std::setprecision(std::numeric_limits<double>::digits10 + 1);
  for(unsigned long i = 1; i <= max; i++) {
    for (size_t j = 0; j < trials; j++)
      hlls[j].insert(std::to_string(i) + "-" + std::to_string(j));

    if (i%inc == 0) {
      std::cerr << i/inc << std::endl;
      long double err = 0;
      for (const auto& h: hlls)
        err += h.measure_error(i);
      long double avg_error = err/trials;

      out << "{" << i + avg_error
        << ", " << avg_error << "}, ";
    }
  }
  out << "}";

  out.close();
}


int main() {

  size_t sample_size = 100'000;
  unsigned int points = 500;
  std::string out_dir = "./src/biases/";

  record_bias<4>(sample_size, points, out_dir);
  record_bias<5>(sample_size, points, out_dir);
  record_bias<6>(sample_size, points, out_dir);
  record_bias<7>(sample_size, points, out_dir);
  record_bias<8>(sample_size, points, out_dir);
  record_bias<9>(sample_size, points, out_dir);
  record_bias<10>(sample_size, points, out_dir);
  record_bias<11>(sample_size, points, out_dir);
  record_bias<12>(sample_size, points, out_dir);
  record_bias<13>(sample_size, points, out_dir);
  record_bias<14>(sample_size, points, out_dir);
  record_bias<15>(sample_size, points, out_dir);
  record_bias<16>(sample_size, points, out_dir);
  record_bias<17>(sample_size, points, out_dir);
  record_bias<18>(sample_size, points, out_dir);

  return 0;
}
