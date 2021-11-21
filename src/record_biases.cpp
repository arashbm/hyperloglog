#include <cmath>
#include <iostream>
#include <fstream>
#include <numeric>
#include <limits>
#include <iomanip>

#include "../include/hll/hyperloglog.hpp"

template <uint8_t precision>
void record_bias(std::size_t trials, std::size_t points,
    std::string out_dir) {
  std::ofstream out;
  out.open(out_dir + std::to_string(precision),
      std::ios::trunc | std::ios::out);

  using hll_t = hll::hyperloglog<precision, static_cast<uint8_t>(precision+1)>;

  std::vector<hll_t> hlls(trials, hll_t(true));

  std::size_t max = (1u << precision)*6;
  std::size_t inc = std::max(1ul, max/points);

  out << "{";
  out << std::setprecision(std::numeric_limits<double>::digits10 + 1);
  for (std::size_t i = 1; i <= max; i++) {
    for (std::size_t j = 0; j < trials; j++)
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
  std::size_t sample_size = 100'000;
  std::size_t points = 500;
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
