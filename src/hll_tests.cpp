#include <catch2/catch.hpp>

#include <iostream>

#include <hll/murmurhash.hpp>

TEST_CASE("MurmurHash3 implementation", "[murmurhash]") {
  std::string data = "The quick brown fox jumps over the lazy dog";
  std::uint64_t seed = 0;
  REQUIRE(
      hll::murmurhash3_x64_128(
        data.c_str(),
        static_cast<int>(data.size()),
        seed) == 16378391709484522348ul);

  std::string data2 = "THE QUICK BROWN FOX JUMPS OVER THE LAZY DOG";
  REQUIRE(
      hll::murmurhash3_x64_128(
        data2.c_str(),
        static_cast<int>(data2.size()),
        seed) == 11970594202964392905ul);

  std::uint64_t number = 350285;
  std::uint64_t seed2 = 0x9E3779B97F4A7C15ul;
  REQUIRE(
      hll::murmurhash3_x64_128(
        &number, sizeof(number), seed2) == 8023538134681085539ul);
}

#include <hll/hyperloglog.hpp>

constexpr uint8_t p = 18, sp = 25;

TEST_CASE("counts small sets", "[sparse]") {
  hll::hyperloglog<std::size_t, p, sp> h;
  SECTION("small cardinalities") {
    for (std::size_t i = 1; i <= 20; i++) {
      h.insert(i);
      double est = h.estimate();
      REQUIRE(est < static_cast<double>(i+1));
      REQUIRE(static_cast<double>(i-1) < est);
    }
  }

  SECTION("only counting distincts") {
    for (std::size_t i = 0; i < 20; i++)
      for (std::size_t j = 1; j <= 20; j++)
        h.insert(j);
    double est = h.estimate();
    REQUIRE(est < 21);
    REQUIRE(19 < est);
  }

  SECTION("merging sparse counters") {
    hll::hyperloglog<std::size_t, p, sp> h2;
    for (std::size_t i = 1; i<= 20; i++) {
      h.insert(i);
      h2.insert(i+5);
    }
    h.merge(h2);
    double est = h.estimate();
    REQUIRE(est < 26);
    REQUIRE(24 < est);
  }

  SECTION("merging dense counter into a sparse counter") {
    hll::hyperloglog<std::size_t, p, sp> h2(true);
    for (std::size_t i = 1; i<= 20; i++) {
      h.insert(i);
      h2.insert(i+5);
    }
    h.merge(h2);
    double est = h.estimate();
    REQUIRE(est < 26);
    REQUIRE(24 < est);
  }
}

TEST_CASE("counts large sets", "[dense]") {
  hll::hyperloglog<std::size_t, p, sp> h(true);
  std::size_t m = (1ul << p);
  std::size_t count = 10*m;
  long double relative_error = 3.0/std::sqrt(m);  // almost 3 sigmas or %99.73
  SECTION("large cardinalities") {
    for (std::size_t i = 1; i <= count ; i++) {
      h.insert(i);
    }
    double est = h.estimate();
    REQUIRE(est < count*(1.0 + relative_error));
    REQUIRE(count*(1.0 - relative_error) < est);
  }
  SECTION("only counting distincts") {
    for (std::size_t i = 0; i < 20; i++)
      for (std::size_t j = 1; j <= count; j++)
        h.insert(j);
    double est = h.estimate();
    REQUIRE(est < count*(1.0 + relative_error));
    REQUIRE(count*(1.0 - relative_error) < est);
  }

  SECTION("merging dense counters") {
    hll::hyperloglog<std::size_t, p, sp> h2(true);
    for (std::size_t i = 1; i<= 20; i++) {
      h.insert(i);
      h2.insert(i+5);
    }
    h.merge(h2);
    double est = h.estimate();
    REQUIRE(est < 26);
    REQUIRE(24 < est);
  }

  SECTION("merging sparse counter into a dense counter") {
    hll::hyperloglog<std::size_t, p, sp> h2;
    for (std::size_t i = 1; i<= 20; i++) {
      h.insert(i);
      h2.insert(i+5);
    }
    h.merge(h2);
    double est = h.estimate();
    REQUIRE(est < 26);
    REQUIRE(24 < est);
  }
}

double estimate_bias(double est) {
  auto bias = hll::biases[p-4];
  constexpr std::ptrdiff_t k = 6;  // K-nn parameter
  std::vector<std::pair<double, double>> keys(k);
  auto est_it = std::lower_bound(bias.begin(), bias.end(),
      std::make_pair(est, 0.0));
  std::ptrdiff_t est_idx = est_it - bias.begin();
  std::ptrdiff_t ssize = static_cast<std::ptrdiff_t>(bias.size());
  std::partial_sort_copy(
      est_idx <= k ? bias.begin() : est_it - k,
      est_idx + k >= ssize ? bias.end() : est_it + k,
      keys.begin(), keys.end(),
      [est] (const std::pair<double, double>& a,
              const std::pair<double, double>& b) {
              return std::abs(a.first - est) < std::abs(b.first - est);
      });

  std::cerr << "keys:" << std::endl;
  for (auto& k: keys)
    std::cerr << "{" << k.first << ", " << k.second << "}, ";
  std::cerr << std::endl;

  double sum = 0;
  double weight_sum = 0;
  std::tie(sum, weight_sum) = std::accumulate(keys.begin(), keys.begin()+k,
      std::make_pair(0.0, 0.0),
      [est] (
        std::pair<double, double> t,  // (sum, sum of weights)
        std::pair<double, double> key) {  // (distance, bias)
      return std::make_pair(
        t.first + key.second*1.0/std::abs(key.first-est),
        t.second + 1.0/std::abs(key.first-est));
      });

  std::cerr << "distances:" << std::endl;
  for (auto& k: keys)
    std::cerr << k.first - est << ", ";
  std::cerr << std::endl;

  std::cerr << "sum: " << sum << std::endl;
  std::cerr << "weight_sum: " << weight_sum << std::endl;

  return sum/weight_sum;
}

TEST_CASE("The weird case of 350285", "[WTF]") {
  hll::hyperloglog<std::size_t, p, sp> h;
  for (std::size_t i = 1; i <= 350285; i++)
    h.insert(i);

  std::cerr << "estimate: " << h.estimate() << std::endl;

  double bias = estimate_bias(409892.0);
  std::cerr << "bias: " << bias << std::endl;
}

TEST_CASE("counts after transitioning from sparse to dense", "[transition]") {
  hll::hyperloglog<std::size_t, p, sp> h;
  std::size_t m = (1ul << p);
  std::size_t count = 10*m;
  long double relative_error = 3.0/std::sqrt(m);  // almost 3 sigmas or %99.73

  SECTION("coorectly transitions") {
    hll::hyperloglog<std::size_t, p, sp> h2(true);
    for (std::size_t i = 1; h.is_sparse(); i++) {
      h.insert(i);
      h2.insert(i);
    }
    REQUIRE(h.estimate() == h2.estimate());

    h2.merge(h);
    REQUIRE(h.estimate() == h2.estimate());
  }

  SECTION("large cardinalities") {
    for (std::size_t i = 1; i <= count/5; i++) {
      h.insert(i);
      if (i > 350280 && i < 350290) {
        double est = h.estimate();
        std::cerr << i << " " << h.is_sparse() << " " <<
          i*(1.0 - relative_error) << " " << est << " " <<
          i*(1.0 + relative_error) <<  std::endl;
        REQUIRE(est < i*(1.0 + relative_error));
        REQUIRE(i*(1.0 - relative_error) < est);
      }
    }
  }
  SECTION("only counting distincts") {
    for (std::size_t i = 0; i < 20; i++)
      for (std::size_t j = 1; j <= count; j++)
        h.insert(j);
    double est = h.estimate();
    REQUIRE(est < count*(1.0 + relative_error));
    REQUIRE(count*(1.0 - relative_error) < est);
  }
}
