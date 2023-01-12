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


std::pair<std::uint64_t, std::uint8_t> get_hash_rank(std::uint64_t hash) {
  std::uint8_t precision;
  precision = p;

  std::uint64_t index = (std::uint64_t)(hash >> (sizeof(hash)*8 - precision));

  std::uint8_t rank = static_cast<std::uint8_t>(sizeof(hash)*8 - precision);
  std::uint64_t h = hash << precision;
  if (h > 0)
    rank = std::min(rank, static_cast<std::uint8_t>(hll_countl_zero(h) + 1));
  return std::make_pair(index, rank);
}

TEST_CASE("The weird case of 350285", "[WTF]") {
  hll::hyperloglog<std::size_t, p, sp> h;
  for (std::size_t i = 1; i < 350284; i++)
    h.insert(i);

  std::cerr << "estimate: " << h.estimate() << std::endl;
  std::uint64_t seed = 0x9E3779B97F4A7C15ul;

  std::size_t item = 350284;
  std::cerr << "item: " << item << std::endl;
  std::uint64_t hash = hll::hash<std::size_t>{}(item, seed);
  std::cerr << "hash: " << hash << std::endl;

  std::uint64_t index;
  std::uint8_t rank;
  std::tie(index, rank) = get_hash_rank(hash);
  std::cerr << "index: " << index <<
    " rank: " << (std::uint64_t)rank << std::endl;
  std::cerr << "current_rank: " <<
    (std::uint64_t)h.dense_vec()[index] << std::endl;
  h.insert(item);
  std::cerr << "current_rank: " <<
    (std::uint64_t)h.dense_vec()[index] << std::endl;
  std::cerr << "estimate: " << h.estimate() << std::endl;

  item = 350285;
  std::cerr << "item: " << item << std::endl;
  hash = hll::hash<std::size_t>{}(item, seed);
  std::cerr << "hash: " << hash << std::endl;
  std::tie(index, rank) = get_hash_rank(hash);
  std::cerr << "index: " << index <<
    " rank: " << (std::uint64_t)rank << std::endl;
  std::cerr << "current_rank: " <<
    (std::uint64_t)h.dense_vec()[index] << std::endl;
  h.insert(item);
  std::cerr << "current_rank: " <<
    (std::uint64_t)h.dense_vec()[index] << std::endl;
  std::cerr << "estimate: " << h.estimate() << std::endl;

  for (std::size_t i = 0; i < h.dense_vec().size(); i++) {
    if (i % 128 == 0)
      std::cerr << std::endl;
    std::cerr << std::hex << (std::uint64_t)h.dense_vec()[i] << std::dec;
  }
  std::cerr << std::endl;
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
