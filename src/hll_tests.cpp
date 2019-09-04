#include "catch.hpp"

#include "../include/hyperloglog.hpp"

constexpr unsigned short p = 18, sp = 25;

TEST_CASE( "counts small sets", "[sparse]" ) {
  auto h  = hll::HyperLogLog<p, sp>();
  SECTION("small cardinalities") {
    for (int i = 1; i <= 20; i++) {
      h.insert(i);
      double est = h.estimate();
      REQUIRE( est < i+1 );
      REQUIRE( i-1 < est );
    }
  }

  SECTION("only counting distincts") {
    for (int i = 0; i < 20; i++)
      for (int j = 1; j <= 20; j++)
        h.insert(j);
    double est = h.estimate();
    REQUIRE( est < 21 );
    REQUIRE( 19 < est );
  }

  SECTION("merging with itself shouldn't deadlock") {
    for (int i = 1; i<= 20; i++) {
      h.insert(i);
    }
    h.merge(h);
    double est = h.estimate();
    REQUIRE( est < 21 );
    REQUIRE( 19 < est );
  }

  SECTION("merging sparse counters") {
    auto h2 = hll::HyperLogLog<p, sp>();
    for (int i = 1; i<= 20; i++) {
      h.insert(i);
      h2.insert(i+5);
    }
    h.merge(h2);
    double est = h.estimate();
    REQUIRE( est < 26 );
    REQUIRE( 24 < est );
  }

  SECTION("merging dense counter into a sparse counter") {
    auto h2 = hll::HyperLogLog<p, sp>(true);
    for (int i = 1; i<= 20; i++) {
      h.insert(i);
      h2.insert(i+5);
    }
    h.merge(h2);
    double est = h.estimate();
    REQUIRE( est < 26 );
    REQUIRE( 24 < est );
  }
}

TEST_CASE( "counts large sets", "[dense]" ) {
  auto h = hll::HyperLogLog<p, sp>(true);
  unsigned long m = (1ul << p);
  unsigned long count = 10*m;
  long double relative_error = 3.0/std::sqrt(m); // almost 3 sigmas or %99.73
  SECTION("large cardinalities") {
    for (unsigned long i = 1; i <= count ; i++) {
      h.insert(i);
    }
    double est = h.estimate();
    REQUIRE( est < count*(1.0 + relative_error) );
    REQUIRE( count*(1.0 - relative_error) < est );
  }
  SECTION("only counting distincts") {
    for (int i = 0; i < 20; i++)
      for (unsigned long j = 1; j <= count; j++)
        h.insert(j);
    double est = h.estimate();
    REQUIRE( est < count*(1.0 + relative_error) );
    REQUIRE( count*(1.0 - relative_error) < est );
  }

  SECTION("merging dense counters") {
    auto h2 = hll::HyperLogLog<p, sp>(true);
    for (int i = 1; i<= 20; i++) {
      h.insert(i);
      h2.insert(i+5);
    }
    h.merge(h2);
    double est = h.estimate();
    REQUIRE( est < 26 );
    REQUIRE( 24 < est );
  }

  SECTION("merging sparse counter into a dense counter") {
    auto h2 = hll::HyperLogLog<p, sp>();
    for (int i = 1; i<= 20; i++) {
      h.insert(i);
      h2.insert(i+5);
    }
    h.merge(h2);
    double est = h.estimate();
    REQUIRE( est < 26 );
    REQUIRE( 24 < est );
  }
}

TEST_CASE( "counts after transitioning from sparse to dense", "[transition]" ) {
  auto h = hll::HyperLogLog<p, sp>();
  unsigned long m = (1ul << p);
  unsigned long count = 10*m;
  long double relative_error = 3.0/std::sqrt(m); // almost 3 sigmas or %99.73
  SECTION("large cardinalities") {
    for (unsigned long i = 1; i <= count ; i++) {
      h.insert(i);
      if (i % (count/100) == 0) {
        double est = h.estimate();
        REQUIRE( est < i*(1.0 + relative_error) );
        REQUIRE( i*(1.0 - relative_error) < est );
      }
    }
  }
  SECTION("only counting distincts") {
    for (int i = 0; i < 20; i++)
      for (unsigned long j = 1; j <= count; j++)
        h.insert(j);
    double est = h.estimate();
    REQUIRE( est < count*(1.0 + relative_error) );
    REQUIRE( count*(1.0 - relative_error) < est );
  }
}
