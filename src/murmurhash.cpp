#include "../include/hll/murmurhash.hpp"
#if defined(_MSC_VER)
#define HLL_FORCE_INLINE __forceinline

#include <stdlib.h>
#define HLL_ROTL64(x, y) _rotl64(x, y)
#else
#define HLL_FORCE_INLINE inline __attribute__((always_inline))
namespace hll {
  namespace detail {
    inline std::uint64_t rotl64(std::uint64_t x, std::int8_t r) {
      return (x << r) | (x >> (64 - r));
    }
  }  // namespace detail
}  // namespace hll
#define HLL_ROTL64(x, y) ::hll::detail::rotl64(x, y)
#endif

namespace hll {
  namespace detail {
    HLL_FORCE_INLINE std::uint64_t getblock64(const std::uint64_t* p, int i) {
      return p[i];
    }


    HLL_FORCE_INLINE std::uint64_t fmix64(std::uint64_t k) {
      k ^= k >> 33;
      k *= 0xff51afd7ed558ccdULL;
      k ^= k >> 33;
      k *= 0xc4ceb9fe1a85ec53ULL;
      k ^= k >> 33;

      return k;
    }
  }  // namespace detail

  // Based on the reference MurmurHash3 implementation by Austin Appleby.
  std::uint64_t murmurhash3_x64_128(
      const void* key, const int len, const std::uint64_t seed) {
    const std::uint8_t* data = (const std::uint8_t*)key;
    const int nblocks = len/16;

    std::uint64_t h1 = seed;
    std::uint64_t h2 = seed;

    const std::uint64_t c1 = 0x87c37b91114253d5ULL;
    const std::uint64_t c2 = 0x4cf5ad432745937fULL;


    const std::uint64_t* blocks = (const std::uint64_t*)(data);

    for (int i = 0; i < nblocks; i++) {
      std::uint64_t k1 = detail::getblock64(blocks, i*2+0);
      std::uint64_t k2 = detail::getblock64(blocks, i*2+1);

      k1 *= c1;
      k1 = HLL_ROTL64(k1, 31);
      k1 *= c2;
      h1 ^= k1;

      h1 = HLL_ROTL64(h1, 27);
      h1 += h2;
      h1 = h1*5+0x52dce729;

      k2 *= c2;
      k2 = HLL_ROTL64(k2, 33);
      k2 *= c1;
      h2 ^= k2;

      h2 = HLL_ROTL64(h2, 31);
      h2 += h1;
      h2 = h2*5+0x38495ab5;
    }

    const std::uint8_t* tail = (const std::uint8_t*)(data + nblocks*16);

    std::uint64_t k1 = 0;
    std::uint64_t k2 = 0;

    switch (len & 15) {
    case 15: k2 ^= ((std::uint64_t)tail[14]) << 48;
    case 14: k2 ^= ((std::uint64_t)tail[13]) << 40;
    case 13: k2 ^= ((std::uint64_t)tail[12]) << 32;
    case 12: k2 ^= ((std::uint64_t)tail[11]) << 24;
    case 11: k2 ^= ((std::uint64_t)tail[10]) << 16;
    case 10: k2 ^= ((std::uint64_t)tail[ 9]) << 8;
    case  9: k2 ^= ((std::uint64_t)tail[ 8]) << 0;
            k2 *= c2; k2  = HLL_ROTL64(k2, 33); k2 *= c1; h2 ^= k2;

    case  8: k1 ^= ((std::uint64_t)tail[ 7]) << 56;
    case  7: k1 ^= ((std::uint64_t)tail[ 6]) << 48;
    case  6: k1 ^= ((std::uint64_t)tail[ 5]) << 40;
    case  5: k1 ^= ((std::uint64_t)tail[ 4]) << 32;
    case  4: k1 ^= ((std::uint64_t)tail[ 3]) << 24;
    case  3: k1 ^= ((std::uint64_t)tail[ 2]) << 16;
    case  2: k1 ^= ((std::uint64_t)tail[ 1]) << 8;
    case  1: k1 ^= ((std::uint64_t)tail[ 0]) << 0;
            k1 *= c1; k1  = HLL_ROTL64(k1, 31); k1 *= c2; h1 ^= k1;
    }


    h1 ^= (std::uint64_t)len;
    h2 ^= (std::uint64_t)len;

    h1 += h2;
    h2 += h1;

    h1 = detail::fmix64(h1);
    h2 = detail::fmix64(h2);

    return h1 + h2;
  }
}  // namespace hll
