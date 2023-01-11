#ifndef INCLUDE_HLL_MURMURHASH_HPP_
#define INCLUDE_HLL_MURMURHASH_HPP_

// Based on the reference MurmurHash3 implementation by Austin Appleby. Unlike
// the reference implementation, this version only returns the first 64bits of
// the result, as well as receiving a 64-bit seed instead of a 32-bit seed.

#include <cstdint>

namespace hll {
  std::uint64_t murmurhash3_x64_128(
      const void* key, int len, std::uint64_t seed);
}  // namespace hll

#endif  // INCLUDE_HLL_MURMURHASH_HPP_
