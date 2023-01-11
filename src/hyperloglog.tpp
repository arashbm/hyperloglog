/**
 *
 * HyperLogLog++ Implementation
 *
 * Arash Badie Modiri <arashbm@gmail.com>
 *
 */

#include <string>
#include <algorithm>
#include <functional>
#include <cmath>
#include <numeric>
#include <type_traits>

// switch to std::countl_zero when only supporting 
#ifdef _MSC_VER
  #define hll_countl_zero __lzcnt64
#else
  #define hll_countl_zero __builtin_clzl
#endif

#define HLL_RANK_BITS 6  // == log2(64)

#include <MurmurHash3.h>


namespace hll {
#define hll_define_integral_hash(Tp)                \
  template<>                                        \
  struct hash<Tp> {                                 \
    uint64_t                                        \
    operator()(                                     \
        const Tp & k, uint32_t seed) const {        \
      uint64_t hash_out[2];                         \
      MurmurHash3_x64_128(&k, sizeof(k),            \
          seed, hash_out);                          \
      return hash_out[1];                           \
    }                                               \
  };

  hll_define_integral_hash(bool)  // NOLINT
  hll_define_integral_hash(char)  // NOLINT
  hll_define_integral_hash(signed char)  // NOLINT
  hll_define_integral_hash(unsigned char)  // NOLINT
  hll_define_integral_hash(wchar_t)  // NOLINT
  hll_define_integral_hash(char16_t)  // NOLINT
  hll_define_integral_hash(char32_t)  // NOLINT
  hll_define_integral_hash(short)  // NOLINT
  hll_define_integral_hash(int)  // NOLINT
  hll_define_integral_hash(long)  // NOLINT
  hll_define_integral_hash(long long)  // NOLINT
  hll_define_integral_hash(unsigned short)  // NOLINT
  hll_define_integral_hash(unsigned int)  // NOLINT
  hll_define_integral_hash(unsigned long)  // NOLINT
  hll_define_integral_hash(unsigned long long)  // NOLINT

#undef hll_define_integral_hash

  // making sure hash(-0.0f) == hash(+0.0f)
#define hll_define_floating_point_hash(Tp)          \
  template<>                                        \
  struct hash<Tp> {                                 \
    uint64_t                                        \
    operator()(                                     \
        const Tp & k, uint32_t seed) const {        \
      Tp zero = 0.0;                                \
      uint64_t hash_out[2];                         \
      MurmurHash3_x64_128(                          \
          ((k == 0.0) ? &zero : &k),                \
          sizeof(k), seed, hash_out);               \
      return hash_out[1];                           \
    }                                               \
  };

  hll_define_floating_point_hash(float)  // NOLINT
  hll_define_floating_point_hash(double)  // NOLINT
  hll_define_floating_point_hash(long double)  // NOLINT

#undef _hll_define_floating_point_hash


  template<> struct hash<std::string> {
    uint64_t
    operator()(const std::string& k, uint32_t seed) const {
      uint64_t hash_out[2];
      MurmurHash3_x64_128(
          k.c_str(),
          static_cast<int>(k.length())+1,
          seed, hash_out);
      return hash_out[1];
    }
  };

  const std::vector<std::pair<double, double>> biases[15] = {
#include "biases/4"
    ,
#include "biases/5"
    ,
#include "biases/6"
    ,
#include "biases/7"
    ,
#include "biases/8"
    ,
#include "biases/9"
    ,
#include "biases/10"
    ,
#include "biases/11"
    ,
#include "biases/12"
    ,
#include "biases/13"
    ,
#include "biases/14"
    ,
#include "biases/15"
    ,
#include "biases/16"
    ,
#include "biases/17"
    ,
#include "biases/18"
  };
}  // namespace hll

template <typename T, uint8_t p, uint8_t sp>
hll::hyperloglog<T, p, sp>::hyperloglog(bool create_dense, uint32_t seed)
  : seed(seed) {
  if (create_dense) {
    sparse = false;
    convert_to_dense();
  } else {
    sparse = true;
    temporary_list.reserve(temporary_list_max);
  }
}

template <typename T, uint8_t p, uint8_t sp>
double
hll::hyperloglog<T, p, sp>::estimate_bias(double est) const {
  constexpr std::size_t k = 6;  // K-nn parameter
  std::vector<std::pair<double, double>> keys(k);
  auto est_it = std::lower_bound(bias.begin(), bias.end(),
      std::make_pair(est, 0.0));
  std::partial_sort_copy(
      std::max(est_it - k, bias.begin()),
      std::min(est_it + k, bias.end()),
      keys.begin(), keys.end(),
      [est] (std::pair<double, double> a,
              std::pair<double, double> b) {
              return std::abs(a.first - est) < std::abs(b.first - est);
      });

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
  return sum/weight_sum;
}


template <typename T, uint8_t p, uint8_t sp>
const std::vector<std::pair<double, double>>
hll::hyperloglog<T, p, sp>::bias = hll::biases[p-4];


template <typename T, uint8_t p, uint8_t sp>
constexpr double
hll::hyperloglog<T, p, sp>::alpha() const {
  if (p == 4)
    return 0.673;
  else if (p == 5)
    return 0.697;
  else if (p == 6)
    return 0.709;
  else
    return 0.7213/(1.0+1.079/(1ul << p));
}

template <typename T, uint8_t p, uint8_t sp>
constexpr double
hll::hyperloglog<T, p, sp>::threshold() const {
  double thresholds[] = { 10,     20,     40,     80,     220,
                          400,    900,    1800,   3100,   6500,
                          11500,  20000,  50000,  120000, 350000};
  return thresholds[p-4];
}

template <typename T, uint8_t p, uint8_t sp>
double hll::hyperloglog<T, p, sp>::linear_estimate(std::size_t non_zero) const {
  double m;
  if (sparse)
    m = static_cast<double>(1ul << sp);
  else
    m = static_cast<double>(1ul << p);

  return m*std::log(m/(m-static_cast<double>(non_zero)));
}

template <typename T, uint8_t p, uint8_t sp>
std::pair<uint64_t, uint8_t>
hll::hyperloglog<T, p, sp>::get_hash_rank(uint64_t hash) const {
  uint8_t precision;
  if (sparse)
    precision = sp;
  else
    precision = p;

  uint64_t index = (uint64_t)(hash >> (sizeof(hash)*8 - precision));
  uint8_t rank = (uint8_t)(std::min(sizeof(hash)*8 - precision,
      (std::size_t)(hll_countl_zero(hash << precision)) + 1));
  return std::make_pair(index, rank);
}

template <typename T, uint8_t p, uint8_t sp>
uint64_t
hll::hyperloglog<T, p, sp>::encode_hash(uint64_t index, uint8_t rank) const {
  return (index << HLL_RANK_BITS) | rank;
}

template <typename T, uint8_t p, uint8_t sp>
std::pair<uint64_t, uint8_t>
hll::hyperloglog<T, p, sp>::decode_hash(uint64_t hash) const {
  uint64_t index = (hash >> HLL_RANK_BITS);
  uint8_t rank = (uint8_t)(((1 << HLL_RANK_BITS)-1) & hash);  // first 6 bits
  return std::make_pair(index, rank);
}



template <typename T, uint8_t p, uint8_t sp>
void hll::hyperloglog<T, p, sp>::merge(
    const hll::hyperloglog<T, p, sp> &other) {
  if (seed != other.seed)
    throw std::invalid_argument(
        "two counters should have the same seed to merge");

  if (other.sparse && sparse) {
      merge_temp();
      std::vector<uint64_t> other_slist = other.merged_temp_list();
      std::vector<uint64_t> merged_slist = merged_sorted_list(other_slist);
      sparse_list.swap(merged_slist);
  } else {
    std::vector<uint8_t> other_converted;
    if (sparse)
      convert_to_dense();

    auto other_dense_begin = other.dense.begin();
    if (other.sparse) {
      other_converted = other.converted_to_dense();
      other_dense_begin = other_converted.begin();
    }

    std::transform(dense.begin(), dense.end(),
        other_dense_begin,
        dense.begin(),
        [] (const uint8_t a, const uint8_t b) {return std::max(a, b);});
  }
}


template <typename T, uint8_t precision, uint8_t sparse_precision>
void hll::hyperloglog<T, precision, sparse_precision>::insert(T item) {
  uint64_t hash = hll::hash<T>{}(item, seed);
  uint64_t index;
  uint8_t rank;
  std::tie(index, rank) = get_hash_rank(hash);

  if (sparse) {
    uint64_t encoded = encode_hash(index, rank);
    temporary_list.push_back(encoded);

    if (temporary_list.size() >= temporary_list_max)
      merge_temp();

    if (sparse_list.size() >= sparse_list_max)
      convert_to_dense();
  } else if (rank > dense[index]) {
    dense[index] = rank;
  }
}


template <typename T, uint8_t precision, uint8_t sparse_precision>
void
hll::hyperloglog<T, precision, sparse_precision>::merge_temp() {
  std::vector<uint64_t> new_sparse_list = merged_temp_list();
  sparse_list.swap(new_sparse_list);
  temporary_list.clear();
}


template <typename T, uint8_t precision, uint8_t sparse_precision>
std::vector<uint64_t>
hll::hyperloglog<T, precision, sparse_precision>::merged_sorted_list(
    const std::vector<uint64_t> sorted_list) const {
  std::vector<uint64_t> new_sparse_list;

  auto it1 = sorted_list.begin();
  auto it2 = sparse_list.begin();

  int i = 0;
  while (it1 != sorted_list.end() && it2 != sparse_list.end()) {
    i++;
    uint64_t index1;
    uint8_t rank1;
    std::tie(index1, rank1) = decode_hash(*it1);
    uint64_t index2;
    uint8_t rank2;
    std::tie(index2, rank2) = decode_hash(*it2);

    if (index1 == index2) {
      new_sparse_list.push_back(encode_hash(index1, std::max(rank1, rank2)));
      ++it1;
      ++it2;
    } else if (index1 > index2) {
      new_sparse_list.push_back(encode_hash(index2, rank2));
      ++it2;
    } else {
      new_sparse_list.push_back(encode_hash(index1, rank1));
      ++it1;
    }
  }


  if (it1 == sorted_list.end() && it2 != sparse_list.end())
    new_sparse_list.insert(new_sparse_list.end(),
        it2, sparse_list.end());
  else if (it2 == sparse_list.end() && it1 != sorted_list.end())
    new_sparse_list.insert(new_sparse_list.end(),
        it1, sorted_list.end());

  return new_sparse_list;
}


template <typename T, uint8_t precision, uint8_t sparse_precision>
std::vector<uint64_t>
hll::hyperloglog<T, precision, sparse_precision>::merged_temp_list() const {
  std::vector<uint64_t> temporary_list_copy = temporary_list;
  std::sort(temporary_list_copy.begin(), temporary_list_copy.end());

  std::reverse(temporary_list_copy.begin(), temporary_list_copy.end());

  auto last = std::unique(
      temporary_list_copy.begin(),
      temporary_list_copy.end(),
      [this] (uint64_t a, uint64_t b) -> bool {
        uint64_t index1, index2;
        std::tie(index1, std::ignore) = decode_hash(a);
        std::tie(index2, std::ignore) = decode_hash(b);
        return index1 == index2;
      });

  temporary_list_copy.erase(last, temporary_list_copy.end());
  std::reverse(temporary_list_copy.begin(), temporary_list_copy.end());

  return merged_sorted_list(temporary_list_copy);
}


template <typename T, uint8_t precision, uint8_t sparse_precision>
std::vector<uint8_t>
hll::hyperloglog<T, precision, sparse_precision>::converted_to_dense()
  const {
  std::vector<uint8_t> new_dense(1ul << precision, (uint8_t)0);

  std::vector<uint64_t> slist = merged_temp_list();
  for (const auto i: slist) {
    uint64_t index;
    uint8_t rank;
    std::tie(index, rank) = decode_hash(i);

    uint64_t dense_index = (index >> (sparse_precision - precision));
    uint64_t betweens = (index & ((1u << (sparse_precision - precision)) - 1));
    uint8_t dense_rank;

    if (betweens == 0)
      dense_rank = (uint8_t)(rank + (sparse_precision - precision));
    else
      dense_rank = (uint8_t)
        (((uint8_t)hll_countl_zero(betweens) -
          (sizeof(betweens)*8 - (sparse_precision - precision))) + 1);
    if (dense_rank > new_dense[dense_index])
      new_dense[dense_index] = dense_rank;
  }
  return new_dense;
}

template <typename T, uint8_t precision, uint8_t sparse_precision>
void
hll::hyperloglog<T, precision, sparse_precision>::convert_to_dense() {
  dense = converted_to_dense();

  temporary_list.clear();
  temporary_list.shrink_to_fit();

  sparse = false;
  sparse_list.clear();
  sparse_list.shrink_to_fit();
}


template <typename T, uint8_t precision, uint8_t sparse_precision>
std::pair<double, std::size_t>
hll::hyperloglog<T, precision, sparse_precision>::raw_estimate() const {
  if (sparse) {
    throw std::logic_error(
        "`raw_estimate()` does not work with sparse representation.");
  } else {
    double sum = 0;
    std::size_t non_zeros = 0;
    for (auto&& m_j: dense) {
      if (m_j > 0) non_zeros += 1;
      sum += 1.0/static_cast<double>(1ul << m_j);
    }

    return std::make_pair(
        alpha()*std::pow((1ul << precision), 2)/sum,
        non_zeros);
  }
}

template <typename T, uint8_t precision, uint8_t sparse_precision>
double
hll::hyperloglog<T, precision, sparse_precision>::estimate() const {
  if (sparse) {
    size_t nonzero = merged_temp_list().size();
    return linear_estimate((std::size_t)nonzero);
  } else {
    double e;
    std::size_t non_zeros;
    std::tie(e, non_zeros) = raw_estimate();

    if (e <= 5*(1ul << precision))
      e = e - estimate_bias(e);

    double h;
    if (non_zeros < (1ul << precision))
      h = linear_estimate(non_zeros);
    else
      h = e;

    if (h <= threshold())
      return h;
    else
      return e;
  }
}


template <typename T, uint8_t precision, uint8_t sparse_precision>
double
hll::hyperloglog<T, precision, sparse_precision>::measure_error(
    std::size_t orig_card) const {
  double e;
  std::tie(e, std::ignore) = raw_estimate();
  return e-static_cast<double>(orig_card);
}
