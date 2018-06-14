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
#include <mutex>
#include <type_traits>

#define MURMUR_SEED (uint32_t)161803 // golden ratio digits
#define RANK_BITS 6 // log2(64)

#include <iostream>


namespace hll {
#define _hll_define_arithmetic_hash(_Tp)  \
  template<>                    \
  uint64_t hash(const _Tp k) {            \
    uint64_t hash_out[2];                 \
    MurmurHash3_x64_128(&k, sizeof(k),    \
        MURMUR_SEED, hash_out);           \
    return hash_out[1];                   \
  };

  _hll_define_arithmetic_hash(bool)
  _hll_define_arithmetic_hash(char)
  _hll_define_arithmetic_hash(signed char)
  _hll_define_arithmetic_hash(unsigned char)
  _hll_define_arithmetic_hash(wchar_t)
  _hll_define_arithmetic_hash(char16_t)
  _hll_define_arithmetic_hash(char32_t)
  _hll_define_arithmetic_hash(short)
  _hll_define_arithmetic_hash(int)
  _hll_define_arithmetic_hash(long)
  _hll_define_arithmetic_hash(long long)
  _hll_define_arithmetic_hash(unsigned short)
  _hll_define_arithmetic_hash(unsigned int)
  _hll_define_arithmetic_hash(unsigned long)
  _hll_define_arithmetic_hash(unsigned long long)

  template<>
  uint64_t hash<::std::string>(const ::std::string k) {
    uint64_t hash_out[2];
    MurmurHash3_x64_128(k.c_str(), (int)k.length()+1,
        MURMUR_SEED, hash_out);
    return hash_out[1];
  }

  const std::map<double, double> biases[15] = {
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
}

template <unsigned short p, unsigned short sp>
hll::HyperLogLog<p, sp>::HyperLogLog(bool create_dense) {
  if (create_dense) {
    sparse = false;
    convert_to_dense();
  } else {
    sparse = true;
    temporary_list.reserve(temporary_list_max);
  }

}

template <unsigned short p, unsigned short sp>
hll::HyperLogLog<p, sp>::HyperLogLog(const hll::HyperLogLog<p, sp>& other)
  : sparse((std::lock_guard<std::mutex>(other.insert_mutex),
        other.sparse))
    , dense((std::lock_guard<std::mutex>(other.insert_mutex),
        other.dense))
    , sparse_list((std::lock_guard<std::mutex>(other.insert_mutex),
        other.sparse_list))
    , temporary_list((std::lock_guard<std::mutex>(other.insert_mutex),
        other.temporary_list)) {
    }


template <unsigned short p, unsigned short sp>
double hll::HyperLogLog<p, sp>::estimate_bias(double est) const {
  constexpr unsigned int k = 6; // K-nn parameter
  std::vector<std::pair<double, double>> keys;

  std::transform(bias.begin(), bias.end(), std::back_inserter(keys),
      [est] (std::pair<double, double> b) -> std::pair<double, double> {
      return std::make_pair(std::abs(b.first - est), b.second);
      });
  std::partial_sort(keys.begin(), keys.begin()+k, keys.end());

  double sum = 0;
  double weight_sum = 0;
  std::tie(sum, weight_sum) = std::accumulate(keys.begin(), keys.begin()+k,
      std::make_pair(0.0, 0.0),
      [] (
        std::pair<double, double> t,  // (sum, sum of weights)
        std::pair<double, double> key // (distance, bias)
        ) -> std::pair<double, double> {
      return std::make_pair(
        t.first + key.second*1.0/key.first,
        t.second + 1.0/key.first);
      });
  return sum/weight_sum;
}


template <unsigned short p, unsigned short sp>
const std::map<double, double> hll::HyperLogLog<p, sp>::bias = hll::biases[p-4];


template <unsigned short p, unsigned short sp>
constexpr double hll::HyperLogLog<p, sp>::alpha() const {
  if (p == 4)
    return 0.673;
  else if (p == 5)
    return 0.697;
  else if (p == 6)
    return 0.709;
  else
    return 0.7213/(1.0+1.079/(1ul << p));
}

template <unsigned short p, unsigned short sp>
constexpr double hll::HyperLogLog<p, sp>::threshold() const {
  double thresholds[] = { 10,     20,     40,     80,     220,
                          400,    900,    1800,   3100,   6500,
                          11500,  20000,  50000,  120000, 350000};
  return thresholds[p-4];
}

template <unsigned short p, unsigned short sp>
double hll::HyperLogLog<p, sp>::linear_estimate(unsigned non_zero) const {
  double m;
  if (sparse)
    m = (double)(1ul << sp);
  else
    m = (double)(1ul << p);

  return m*std::log(m/(m-non_zero));
}

template <unsigned short p, unsigned short sp>
std::pair<uint64_t, uint8_t>
hll::HyperLogLog<p, sp>::get_hash_rank(uint64_t hash) const {
  unsigned short precision;
  if (sparse)
    precision = sp;
  else
    precision = p;

  uint64_t index = (uint64_t)(hash >> (sizeof(hash)*8 - precision));
  uint8_t rank = (uint8_t)(std::min(sizeof(hash)*8 - precision,
      (unsigned long)(__builtin_clzl(hash << precision)) + 1));
  return std::make_pair(index, rank);
}

template <unsigned short p, unsigned short sp>
uint64_t
hll::HyperLogLog<p, sp>::encode_hash(uint64_t index, uint8_t rank) const {
  return (index << RANK_BITS) | rank;
}

template <unsigned short p, unsigned short sp>
std::pair<uint64_t, uint8_t>
hll::HyperLogLog<p, sp>::decode_hash(uint64_t hash) const {
  uint64_t index = (hash >> RANK_BITS);
  uint8_t rank = (uint8_t)(((1 << RANK_BITS)-1) & hash); // first 6 bits
  return std::make_pair(index, rank);
}


template <unsigned short precision, unsigned short sparse_precision>
template <typename T>
void hll::HyperLogLog<precision, sparse_precision>::insert(T item) {

  uint64_t hash = hll::hash(item);
  uint64_t index;
  uint8_t rank;
  std::tie(index, rank) = get_hash_rank(hash);

  std::lock_guard<std::mutex> lock(insert_mutex);

  if (sparse) {
    uint64_t encoded = encode_hash(index, rank);
    temporary_list.push_back(encoded);

    if (temporary_list.size() >= temporary_list_max) {
      std::vector<uint64_t> new_sparse_list = merged_temp_list();
      sparse_list.swap(new_sparse_list);
      temporary_list.clear();
    }

    if (sparse_list.size() >= sparse_list_max)
      convert_to_dense();

  } else
    if (rank > dense[index]) dense[index] = rank;
}


template <unsigned short precision, unsigned short sparse_precision>
std::vector<uint64_t>
hll::HyperLogLog<precision, sparse_precision>::merged_temp_list() const{
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

  std::vector<uint64_t> new_sparse_list;

  auto it1 = temporary_list_copy.begin();
  auto it2 = sparse_list.begin();

  int i = 0;
  while (it1 != temporary_list_copy.end() and it2 != sparse_list.end()) {
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


  if (it1 == temporary_list_copy.end() && it2 != sparse_list.end())
    new_sparse_list.insert(new_sparse_list.end(),
        it2, sparse_list.end());
  else if (it2 == sparse_list.end() && it1 != temporary_list_copy.end())
    new_sparse_list.insert(new_sparse_list.end(),
        it1, temporary_list_copy.end());

  return new_sparse_list;
}

template <unsigned short precision, unsigned short sparse_precision>
void hll::HyperLogLog<precision, sparse_precision>::convert_to_dense() {

  std::vector<uint64_t> new_sparse_list = merged_temp_list();
  sparse_list.swap(new_sparse_list);
  temporary_list.clear();
  temporary_list.shrink_to_fit();

  dense.resize(1ul << precision, (uint8_t)0);
  for (const auto i: sparse_list) {
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
        (((uint8_t)__builtin_clzl(betweens) -
          (sizeof(betweens)*8 - (sparse_precision - precision))) + 1);
    if (dense_rank > dense[dense_index]) dense[dense_index] = dense_rank;
  }
  sparse = false;
  sparse_list.clear();
  sparse_list.shrink_to_fit();
}



template <unsigned short precision, unsigned short sparse_precision>
std::pair<double, unsigned>
hll::HyperLogLog<precision, sparse_precision>::raw_estimate() const {
  if (sparse)
    throw std::logic_error(
        "`raw_estimate()` does not work with sparse representation.");
  else {
    double sum = 0;
    unsigned non_zeros = 0;
    for(const auto& m_j: dense) {
      if (m_j > 0) non_zeros += 1;
      sum += std::pow(2, -m_j);
    }

    return std::make_pair(
        alpha()*std::pow((1ul << precision), 2)/sum,
        non_zeros);
  }
}

template <unsigned short precision, unsigned short sparse_precision>
double hll::HyperLogLog<precision, sparse_precision>::estimate() const {
  if (sparse) {
    size_t nonzero = merged_temp_list().size();
    return linear_estimate((unsigned)nonzero);
  } else {
    double e;
    unsigned non_zeros;
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


template <unsigned short p, unsigned short sp>
double hll::HyperLogLog<p, sp>::measure_error(unsigned long orig_card) const {
  double e;
  std::tie(e, std::ignore) = raw_estimate();
  return e-(double)orig_card;
}
