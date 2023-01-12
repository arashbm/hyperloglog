#ifndef INCLUDE_HLL_HYPERLOGLOG_HPP_
#define INCLUDE_HLL_HYPERLOGLOG_HPP_

#include <vector>
#include <map>
#include <cstdint>

namespace hll {
  template<typename T>
  struct hash {
    std::uint64_t operator()(const T&, std::uint64_t seed) const;
  };

  template <
    typename T,
    std::uint8_t precision = 14,
    std::uint8_t sparse_precision = 24>
  class hyperloglog {
    static_assert(precision > 3,
        "Precision should be 4 or greater");
    static_assert(precision <= 18,
        "precision should be 18 or less");
    static_assert(sparse_precision <= 58,
        "Sparse precision should be 58 or less");
    static_assert(precision < sparse_precision,
        "Precision should be less than sparse_precision");

  public:
    static constexpr std::uint8_t dense_prec = precision;
    static constexpr std::uint8_t sparse_prec = sparse_precision;

    explicit hyperloglog(
        bool create_dense = false,
        std::uint64_t seed = 0x9E3779B97F4A7C15);

    void insert(T item);
    void merge(const hll::hyperloglog<T, precision, sparse_precision>& other);

    bool is_sparse() const;
    double estimate() const;
    double measure_error(std::size_t original_cardinality) const;
    const std::vector<std::uint8_t>& dense_vec() const;

  private:
    constexpr static std::size_t sparse_list_max
      = (1ul << precision)/sizeof(std::uint64_t);
    constexpr static std::size_t temporary_list_max
      = sparse_list_max/10;
    constexpr static int rank_bits = 6;  // == log2(64)

    bool sparse;
    std::uint64_t seed;
    std::vector<std::uint8_t> dense;
    std::vector<std::uint64_t> sparse_list;
    std::vector<std::uint64_t> temporary_list;

    std::vector<std::uint8_t> converted_to_dense() const;
    void convert_to_dense();
    void merge_temp();

    std::vector<std::uint64_t> merged_temp_list() const;
    std::vector<std::uint64_t> merged_sorted_list(
        const std::vector<std::uint64_t> other) const;

    double estimate_bias(double biassed_estimate) const;

    std::pair<std::uint64_t, std::uint8_t>
    get_hash_rank(std::uint64_t hash) const;

    std::pair<std::uint64_t, std::uint8_t>
    decode_hash(std::uint64_t hash) const;

    std::uint64_t encode_hash(uint64_t index, uint8_t rank) const;

    std::pair<double, std::size_t> raw_estimate() const;
    double linear_estimate(std::size_t non_zero) const;

    constexpr double threshold() const;
    constexpr double alpha() const;
  };
}  // namespace hll


#include "../../src/hyperloglog.tpp"

#endif  // INCLUDE_HLL_HYPERLOGLOG_HPP_
