#include <vector>
#include <map>
#include <mutex>

#include "MurmurHash3.h"


namespace hll {

  template<typename T>
  uint64_t hash(const T k);

  template <unsigned short int precision=14,
           unsigned short int sparse_precision=24>
  class HyperLogLog {
    static_assert(precision > 3,
        "Precision should be 4 or greater");
    static_assert(precision <= 18,
        "precision should be 18 or less");
    static_assert(sparse_precision <= 58,
        "Sparse precision should be 58 or less");
    static_assert(precision < sparse_precision,
        "Precision should be less than sparse_precision");

    const static std::map<double, double> bias;
    constexpr static unsigned long sparse_list_max
      = (1ul << precision)/sizeof(uint64_t);
    constexpr static unsigned long temporary_list_max
      = sparse_list_max/10;

    mutable std::mutex insert_mutex;

    bool sparse;
    std::vector<uint8_t> dense;
    std::vector<uint64_t> sparse_list;
    std::vector<uint64_t> temporary_list;


    void convert_to_dense();
    std::vector<uint64_t> merged_temp_list() const;

    double estimate_bias(double biassed_estimate) const;

    std::pair<uint64_t, uint8_t> get_hash_rank(uint64_t hash) const;

    std::pair<uint64_t, uint8_t> decode_hash(uint64_t hash) const;
    uint64_t encode_hash(uint64_t index, uint8_t rank) const;

    std::pair<double, unsigned> raw_estimate() const;
    double linear_estimate(unsigned non_zero) const;

    constexpr double threshold() const;
    constexpr double alpha() const;

  public:
    HyperLogLog(bool create_dense=false);
    HyperLogLog(const hll::HyperLogLog<precision, sparse_precision> &other);

    template <typename T> void insert(T item);
    void merge(const HyperLogLog<precision, sparse_precision>& other);

    double estimate() const;
    double measure_error(unsigned long original_cardinality) const;
  };
}


#include "hyperloglog.tpp"
