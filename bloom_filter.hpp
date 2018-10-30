#ifndef _BLOOM_FILTER_HPP
#define _BLOOM_FILTER_HPP

#include <algorithm>
#include <array>
#include <sdsl/bit_vectors.hpp>
#include <string>

#include "MurmurHash3.hpp"
#include "kmc_api/kmc_file.h"

// using namespace std;
using namespace sdsl;

static const char RCN[128] = {
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   //  0
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 10
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 20
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 30
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 40
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 50
    0,   0,   0, 0,   0,   'T', 0,   'G', 0,   0,   // 60
    0,   'C', 0, 0,   0,   0,   0,   0,   'N', 0,   // 70
    0,   0,   0, 0,   'A', 0,   0,   0,   0,   0,   // 80
    0,   0,   0, 0,   0,   0,   0,   'T', 0,   'G', // 90
    0,   0,   0, 'G', 0,   0,   0,   0,   0,   0,   // 100
    'N', 0,   0, 0,   0,   0,   'A', 0,   0,   0,   // 110
    0,   0,   0, 0,   0,   0,   0,   0              // 120
};

class BF {

private:
  static const char _rcc(const char &c) { return RCN[c]; }

  static const std::string _rc(const std::string &kmer) {
    std::string rc(kmer);
    transform(rc.begin(), rc.end(), rc.begin(), _rcc);
    reverse(rc.begin(), rc.end());
    return rc;
  }

  const std::string _minrc(const std::string &kmer) const { return min(kmer, _rc(kmer)); }

  uint64_t _get_hash(const std::string &kmer) const {
    std::string k = _minrc(kmer);
    array<uint64_t, 2> hashes;
    MurmurHash3_x64_128(k.c_str(), k.size(), 0,
                        reinterpret_cast<void *>(&hashes));
    return hashes[0];
  }

public:
  BF(const size_t size) : _mode(false), _bf(size, 0) { _size = size; }
  ~BF() {}

  void add_key(const std::string &kmer) {
    std::string uc_kmer(kmer);
    transform(uc_kmer.begin(), uc_kmer.end(), uc_kmer.begin(), ::toupper);
    uint64_t hash = _get_hash(uc_kmer);
    _bf[hash % _size] = 1;
  }

  void add_refkey(const std::string &kmer) {
    std::string uc_kmer(kmer);
    transform(uc_kmer.begin(), uc_kmer.end(), uc_kmer.begin(), ::toupper);
    uint64_t hash = _get_hash(uc_kmer);
    _bf[hash % _size] = 0;
  }

  bool test_key(const std::string &kmer) const {
    std::string uc_kmer(kmer);
    transform(uc_kmer.begin(), uc_kmer.end(), uc_kmer.begin(), ::toupper);
    uint64_t hash = _get_hash(uc_kmer);
    return _bf[hash % _size];
  }

  void switch_mode() {
    _mode = true;
    _brank = rank_support_v<1>(&_bf);
    _counts = int_vector<8>(_brank(_size), 0, 8);
  }

  bool increment(const std::string &kmer) {
    if (!_mode)
      return false;
    std::string uc_kmer(kmer);
    transform(uc_kmer.begin(), uc_kmer.end(), uc_kmer.begin(), ::toupper);
    uint64_t hash = _get_hash(uc_kmer);
    size_t bf_idx = hash % _size;
    if (_bf[bf_idx]) {
      size_t cnts_idx = _brank(bf_idx);
      _counts[cnts_idx] = _counts[cnts_idx] < 250 ? _counts[cnts_idx] + 1 : 250;
    }
    return true;
  }

  bool increment(const std::string &kmer, const uint32 counter) {
    if (!_mode)
      return false;
    std::string uc_kmer(kmer);
    transform(uc_kmer.begin(), uc_kmer.end(), uc_kmer.begin(), ::toupper);
    uint64_t hash = _get_hash(uc_kmer);
    size_t bf_idx = hash % _size;
    if (_bf[bf_idx]) {
      size_t cnts_idx = _brank(bf_idx);
      uint32 value = _counts[cnts_idx] + counter;
      _counts[cnts_idx] = value < 250 ? value : 250;
    }
    return true;
  }

  uint8_t get_count(const std::string &kmer) const {
    if (_mode) {
      std::string uc_kmer(kmer);
      transform(uc_kmer.begin(), uc_kmer.end(), uc_kmer.begin(), ::toupper);
      uint64_t hash = _get_hash(uc_kmer);
      size_t bf_idx = hash % _size;
      if (_bf[bf_idx])
        return _counts[_brank(bf_idx)];
    }
    return 0;
  }

private:
  BF() {}
  const BF &operator=(const BF &other) { return *this; }
  const BF &operator=(const BF &&other) { return *this; }

  bool _mode; // false = write, true = read
  size_t _size;
  bit_vector _bf;
  rank_support_v<1> _brank;
  int_vector<8> _counts;
};

#endif
