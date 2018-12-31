#ifndef _BLOOM_FILTER_HPP
#define _BLOOM_FILTER_HPP

#include <algorithm>
#include <array>
#include <cstring>
#include <sdsl/bit_vectors.hpp>

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
  static const char _compl(const char &c) { return RCN[c]; }

  void _canonical(const char *kmer, char *ckmer, const int &k) const {
    strcpy(ckmer, kmer);
    std::transform(ckmer, ckmer + k, ckmer, _compl);
    std::reverse(ckmer, ckmer + k);
    if (strcmp(kmer, ckmer) < 0)
      memmove(ckmer, kmer, k);
  }

  uint64_t _get_hash(const char *kmer) const {
    uint k = strlen(kmer);
    char ckmer[k + 1];
    _canonical(kmer, ckmer, k);
    array<uint64_t, 2> hashes;
    MurmurHash3_x64_128(ckmer, k, 0, reinterpret_cast<void *>(&hashes));
    return hashes[0];
  }

public:
  BF(const size_t size) : _mode(false), _bf(size, 0) { _size = size; }
  ~BF() {}

  void add_key(const char *kmer) {
    uint64_t hash = _get_hash(kmer);
    _bf[hash % _size] = 1;
  }

  void add_refkey(const char *kmer) {
    uint64_t hash = _get_hash(kmer);
    _bf[hash % _size] = 0;
  }

  bool test_key(const char *kmer) const {
    uint64_t hash = _get_hash(kmer);
    return _bf[hash % _size];
  }

  int get_times(const char *kmer) const {
    uint64_t hash = _get_hash(kmer);
    size_t bf_idx = hash % _size;
    size_t cnts_idx = _brank(bf_idx);
    return _times[cnts_idx];
  }

  void switch_mode() {
    _mode = true;
    _brank = rank_support_v<1>(&_bf);
    _counts = int_vector<8>(_brank(_size), 0, 8);
    _times = int_vector<8>(_brank(_size), 0, 8);
  }

  bool increment(const char *kmer, const uint32 counter) {
    if (!_mode)
      return false;
    uint64_t hash = _get_hash(kmer);
    size_t bf_idx = hash % _size;
    if (_bf[bf_idx]) {
      size_t cnts_idx = _brank(bf_idx);
      uint32 new_value = _counts[cnts_idx] + counter;
      _counts[cnts_idx] = new_value < 250 ? new_value : 250;
      ++_times[cnts_idx];
    }
    return true;
  }

    bool increment_with_average(const char *kmer, const uint32 counter) {
    if (!_mode)
      return false;
    uint64_t hash = _get_hash(kmer);
    size_t bf_idx = hash % _size;
    if (_bf[bf_idx]) {
      size_t cnts_idx = _brank(bf_idx);
      uint32 new_value = (_counts[cnts_idx] * _times[cnts_idx] + counter) / (_times[cnts_idx]+1);
      _counts[cnts_idx] = new_value < 250 ? new_value : 250;
      ++_times[cnts_idx];
    }
    return true;
  }

  uint8_t get_count(const char *kmer) const {
    if (_mode) {
      uint64_t hash = _get_hash(kmer);
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
  int_vector<8> _times;
};

#endif
