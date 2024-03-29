/**
 * MALVA - genotyping by Mapping-free ALternate-allele detection of known VAriants
 * Copyright (C) 2019  Giulia Bernardini, Luca Denti, Marco Previtali
 *
 * This file is part of MALVA.
 *
 * MALVA is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MALVA is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MALVA; see the file LICENSE. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#ifndef _BLOOM_FILTER_HPP
#define _BLOOM_FILTER_HPP

#include <algorithm>
#include <array>
#include <cstring>
#include <sdsl/bit_vectors.hpp>

#include "xxhash.h"
#include "kmc_file.h"

using namespace std;
using namespace sdsl;

static const char RCN[128] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     //  0
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     // 10
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     // 20
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     // 30
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     // 40
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     // 50
    0, 0, 0, 0, 0, 'T', 0, 'G', 0, 0, // 60
    0, 'C', 0, 0, 0, 0, 0, 0, 'N', 0, // 70
    0, 0, 0, 0, 'A', 0, 0, 0, 0, 0,   // 80
    0, 0, 0, 0, 0, 0, 0, 'T', 0, 'G', // 90
    0, 0, 0, 'G', 0, 0, 0, 0, 0, 0,   // 100
    'N', 0, 0, 0, 0, 0, 'A', 0, 0, 0, // 110
    0, 0, 0, 0, 0, 0, 0, 0            // 120
};

class BF
{

private:
  static const char _compl(const char &c) { return RCN[c]; }

  void _canonical(const char *kmer, char *ckmer, const int &k) const
  {
    strcpy(ckmer, kmer);
    transform(ckmer, ckmer + k, ckmer, _compl);
    reverse(ckmer, ckmer + k);
    if (strcmp(kmer, ckmer) < 0)
      memmove(ckmer, kmer, k);
  }

  uint64_t _get_hash(const char *kmer) const
  {
    uint k = strlen(kmer);
    char ckmer[k + 1];
    _canonical(kmer, ckmer, k);
    uint64_t hashes = XXH3_64bits(ckmer, k);
    return hashes;
  }

public:
  BF() : _mode(false), _bf(0, 0) { _size = 0; };
  BF(const size_t size) : _mode(false), _bf(size, 0) { _size = size; }
  ~BF() {}

  void add_key(const char *kmer)
  {
    uint64_t hash = _get_hash(kmer);
    _bf[hash % _size] = 1;
  }

  bool test_key(const char *kmer) const
  {
    uint64_t hash = _get_hash(kmer);
    return _bf[hash % _size];
  }

  void switch_mode()
  {
    _mode = true;
    _brank = rank_support_v<1>(&_bf);
    _counts = int_vector<16>(_brank(_size), 0, 16);
  }

  bool increment(const char *kmer, const uint32 counter)
  {
    if (!_mode)
      return false;
    uint64_t hash = _get_hash(kmer);
    size_t bf_idx = hash % _size;
    if (_bf[bf_idx])
    {
      size_t cnts_idx = _brank(bf_idx);
      uint32 new_value = _counts[cnts_idx] + counter;
      _counts[cnts_idx] = new_value;
    }
    return true;
  }

  uint16_t get_count(const char *kmer) const
  {
    if (_mode)
    {
      uint64_t hash = _get_hash(kmer);
      size_t bf_idx = hash % _size;
      if (_bf[bf_idx])
        return _counts[_brank(bf_idx)];
    }
    return 0;
  }

  ostream &operator>>(ostream &stream)
  {
    stream.write(reinterpret_cast<const char *>(&_mode), sizeof(bool));
    stream.write(reinterpret_cast<const char *>(&_size), sizeof(size_t));
    _bf.serialize(stream);
    // We don't serialize _brank since loading it will crash.
    // TODO: this need some further investigation.
    _counts.serialize(stream);
    return stream;
  }

  istream &operator<<(istream &stream)
  {
    stream.read(reinterpret_cast<char *>(&_mode), sizeof(bool));
    stream.read(reinterpret_cast<char *>(&_size), sizeof(size_t));
    _bf.load(stream);
    _brank = rank_support_v<1>(&_bf);
    _counts.load(stream);
    return stream;
  }

private:
  const BF &operator=(const BF &other) { return *this; }
  const BF &operator=(const BF &&other) { return *this; }

  bool _mode; // false = write, true = read
  size_t _size;
  bit_vector _bf;
  rank_support_v<1> _brank;
  int_vector<16> _counts;
};

#endif
