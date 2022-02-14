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

#ifndef _KMAP_HPP_
#define _KMAP_HPP_

#include <string>
#include <unordered_map>

using namespace std;

// static const char RCN[128] = {
//     0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   //  0
//     0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 10
//     0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 20
//     0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 30
//     0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 40
//     0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 50
//     0,   0,   0, 0,   0,   'T', 0,   'G', 0,   0,   // 60
//     0,   'C', 0, 0,   0,   0,   0,   0,   'N', 0,   // 70
//     0,   0,   0, 0,   'A', 0,   0,   0,   0,   0,   // 80
//     0,   0,   0, 0,   0,   0,   0,   'T', 0,   'G', // 90
//     0,   0,   0, 'G', 0,   0,   0,   0,   0,   0,   // 100
//     'N', 0,   0, 0,   0,   0,   'A', 0,   0,   0,   // 110
//     0,   0,   0, 0,   0,   0,   0,   0              // 120
// };

struct KMAP
{
  unordered_map<string, int> kmers;

  KMAP() {}

  ostream &operator>>(ostream &stream)
  {
    size_t size = kmers.size();
    stream.write(reinterpret_cast<const char *>(&size), sizeof(size_t));
    for (const auto &pair : kmers)
    {
      string::size_type l = pair.first.length();
      stream.write(reinterpret_cast<const char *>(&l), sizeof(string::size_type));
      stream.write(reinterpret_cast<const char *>((char *)pair.first.data()), l);
      stream.write(reinterpret_cast<const char *>(&pair.second), sizeof(int));
    }
    return stream;
  }

  istream &operator<<(istream &stream)
  {
    size_t size;
    stream.read(reinterpret_cast<char *>(&size), sizeof(size_t));
    for (size_t i = 0; i < size; ++i)
    {
      string::size_type l;
      string k;
      int v;
      stream.read(reinterpret_cast<char *>(&l), sizeof(string::size_type));
      k.resize(l);
      stream.read(reinterpret_cast<char *>((char *)k.data()), l);
      stream.read(reinterpret_cast<char *>(&v), sizeof(int));
      kmers[k] = v;
    }
    return stream;
  }

  static const char _compl(const char &c) { return RCN[c]; }

  string canonical(const char *kmer)
  {
    uint k = strlen(kmer);
    char ckmer[k + 1];
    strcpy(ckmer, kmer);
    transform(ckmer, ckmer + k, ckmer, _compl);
    reverse(ckmer, ckmer + k);
    if (strcmp(kmer, ckmer) < 0)
      memmove(ckmer, kmer, k);
    string kmer_string(ckmer);
    return kmer_string;
  }

  bool test_key(const char *kmer)
  {
    string ckmer = canonical(kmer);
    if (kmers.find(ckmer) == kmers.end())
      return false;
    else
      return true;
  }

  void add_key(const char *kmer)
  {
    string ckmer = canonical(kmer);
    kmers[ckmer] = 0;
  }

  void increment(const char *kmer, int counter)
  {
    string ckmer = canonical(kmer);
    if (kmers.find(ckmer) != kmers.end())
    {
      uint32 new_value = kmers[ckmer] + counter;
      kmers[ckmer] = new_value;
    }
  }

  int get_count(const char *kmer)
  {
    string ckmer = canonical(kmer);
    if (kmers.find(ckmer) != kmers.end())
      return kmers[ckmer];
    else
      return 0;
  }
};

#endif