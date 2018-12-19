#ifndef _KMAP_HPP_
#define _KMAP_HPP_

#include <string>
#include <unordered_map>

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

struct KMAP {
  std::unordered_map<std::string, int> kmers;

  KMAP() {}

  static const char _compl(const char &c) { return RCN[c]; }

  std::string canonical(const char* kmer) {
    uint k = strlen(kmer);
    char ckmer[k + 1];
    strcpy(ckmer, kmer);
    std::transform(ckmer, ckmer + k, ckmer, _compl);
    std::reverse(ckmer, ckmer + k);
    if (strcmp(kmer, ckmer) < 0)
      memmove(ckmer, kmer, k);
    std::string kmer_string (ckmer);
    return kmer_string;
  }

  bool test_key(const char* kmer) {
    std::string ckmer = canonical(kmer);
    if(kmers.find(ckmer) == kmers.end())
      return false;
    else
      return true;
  }

  void add_key(const char* kmer) {
    std::string ckmer = canonical(kmer);
    kmers[ckmer] = 0;
  }

  void increment_with_average(const char* kmer, int counter) {
    std::string ckmer = canonical(kmer);
    if(kmers.find(ckmer) != kmers.end()) {
      uint32 old_value = kmers[ckmer];
      uint32 new_value;
      if (old_value == 0)
        new_value = counter;
      else
        new_value = round((old_value + counter) / 2);
      kmers[ckmer] = new_value < 250 ? new_value : 250;
    }
  }

  int get_count(const char* kmer) {
    std::string ckmer = canonical(kmer);
    if(kmers.find(ckmer) != kmers.end())
      return kmers[ckmer];
    else
      return 0;
  }
};

#endif
