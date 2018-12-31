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
  std::unordered_map<std::string, int> _times;

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

  void increment(const char* kmer, int counter) {
    std::string ckmer = canonical(kmer);
    if(kmers.find(ckmer) != kmers.end()) {
      uint32 new_value = kmers[ckmer] + counter;
      kmers[ckmer] = new_value < 250 ? new_value : 250;
      ++_times[ckmer];
    }
  }

  void increment_with_average(const char* kmer, int counter) {
    std::string ckmer = canonical(kmer);
    if(kmers.find(ckmer) != kmers.end()) {
      uint32 new_value = (kmers[ckmer] * _times[ckmer] + counter) / (_times[ckmer]+1);
      kmers[ckmer] = new_value < 250 ? new_value : 250;
      ++_times[ckmer];
    }
  }

  int get_count(const char* kmer) {
    std::string ckmer = canonical(kmer);
    if(kmers.find(ckmer) != kmers.end())
      return kmers[ckmer];
    else
      return 0;
  }

  int get_times(const char* kmer) {
    std::string ckmer = canonical(kmer);
    if(kmers.find(ckmer) != kmers.end())
      return _times[ckmer];
    else
      return 0;
  }
};

#endif
