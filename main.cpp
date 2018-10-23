#include <fstream>
#include <iostream>
#include <utility>
#include <algorithm>
#include <chrono>

#include <string>
#include <vector>
#include <map>
#include <set>

#include <zlib.h>
#include <math.h>

#include "kseq.h"
#include "hts_log.h"
#include "vcf.h"

#include "argument_parser.hpp"
#include "var_block.hpp"
#include "bloom_filter.hpp"
// #include "string_3bits.hpp"

auto start_t = chrono::high_resolution_clock::now();

void pelapsed(const string &s = "") {
  auto now_t = chrono::high_resolution_clock::now();
  cerr << "[malva/" << s << "] Time elapsed "
       << chrono::duration_cast<chrono::milliseconds>(now_t - start_t).count()
       << endl;
}

KSEQ_INIT(gzFile, gzread)

// Useless methods -----------------------------------------------------------
void print_kmers(const VK_GROUP &kmers) {
  for (VK_GROUP::const_iterator it = kmers.begin(); it != kmers.end(); ++it) {
    int vID = it->first;
    std::cout << vID << std::endl;
    for (const auto &p : it->second) {
      int altID = p.first;
      std::cout << " " << altID << ": ";
      for (const auto k : p.second) {
        std::cout << k << " ";
      }
      std::cout << std::endl;
    }
  }
}

int count_kmers(const VK_GROUP &kmers) {
  int N = 0;
  for (VK_GROUP::const_iterator it = kmers.begin(); it != kmers.end(); ++it) {
    for (const auto &p : it->second) {
      N += p.second.size();
    }
  }
  return N;
}
// ---------------------------------------------------------------------------

void add_kmers_to_bf(VK_GROUP &kmers, BF &bf) {
  for(const auto &v : kmers) {
    // For each variant
    for(const auto &p : v.second) {
      // For each allele of the variant
      for(const auto &kmer : p.second) {
        // For each kmer of the allele of the variant
        bf.add_key(kmer);
      }
    }
  }
}

int main(int argc, char *argv[]) {
  hts_set_log_level(HTS_LOG_OFF);

  parse_arguments(argc, argv);

  gzFile fasta_in = gzopen(opt::fasta_path.c_str(), "r");
  kseq_t *reference = kseq_init(fasta_in);
  kseq_read(reference);

  htsFile *vcf = bcf_open(opt::vcf_path.c_str(), "r");
  bcf_hdr_t *vcf_header = bcf_hdr_read(vcf);

  // STEP 1: add VCF kmers to bloom filter
  BF bf (opt::bf_size);
  VB vb (reference->seq.s, opt::k);
  int tot_kmers = 0;

  bcf1_t *vcf_record = bcf_init();
  while (bcf_read(vcf, vcf_header, vcf_record) == 0) {
    bcf_unpack(vcf_record, BCF_UN_STR);

    Variant v (vcf_header, vcf_record);
    if(!v.is_good)
      continue;

    if(vb.empty()) {
      vb.add_variant(v);
      continue;
    }

    if(!vb.is_near_to_last(v)) {
      /***
       * 1. extract k-mers
       * 2. add k-mers to BF
       * 3. clear block
       ***/
      VK_GROUP kmers = vb.extract_kmers();
      tot_kmers += count_kmers(kmers);
      // print_kmers(kmers);
      add_kmers_to_bf(kmers, bf);
      vb.clear();
    }
    vb.add_variant(v);
  }
  if(!vb.empty()) {
    /***
     * 1. extract k-mers
     * 2. add k-mers to BF
     * 3. clear block
     ***/
    VK_GROUP kmers = vb.extract_kmers();
    tot_kmers += count_kmers(kmers);
    // print_kmers(kmers);
    add_kmers_to_bf(kmers, bf);
    vb.clear();
  }

  pelapsed("BF creation complete (" + std::to_string(tot_kmers) + ")");

  kseq_destroy(reference);
  gzclose(fasta_in);
  bcf_hdr_destroy(vcf_header);
  bcf_destroy(vcf_record);
  bcf_close(vcf);

  return 0;
}
