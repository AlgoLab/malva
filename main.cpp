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
#include "kmc_api/kmc_file.h"

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

// Useless debug methods -----------------------------------------------------
void check_kmers(const VK_GROUP &kmers) {
  for (VK_GROUP::const_iterator it = kmers.begin(); it != kmers.end(); ++it) {
    for (const auto &p : it->second) {
      for (const auto k : p.second) {
        if(k.size() != opt::k)
          std::cout << "#" << std::endl;
      }
    }
  }
}

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

/**
 * Method to add kmers to the bloom filter
 **/
void add_kmers_to_bf(BF &bf, const VK_GROUP &kmers) {
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

/**
 * Method to check if variants are well covered.  It returns a map
 * <var_index, {allele_index}> (a set to avoid duplicates).
 **/
std::map<int, std::set<int>> get_well_covered_variants(BF &bf, const VK_GROUP &kmers) {
  std::map<int, std::set<int>> wcvs;
  for(const auto &v : kmers) {
    // For each variant
    bool is_good = false;
    for(const auto &p : v.second) {
      // For each allele of the variant
      for(const auto &kmer : p.second) {
        // For each kmer of the allele of the variant
        uint w = bf.get_count(kmer);
        if(w >= opt::min_coverage) { //&& w <= opt::max_coverage) {
          wcvs[v.first].insert(p.first);
	  is_good = true;
	}
      }
    }
    if(opt::all_variants && !is_good) {
      wcvs[v.first].insert(0);
    }
  }
  return wcvs;
}

int main(int argc, char *argv[]) {
  hts_set_log_level(HTS_LOG_OFF);

  parse_arguments(argc, argv);

  // STEP 1: add VCF kmers to bloom filter
  gzFile fasta_in = gzopen(opt::fasta_path.c_str(), "r");
  kseq_t *reference = kseq_init(fasta_in);
  kseq_read(reference);

  htsFile *vcf = bcf_open(opt::vcf_path.c_str(), "r");
  bcf_hdr_t *vcf_header = bcf_hdr_read(vcf);

  BF bf (opt::bf_size);
  VB vb (reference->seq.s, opt::k);
  int tot_vcf_kmers = 0;

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
      tot_vcf_kmers += count_kmers(kmers);
      add_kmers_to_bf(bf, kmers);
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
    tot_vcf_kmers += count_kmers(kmers);
    add_kmers_to_bf(bf, kmers);
    vb.clear();
  }

  bcf_hdr_destroy(vcf_header);
  bcf_destroy(vcf_record);
  bcf_close(vcf);

  bf.switch_mode();

  BF ref_bf(opt::bf_size);
  if (opt::strict_mode) {
    for (size_t p = (opt::ref_k - opt::k) / 2; p < reference->seq.l - opt::ref_k; ++p) {
      std::string ref_ksub(reference->seq.s + p, opt::k);
      if (bf.test_key(ref_ksub)) {
	std::string context(reference->seq.s + p - ((opt::ref_k - opt::k) / 2), opt::ref_k);
        ref_bf.add_key(context);
      }
    }
  }

  pelapsed("BF creation complete (" + std::to_string(tot_vcf_kmers) + ")");

  // STEP 2: test variants present in read sample
  CKMCFile kmer_db;
  uint32 klen, mode, min_counter, pref_len, sign_len, min_c, counter;
  uint64 tot_kmers, max_c;
  if (!kmer_db.OpenForListing(opt::kmc_sample_path)) {
    std::cerr << "ERROR: cannot open " << opt::kmc_sample_path << std::endl;
    return 1;
  }
  kmer_db.Info(klen, mode, min_counter, pref_len, sign_len, min_c, max_c, tot_kmers);
  CKmerAPI kmer_obj(klen);

  std::string context;
  while(kmer_db.ReadNextKmer(kmer_obj, counter)) {
    kmer_obj.to_string(context);
    if(!ref_bf.test_key(context)) {
      std::string kmer(context.c_str() + ((opt::ref_k - opt::k) / 2), opt::k);
      if (bf.test_key(kmer))
	bf.increment(kmer, counter);
    }
  }

  /**
  // OLD - without context
  std::string kmer;
  while(kmer_db.ReadNextKmer(kmer_obj, counter)) {
    kmer_obj.to_string(kmer);
    if (bf.test_key(kmer))
      bf.increment(kmer, counter);
  }
  **/

  pelapsed("BF weights created");

  // STEP 3: check if variants in vcf are covered enough
  vcf = bcf_open(opt::vcf_path.c_str(), "r");
  vcf_header = bcf_hdr_read(vcf);

  kstring_t htxt = {0, 0, 0};
  bcf_hdr_format(vcf_header, 0, &htxt);
  cout << htxt.s;
  free(htxt.s);

  vcf_record = bcf_init();
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
       * 2. check if variants are covered 
       * 3. output covered variants
       * 4. clear block
       ***/
      VK_GROUP kmers = vb.extract_kmers();
      std::map<int, std::set<int>> well_covered_variants = get_well_covered_variants(bf, kmers);
      vb.output_variants(well_covered_variants);
      vb.clear();
    }
    vb.add_variant(v);
  }
  if(!vb.empty()) {
    /***
     * 1. extract k-mers
     * 2. check if variants are covered 
     * 3. output covered variants
     * 4. clear block
     ***/
    VK_GROUP kmers = vb.extract_kmers();
    std::map<int, std::set<int>> well_covered_variants = get_well_covered_variants(bf, kmers);
    vb.output_variants(well_covered_variants);
    vb.clear();
  }

  bcf_hdr_destroy(vcf_header);
  bcf_destroy(vcf_record);
  bcf_close(vcf);

  kseq_destroy(reference);
  gzclose(fasta_in);

  cout.flush();

  pelapsed("Execution completed");
  
  return 0;
}
