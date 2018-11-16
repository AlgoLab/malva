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

auto start_t = chrono::high_resolution_clock::now();

void pelapsed(const string &s = "") {
  auto now_t = chrono::high_resolution_clock::now();
  cerr << "[malva/" << s << "] Time elapsed "
       << chrono::duration_cast<chrono::milliseconds>(now_t - start_t).count()
       << endl;
}

KSEQ_INIT(gzFile, gzread)

// Useless debug methods -----------------------------------------------------
// !!! (maybe) Deprecated function !!!
void check_kmers(const VK_GROUP &kmers) {
  for (VK_GROUP::const_iterator it = kmers.begin(); it != kmers.end(); ++it) {
    for (const auto &p : it->second) {
      for (const auto &Ks : p.second) {
	for(const auto &k : Ks) {
	  if(k.size() != opt::k)
	    std::cout << "#" << std::endl;
	}
      }
    }
  }
}

// !!! (maybe) Deprecated function !!!
void print_kmers(const VK_GROUP &kmers) {
  for (VK_GROUP::const_iterator it = kmers.begin(); it != kmers.end(); ++it) {
    int vID = it->first;
    std::cout << vID << std::endl;
    for (const auto &p : it->second) {
      int altID = p.first;
      std::cout << " " << altID << ": ";
      for (const auto &Ks : p.second) {
	for (const auto &k : Ks)
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
      for (const auto &Ks : p.second)
	N += Ks.size();
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
      for(const auto &Ks : p.second) {
	// For each list of kmers of the allele
	for(const auto &kmer : Ks) {
	  // For each kmer in the kmer list
	  bf.add_key(kmer.c_str());
	}
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
      for(const auto &Ks : p.second) {
        // For each list of kmers of the allele
	uint kpt = 0; // kmers passing threshold
	for(const auto &kmer : Ks) {
	  uint w = bf.get_count(kmer.c_str());
	  if(w >= opt::min_coverage) //&& w <= opt::max_coverage)
	    ++kpt;
	}
	/**
	 * !!! By asking that *all* kmers must be covered, we could
	 * miss some alleles longer than k
	 **/
	if (kpt == Ks.size()) {
	  wcvs[v.first].insert(p.first);
	  is_good = true;
	  break; // LD: I think this break is good
	}
      }
    }
    if(opt::all_variants && !is_good)
      wcvs[v.first].insert(0);
  }
  return wcvs;
}

/**
 * Method to clean and print VCF header. It adds GT and GQ FORMAT,
 * removes all samples, and adds donor sample.
 **/
void print_cleaned_header(bcf_hdr_t *vcf_header) {
  // Adding format - if already present, they won't be added
  bcf_hdr_append(vcf_header, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
  bcf_hdr_append(vcf_header, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");

  // Adding donor sample and removing all other samples
  const char *new_sample = "DONOR";
  bcf_hdr_add_sample(vcf_header, new_sample);
  bcf_hdr_sync(vcf_header);
  bcf_hdr_set_samples(vcf_header, new_sample, 0);

  // Formatting and printing header
  kstring_t htxt = {0, 0, 0};
  bcf_hdr_format(vcf_header, 0, &htxt);
  std::cout << htxt.s;
  free(htxt.s);
}

int main(int argc, char *argv[]) {
  hts_set_log_level(HTS_LOG_OFF);

  parse_arguments(argc, argv);

  // STEP 0: open and check input files
  gzFile fasta_in = gzopen(opt::fasta_path.c_str(), "r");
  kseq_t *reference = kseq_init(fasta_in);
  kseq_read(reference);

  htsFile *vcf = bcf_open(opt::vcf_path.c_str(), "r");
  bcf_hdr_t *vcf_header = bcf_hdr_read(vcf);
  bcf1_t *vcf_record = bcf_init();

  CKMCFile kmer_db;
  if (!kmer_db.OpenForListing(opt::kmc_sample_path)) {
    std::cerr << "ERROR: cannot open " << opt::kmc_sample_path << std::endl;
    return 1;
  }

  // STEP 1: add VCF kmers to bloom filter
  BF bf (opt::bf_size);
  VB vb (reference->seq.s, opt::k);
  int tot_vcf_kmers = 0;

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

  pelapsed("BF creation complete (" + std::to_string(tot_vcf_kmers) + ")");

  BF ref_bf(opt::bf_size);
  if (opt::strict_mode) {
    std::string ref_ksub(reference->seq.s + (opt::ref_k-opt::k)/2, opt::k);
    std::string context(reference->seq.s, opt::ref_k);
    transform(ref_ksub.begin(), ref_ksub.end(), ref_ksub.begin(), ::toupper);
    transform(context.begin(), context.end(), context.begin(), ::toupper);
    if (bf.test_key(ref_ksub.c_str())) {
      ref_bf.add_key(context.c_str());
    }
    for(uint p = opt::ref_k; p < reference->seq.l; ++p) {
      char c1 = toupper(reference->seq.s[p]);
      context.erase(0,1);
      context += c1;
      char c2 = toupper(reference->seq.s[p - (opt::ref_k-opt::k)/2]);
      ref_ksub.erase(0,1);
      ref_ksub += c2;
      if (bf.test_key(ref_ksub.c_str()))
        ref_bf.add_key(context.c_str());
    }
  }

  pelapsed("Reference BF creation complete");

  // STEP 2: test variants present in read sample
  uint32 klen, mode, min_counter, pref_len, sign_len, min_c, counter;
  uint64 tot_kmers, max_c;
  kmer_db.Info(klen, mode, min_counter, pref_len, sign_len, min_c, max_c, tot_kmers);
  CKmerAPI kmer_obj(klen);

  char context[opt::ref_k+1];
  while(kmer_db.ReadNextKmer(kmer_obj, counter)) {
    kmer_obj.to_string(context);
    transform(context, context + opt::ref_k, context, ::toupper);
    if(!ref_bf.test_key(context)) {
      char kmer[opt::k+1];
      strncpy(kmer, context + ((opt::ref_k - opt::k) / 2), opt::k);
      kmer[opt::k] = '\0';
      bf.increment(kmer, counter);
    }
  }

  pelapsed("BF weights created");

  // STEP 3: check if variants in vcf are covered enough
  vcf = bcf_open(opt::vcf_path.c_str(), "r");
  vcf_header = bcf_hdr_read(vcf);
  print_cleaned_header(vcf_header);
  bcf_hdr_destroy(vcf_header);
  bcf_close(vcf);

  vcf = bcf_open(opt::vcf_path.c_str(), "r");
  vcf_header = bcf_hdr_read(vcf);
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
