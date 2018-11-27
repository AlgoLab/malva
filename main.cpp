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

float binomial(const int &x, const int &y) {
  return tgamma((float)x+1.0)/(tgamma((float)y+1.0)*tgamma((float)x-(float)y+1.0));
}

GT select_genotype(const Variant &v, const std::vector<int> covs, const float &error_rate) {
  if(covs.size() == 1)
    // The variant wasn't present in any sample: we have only the
    // coverage of the reference allele
    return std::make_pair("0/0", 1);

  float max_prob = 0.0;
  std::string best_geno = "0/0";

  for(uint g1 = 0; g1 < covs.size(); ++g1) {
    for(uint g2 = g1; g2 < covs.size(); ++g2) {
      float prior;
      float posterior;
      int total_sum = accumulate(covs.begin(), covs.end(), 0);
      if(g1 == g2) {
        prior = std::pow(v.frequencies[g1], 2);
        int truth = covs[g1];
        int error = total_sum - truth;
        
        posterior = binomial(truth + error, truth) * pow(1-error_rate, truth) * pow(error_rate, error);
      } else {
        prior = 2 * v.frequencies[g1] * v.frequencies[g2];
        int truth1 = covs[g1];
        int truth2 = covs[g2];
        int error = total_sum - truth1 - truth2;
        posterior = binomial(truth1 + truth2 + error, truth1 + truth2) * binomial(truth1 + truth2, truth1) * pow((1-error_rate)/2, truth1) * pow((1-error_rate)/2, truth2) * pow(error_rate, error);;
      }

      float prob = prior * posterior;
      if(prob>max_prob) {
        max_prob = prob;
        best_geno = to_string(g1) + "/" + to_string(g2);
      }
    }
  }
  return std::make_pair(best_geno, max_prob);
}

std::vector<GT> genotype(BF &bf, const VB &vb, const VK_GROUP &kmers, const int &cov, const float &error_rate) {
  std::vector<GT> genotypes;
  for(const auto &var : kmers) {
    // For each variant
    Variant v = vb.get_variant(var.first);
    std::vector<int> allele_covs (var.second.size()+1, 0); //+1 for the reference allele
    int i = 0;
    for(const auto &p : var.second) {
      // For each allele of the variant
      int all_cov = 0;
      for(const auto &Ks : p.second) {
        // For each list of kmers of the allele
        /**
         * Here we can have:
         * - a single list of multiple kmers (ie allele is longer than k)
         * - multiple (>=1) lists of length 1
         * In both the cases, we take as allele coverage the mean of the
         * coverages
         **/
        if(Ks.size() == 1) {
          uint w = bf.get_count(Ks[0].c_str());
          if(w != 0) {
            // If we have only a single kmer, only if it is covered,
            // we can use its coverage. Otherwise, we do not consider
            // it at all (possible there will be other kmers)
            if(all_cov == 0)
              all_cov = w;
            else
              all_cov = round(all_cov + w)/2;
          }
        } else {
          // We have more kmers: we iterate to compute the average
          // coverage
          for(const auto &kmer : Ks) {
            uint w = bf.get_count(kmer.c_str());
            if(all_cov == 0)
              all_cov = w;
            else
              all_cov = round(all_cov + w)/2;
          }
        }
      }
      // we can now set the allele coverage
      allele_covs[++i] = all_cov;
    }
    // We now compute the coverage of the reference allele as
    // (expected_coverage - average_coverage_among_alt_alleles)
    // We need the -1 here because the first element of the vector is
    // 0, reserved for reference allele
    float average_cov_alts = 0;
    if(allele_covs.size() > 1)
      average_cov_alts = accumulate(allele_covs.begin(), allele_covs.end(), 0.0)/(allele_covs.size()-1);
    allele_covs[0] = cov - average_cov_alts;

    GT gt = select_genotype(v, allele_covs, error_rate);
    genotypes.push_back(gt);
  }
  return genotypes;
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

    Variant v (vcf_header, vcf_record, opt::pop);

    // We do not consider variants with <CN> or not present in
    // considered samples, i.e. is 0|0 for all samples
    if(!v.has_alts or !v.is_present)
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
  uint32 kmer_sum = 0;
  uint32 klen, mode, min_counter, pref_len, sign_len, min_c, counter;
  uint64 tot_kmers, max_c;
  kmer_db.Info(klen, mode, min_counter, pref_len, sign_len, min_c, max_c, tot_kmers);
  CKmerAPI kmer_obj(klen);

  char context[opt::ref_k+1];
  while(kmer_db.ReadNextKmer(kmer_obj, counter)) {
    kmer_sum += counter;
    kmer_obj.to_string(context);
    transform(context, context + opt::ref_k, context, ::toupper);
    if(!ref_bf.test_key(context)) {
      char kmer[opt::k+1];
      strncpy(kmer, context + ((opt::ref_k - opt::k) / 2), opt::k);
      kmer[opt::k] = '\0';
      bf.increment_with_average(kmer, counter);
    }
  }

  uint32 kmer_cov = kmer_sum / tot_kmers;
  uint32 cov = (kmer_cov * opt::read_len) / ((opt::read_len - opt::ref_k + 1) * (1 - opt::error_rate * opt::ref_k));
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

    Variant v (vcf_header, vcf_record, opt::pop);
    // In this step, we must consider the variants not present in the
    // samples: their genotype is 0/0
    if(!v.has_alts)
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
      std::vector<GT> gts = genotype(bf, vb, kmers, cov, opt::error_rate);
      vb.output_variants(gts);
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
    std::vector<GT> gts = genotype(bf, vb, kmers, cov, opt::error_rate);
    vb.output_variants(gts);
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
