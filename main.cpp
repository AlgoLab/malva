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

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <utility>

#include <unordered_map>
#include <set>
#include <string>
#include <vector>

#include <math.h>
#include <zlib.h>

#include "hts_log.h"
#include "kmc_api/kmc_file.h"
#include "kseq.h"
#include "vcf.h"

#include "argument_parser.hpp"
#include "bloom_filter.hpp"
#include "variant.hpp"
#include "var_block.hpp"
#include "kmap.hpp"

using namespace std;

auto start_t = chrono::high_resolution_clock::now();

void pelapsed(const string &s = "", const bool rollback = false) {
  auto now_t = chrono::high_resolution_clock::now();
  cerr << "[malva-geno/" << s << "] Time elapsed "
	    << chrono::duration_cast<chrono::milliseconds>(now_t - start_t).count()/1000 << "s";
  if(rollback) cerr << "\r";
  else cerr << endl;
}

KSEQ_INIT(gzFile, gzread)

/**
 * Method to compute and store the coverages of the alleles of the
 * variants. It uses the coverages stored in the bloom filters/map.
 **/
void set_coverages(vector<VB> &vbs, BF &bf, KMAP &ref_bf) {
  for(uint vb_idx=0; vb_idx<vbs.size(); ++vb_idx) {
    for(Variant &v : vbs[vb_idx]) {
      for(uint all_idx = 0; all_idx<v.alts.size()+1; ++all_idx) {
	float all_cov = 0.0;
	for(const signature S : v.SIGNS[all_idx]) {
	  float curr_cov = 0.0;
	  int n = 0; // Number of kmers in the signature
	  for(const string &kmer : S) {
	    int w = 0;
	    if(all_idx == 0)
	      w = ref_bf.get_count(kmer.c_str());
	    else
	      w = bf.get_count(kmer.c_str());
	    if(w>0) { //maybe useless
	      curr_cov = (curr_cov * n + w) / (n + 1);
	      ++n;
	    }
	  }
	  if(curr_cov > all_cov)
	    all_cov = curr_cov;
	}
	v.set_coverage(all_idx, all_cov);
      }
    }
  }
}

/**
 * Method to clean and print VCF header. It adds GT and GQ FORMAT,
 * removes all samples, and adds donor sample.
 **/
void print_cleaned_header(bcf_hdr_t *vcf_header) {
  // Adding format fields - if already present, they won't be added
  bcf_hdr_append(vcf_header, "##FORMAT=<ID=GT,Number=1,Type=String,"
                 "Description=\"Genotype\">");
  bcf_hdr_append(vcf_header, "##FORMAT=<ID=GQ,Number=1,Type=Integer,"
                 "Description=\"Genotype Quality\">");

  // Adding donor sample and removing all other samples
  const char *new_sample = "DONOR";
  bcf_hdr_add_sample(vcf_header, new_sample);
  bcf_hdr_sync(vcf_header);
  bcf_hdr_set_samples(vcf_header, new_sample, 0);

  // Formatting and printing header
  kstring_t htxt = {0, 0, 0};
  bcf_hdr_format(vcf_header, 0, &htxt);
  cout << htxt.s;
  free(htxt.s);
}

// ---------------------------------------------------------------------------

unordered_map<string, string> extract_chrom_sequences(kseq_t *reference) {
  unordered_map<string, string> refs;
  int l;
  while ((l = kseq_read(reference)) >= 0) {
    string id = reference->name.s;
    if(id.compare(0,3,"chr") == 0)
      id = id.substr(3);
    string seq (reference->seq.s);
    refs[id] = seq;
  }
  return refs;
}

// We can do better here...
void fill_bf(BF &bf, const vector<signature> &Ss) {
  for(const signature &S : Ss) {
    for(const string &kmer : S) {
      bf.add_key(kmer.c_str());
    }
  }
}

// ...and here
void fill_rbf(KMAP &ref_bf, const vector<signature> &Ss) {
  for(const signature &S : Ss) {
    for(const string &kmer : S) {
      ref_bf.add_key(kmer.c_str());
    }
  }
}

void fill_bfs(vector<VB> &vbs, BF &bf, KMAP &ref_bf) {
  for(uint vb_idx=0; vb_idx<vbs.size(); ++vb_idx) {
    for(const Variant &v : vbs[vb_idx]) {
      fill_rbf(ref_bf, v.SIGNS[0]);
      for(uint all_idx = 1; all_idx<v.alts.size()+1; ++all_idx) {
	fill_bf(bf, v.SIGNS[all_idx]);
      }
    }
  }
}

void build_signatures(vector<VB> &vbs, const unordered_map<string, string> &refs) {
#pragma omp parallel for num_threads (opt::n_threads) shared (vbs, refs)
  for(uint i = 0; i<vbs.size(); ++i) {
    for(Variant &v : vbs[i]) {
      v.fill_genotypes();
      v.free();
    }
    vbs[i].store_signatures(refs.at(vbs[i].front().seq_name));
  }
}

void compute_genotypes(vector<VB> &vbs) {
#pragma omp parallel for num_threads (opt::n_threads) shared (vbs)
  for(uint i = 0; i<vbs.size(); ++i)
    vbs[i].genotype(opt::max_coverage);
}

void print_variants(vector<VB> &vbs) {
  for(uint i = 0; i<vbs.size(); ++i)
    vbs[i].print();
}

int main(int argc, char *argv[]) {
  hts_set_log_level(HTS_LOG_OFF);

  parse_arguments(argc, argv);

  // STEP 0a: open and check input files
  gzFile fasta_in = gzopen(opt::fasta_path.c_str(), "r");
  kseq_t *reference = kseq_init(fasta_in);

  htsFile *vcf = bcf_open(opt::vcf_path.c_str(), "r");
  bcf_hdr_t *vcf_header = bcf_hdr_read(vcf);
  // This functionality is now deprecated
  // int is_file_flag = 0;
  // if(opt::samples != "-")
  //   is_file_flag = 1;
  // int set_samples_code = bcf_hdr_set_samples(vcf_header, opt::samples.c_str(), is_file_flag);
  // if(set_samples_code != 0) {
  //   cerr << "ERROR: VCF samples subset (code: " << set_samples_code << ")" << endl;
  //   return 1;
  // }
  bcf1_t *vcf_record = bcf_init();

  CKMCFile kmer_db;
  if (!kmer_db.OpenForListing(opt::kmc_sample_path)) {
    cerr << "ERROR: cannot open " << opt::kmc_sample_path << endl;
    return 1;
  }

  // STEP 0b: chromosomes are stored in a map
  unordered_map<string, string> refs = extract_chrom_sequences(reference);
  kseq_destroy(reference);
  gzclose(fasta_in);
  pelapsed("Reference processed");

  // STEP 1: add VCF kmers to bloom filter
  BF bf(opt::bf_size);
  KMAP ref_bf;
  BF context_bf(opt::bf_size);

  vector<string> used_seq_names;
  vector<VB> vbs (1, VB(opt::k, opt::error_rate));

  vcf_record->max_unpack = BCF_UN_INFO;
  while (bcf_read(vcf, vcf_header, vcf_record) == 0) {
    bcf_unpack(vcf_record, BCF_UN_INFO); // Here we unpack till INFO field (no FORMAT and SAMPLES)
    Variant v(vcf, vcf_header, vcf_record, opt::freq_key);

    // We do not consider variants with <CN> or not present in
    // considered samples, i.e. a priori frequency of ref allele = 1
    // (i.e. 0|0 for all samples)
    if (!v.has_alts or !v.is_present)
      continue;

    // We store the contig name (we assume VCF to be ordered, i.e.,
    // contigs are not mixed up)
    if(used_seq_names.empty() || v.seq_name.compare(used_seq_names.back()) != 0)
      used_seq_names.push_back(v.seq_name);

    VB *vb = &(vbs.back());
    if(!vb->empty() && (!vb->back().is_rknear_to(v, opt::k) || v.seq_name.compare(used_seq_names.back()) != 0)) {
      if(vbs.size() == 5000) { // we can process the blocks
	build_signatures(vbs, refs);
	fill_bfs(vbs, bf, ref_bf);
	vbs.clear();
      }

      vbs.push_back(VB(opt::k, opt::error_rate));
      if(v.seq_name.compare(used_seq_names.back()) != 0)
	used_seq_names.push_back(v.seq_name);
    }
    // Here we cannot use pointer vb since we could have add a new vb to the vector
    vbs.back().add_variant(v);
  }
  if(!vbs.empty()) {
    build_signatures(vbs, refs);
    fill_bfs(vbs, bf, ref_bf);
    vbs.clear();
    vbs.push_back(VB(opt::k, opt::error_rate)); // for step 3
  }

  bcf_hdr_destroy(vcf_header);
  bcf_destroy(vcf_record);
  bcf_close(vcf);

  bf.switch_mode();

  pelapsed("BF creation complete");

  for(const auto &seq_name : used_seq_names) {
    string reference = refs[seq_name];
    string ref_ksub(reference, (opt::ref_k - opt::k) / 2, opt::k);
    string context(reference, 0, opt::ref_k);
    transform(ref_ksub.begin(), ref_ksub.end(), ref_ksub.begin(), ::toupper);
    transform(context.begin(), context.end(), context.begin(), ::toupper);
    if (bf.test_key(ref_ksub.c_str()))
      context_bf.add_key(context.c_str());
    for (uint p = opt::ref_k; p < reference.size(); ++p) {
      char c1 = toupper(reference[p]);
      context.erase(0, 1);
      context += c1;
      char c2 = toupper(reference[p - (opt::ref_k - opt::k) / 2]);
      ref_ksub.erase(0, 1);
      ref_ksub += c2;
      if (bf.test_key(ref_ksub.c_str()))
        context_bf.add_key(context.c_str());
    }
  }
  pelapsed("Reference BF creation complete");

  context_bf.switch_mode();

  // STEP 2: test variants present in read sample
  uint32 klen, mode, min_counter, pref_len, sign_len, min_c, counter;
  uint64 tot_kmers, max_c;
  kmer_db.Info(klen, mode, min_counter, pref_len, sign_len, min_c, max_c, tot_kmers);
  CKmerAPI kmer_obj(klen);

  char context[opt::ref_k + 1];
  while (kmer_db.ReadNextKmer(kmer_obj, counter)) {
    kmer_obj.to_string(context);
    transform(context, context + opt::ref_k, context, ::toupper);
    char kmer[opt::k + 1];
    strncpy(kmer, context + ((opt::ref_k - opt::k) / 2), opt::k);
    kmer[opt::k] = '\0';
    ref_bf.increment(kmer, counter);
    if (!context_bf.test_key(context)) {
      bf.increment(kmer, counter);
    }
  }
  pelapsed("BF weights created");

  // STEP 3: check if variants in vcf are covered enough
  vcf = bcf_open(opt::vcf_path.c_str(), "r");
  vcf_header = bcf_hdr_read(vcf);
  int set_samples_code = bcf_hdr_set_samples(vcf_header, NULL, 0);
  print_cleaned_header(vcf_header);
  bcf_hdr_destroy(vcf_header);
  bcf_close(vcf);

  vcf = bcf_open(opt::vcf_path.c_str(), "r");
  vcf_header = bcf_hdr_read(vcf);
  // set_samples_code = bcf_hdr_set_samples(vcf_header, opt::samples.c_str(), is_file_flag);
  vcf_record = bcf_init();

  string last_seq_name = "";
  vcf_record->max_unpack = BCF_UN_INFO;
  while (bcf_read(vcf, vcf_header, vcf_record) == 0) {
    bcf_unpack(vcf_record, BCF_UN_INFO); // Here we unpack till INFO field (no FORMAT and SAMPLES)
    Variant v(vcf, vcf_header, vcf_record, opt::freq_key);

    // We do not consider variants with <CN>
    if (!v.has_alts)
      continue;

    // We store the contig name (we assume VCF to be ordered, i.e.,
    // contigs are not mixed up)
    if(used_seq_names.empty() || v.seq_name.compare(used_seq_names.back()) != 0)
      last_seq_name = v.seq_name;

    VB *vb = &(vbs.back());
    if(!vb->empty() && (!vb->back().is_rknear_to(v, opt::k) || v.seq_name.compare(last_seq_name) != 0)) {
      if(vbs.size() == 5000) { // we can process the blocks
	build_signatures(vbs, refs);
	set_coverages(vbs, bf, ref_bf);
	compute_genotypes(vbs);
	print_variants(vbs);
	vbs.clear();
      }

      vbs.push_back(VB(opt::k, opt::error_rate));
      if(v.seq_name.compare(last_seq_name) != 0)
	last_seq_name = v.seq_name;
    }
    // Here we cannot use pointer vb since we could have add a new vb to the vector
    vbs.back().add_variant(v);
  }
  if(!vbs.empty()) {
    build_signatures(vbs, refs);
    set_coverages(vbs, bf, ref_bf);
    compute_genotypes(vbs);
    print_variants(vbs);
    vbs.clear();
  }


  bcf_hdr_destroy(vcf_header);
  bcf_destroy(vcf_record);
  bcf_close(vcf);

  cout.flush();

  pelapsed("Execution completed");

  return 0;
}
