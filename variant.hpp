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

#ifndef _VARIANT_HPP_
#define _VARIANT_HPP_

#include <numeric>
#include <string>
#include <vector>

using namespace std;

typedef pair<string, long double> GT;

struct Variant {
  string seq_name;
  int ref_pos;                        // Variant position 0-based
  string idx;                         // ID
  string ref_sub;                     // Reference base{s}
  vector<string> alts;                // List of alternatives
  float quality;                      // Quality field
  string filter;                      // Filter field
  string info;                        // Info field
  vector<pair<int, int>> genotypes;   // full list of genotypes
  vector<bool> phasing;               // true if genotype i-th is phased, false otherwise
  int ref_size;                       // Len of reference base{s}
  int min_size;                       // Length of the shortest string (ref and alts)
  int max_size;                       // Length of the longest string (ref and alts)
  bool has_alts = true;               // false if no alternatives, i.e. only <CN>
  bool is_present = true;             // false if no sample has this variant
  vector<float> frequencies;          // Allele frequency in the considered population
  vector<uint> coverages;             // Allele coverages (computed from input sample)
  vector<GT> computed_gts;            // Computed genotypes

  Variant() {}

  Variant(bcf_hdr_t *vcf_header, bcf1_t *vcf_record, const string &freq_key, const bool uniform) {
    seq_name = bcf_hdr_id2name(vcf_header, vcf_record->rid);
    ref_pos = vcf_record->pos;
    idx = vcf_record->d.id;
    ref_sub = vcf_record->d.allele[0];
    /**
     * There is vcf_record->rlen for the ref_size but sometimes it
     * returns a wrong value. The only example I found is when all
     * alternatives are <CN> but maybe there could be more cases...
     **/
    ref_size = (int)ref_sub.size();
    for (int i = 1; i < vcf_record->n_allele; ++i) {
      char *curr_alt = vcf_record->d.allele[i];
      if (curr_alt[0] != '<')
        alts.push_back(string(curr_alt));
    }
    coverages.resize(alts.size() + 1, 0); // +1 for the reference allele
    quality = vcf_record->qual;
    filter = "PASS"; // TODO: get filter string from VCF
    info = ".";      // TODO: get info string from VCF
    // Set sizes and has_alts flag
    set_sizes();
    if (has_alts) {
      // Populate frequencies vector
      extract_frequencies(vcf_header, vcf_record, freq_key, uniform);
      if (is_present)
        // Populate genotypes and phasing
        extract_genotypes(vcf_header, vcf_record);
    }
  }

  /**
   * Set the has_alts flag and the min/max size of the variant
   **/
  void set_sizes() {
    if (alts.size() == 0)
      has_alts = false;
    else {
      min_size = ref_size;
      max_size = ref_size;
      for (uint i = 0; i < alts.size(); i++) {
        if ((int)alts[i].size() < min_size)
          min_size = alts[i].size();
        if ((int)alts[i].size() > max_size)
          max_size = alts[i].size();
      }
    }
  }

  void extract_frequencies(bcf_hdr_t *vcf_header, bcf1_t *vcf_record,
                           const string &freq_key, const bool uniform) {
    int ndst = 0;
    float *altall_freqs = NULL;
    bcf_get_info_float(vcf_header, vcf_record, freq_key.c_str(),
                       &altall_freqs, &ndst);

    if(!uniform) {
      frequencies.push_back(0); // First element is reserved for reference allele
      for (uint i = 0; i < alts.size(); ++i) {
	float *freq = altall_freqs + i;
	frequencies.push_back(freq[0]);
      }
      // Here we compute the frequency of the reference allele
      frequencies[0] = 1.0 - accumulate(frequencies.begin(), frequencies.end(), 0.0);
      if(frequencies[0] < 0) // this to avoid bad approximation
	frequencies[0] = 0.0;
    } else {
      float uniform_freq = 1.0/(alts.size()+1);
      for (uint i = 0; i<alts.size()+1; ++i)
	frequencies.push_back(uniform_freq);
    }

    if (frequencies[0] == 1.0)
      is_present = false;
  }

  void extract_genotypes(bcf_hdr_t *vcf_header, bcf1_t *vcf_record) {
    // number of individuals from header
    int n_individuals = bcf_hdr_nsamples(vcf_header);
    int32_t *gt_arr = NULL, ngt = 0;
    int ngt_ret_value =
        bcf_get_genotypes(vcf_header, vcf_record, &gt_arr, &ngt);
    /***
     * If the record contains GT fields, ngt == ngt_ret_value == #GT.
     * Otherwise, ngt_ret_value is <= 0 that means an error occurred.
     ***/
    if (ngt_ret_value <= 0) {
      // cout << "The record doesn't contain GT information" << endl;
      has_alts = false;
      return;
    }

    int ploidy = ngt / n_individuals;
    for (int i = 0; i < n_individuals; ++i) {
      int32_t *curr_gt = gt_arr + i * ploidy;
      // Here, I'm assuming ploidy = 2. Otherwise, loop til ploidy
      int all_1;
      int all_2;
      bool is_phased = false;
      if (curr_gt[1] == bcf_int32_vector_end) {
        all_1 = bcf_gt_allele(curr_gt[0]);
        all_2 = bcf_gt_allele(curr_gt[0]);
	is_phased = true;
      } else {
        all_1 = bcf_gt_allele(curr_gt[0]);
        all_2 = bcf_gt_allele(curr_gt[1]);
        if (bcf_gt_is_phased(curr_gt[1]))
          // this works, but I'm not 100% sure it's sufficient
          is_phased = true;
      }
      if(all_1 < 0) all_1 = 0;
      if(all_2 < 0) all_2 = 0;

      // FIXME: we store two alleles per sample even when we have an
      // haploid VCF. We manage haploid VCF in a special way in
      // var_block.hpp
      genotypes.push_back(make_pair(all_1, all_2));
      phasing.push_back(is_phased);
    }
    free(gt_arr);
  }

  /**
   * Return the i-th allele, 0 is the reference
   **/
  string get_allele(const int &i) const {
    if (i == 0)
      return ref_sub;
    else
      return alts[i - 1];
  }

  /**
   * Given an allele, return its position in the list (1-based since 0 is the
   *reference)
   **/
  int get_allele_index(const string &a) const {
    if(ref_sub.compare(a) == 0)
      return 0;
    int i = 1;
    for (const string &all : alts) {
      if (all.compare(a) == 0)
        return i;
      ++i;
    }
    return -1;
  }

  void set_coverage(const int i, const uint cov) {
    // maybe we can add some control here
    coverages[i] = cov;
  }

  void add_genotype(const GT &gt) {
    // maybe we can add some control here
    computed_gts.push_back(gt);
  }
};

#endif
