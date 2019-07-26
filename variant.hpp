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

#include "kstring.h"

typedef vector<string> signature;

struct GT { // This name is already used in variant.hpp
  uint8_t a1, a2;
  bool phased;
  long double quality;

  GT() {
    a1 = 0;
    a2 = 0;
    phased = true;
    quality = 0;
  }

  GT(uint8_t _a1, uint8_t _a2, bool _phased, long double _quality = 0) {
    a1 = _a1;
    a2 = _a2;
    phased = _phased; //|| a1 == a2;
    quality = _quality;
  }

  // we could overload operator<<
  string to_str() {
    return phased ? to_string(a1)+"|"+to_string(a2)+":"+to_string((int)round(quality*100)) : to_string(a1)+"/"+to_string(a2)+":"+to_string((int)round(quality*100));
  }

  void print() {
    if(phased)
      cout << (int)a1 << "|" << (int)a2 << endl;
    else
      cout << (int)a1 << "/" << (int)a2 << endl;
  }
};

struct Variant {
  kstring_t vcf_line;                   // kstring_t from the vcf
  string seq_name;                      // Chromosome/Contig
  int ref_pos;                          // Variant position 0-based
  string idx;                           // ID
  string ref_sub;                       // Reference base{s}
  vector<string> alts;                  // List of alternatives
  float quality;                        // Quality field
  string filter;                        // Filter field
  string info;                          // Info field
  vector<GT> genotypes;                 // Genotypes/Samples
  int n_individuals;                    // Number of samples
  int ref_size;                         // Len of reference base{s}
  int min_size;                         // Length of the shortest string (ref and alts)
  int max_size;                         // Length of the longest string (ref and alts)
  bool has_alts = true;                 // false if no alternatives, i.e. only <CN>
  bool is_present = true;               // false if no sample has this variant
  vector<float> frequencies;            // Allele frequency in the considered population
  vector<float> coverages;              // Allele coverages (computed from input sample)
  vector<vector<signature>> SIGNS;      // Alleles' signatures
  GT computed_gt;                       // Computed genotypes (best one)

  Variant() {}

  ~Variant() {}

  Variant(htsFile *vcf, bcf_hdr_t *vcf_header, bcf1_t *vcf_record, const string &freq_key) {
    // vcf_line = new kstring_t();
    ks_initialize(&vcf_line);
    kputsn(vcf->line.s, vcf->line.l, ks_clear(&vcf_line));

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
    extract_alternative_alleles(vcf_record);
    coverages.resize(alts.size() + 1, 0); //+1 for the reference allele
    SIGNS.resize(alts.size() + 1, vector<signature>()); //+1 for the reference allele
    quality = vcf_record->qual;
    filter = "PASS"; // TODO: get filter string from VCF
    info = ".";      // TODO: get info string from VCF
    // Set sizes and has_alts flag
    set_sizes();
    n_individuals = bcf_hdr_nsamples(vcf_header);
    if (has_alts) // Populate frequencies vector
      extract_frequencies(vcf_header, vcf_record, freq_key);
  }

  // If we want, we can use pointers in main to avoid this function
  // and call ks_free directly from the destructor (I think). Right
  // now, we cannot call it from there since the kstring would be
  // freed more than once ("double free or corruption"). Now we have
  // to call this function manually when we know that we would not use
  // the variant anymore.
  void free() {
    ks_free(&vcf_line);
  }
  
  /**
   * Parse the alternative alleles and fill the corresponding vector.
   **/
  void extract_alternative_alleles(bcf1_t *vcf_record) {
    for (int i = 1; i < vcf_record->n_allele; ++i) {
      char *curr_alt = vcf_record->d.allele[i];
      if (curr_alt[0] != '<')
        alts.push_back(string(curr_alt));
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
                           const string &freq_key) {
    int ndst = 0;
    float *altall_freqs = NULL;
    bcf_get_info_float(vcf_header, vcf_record, freq_key.c_str(),
                       &altall_freqs, &ndst);

    frequencies.push_back(0); // First element is reserved for reference allele
    for (uint i = 0; i < alts.size(); ++i) {
      float *freq = altall_freqs + i;
      frequencies.push_back(freq[0]);
    }
    // Here we compute the frequency of the reference allele
    frequencies[0] = 1.0 - accumulate(frequencies.begin(), frequencies.end(), 0.0);

    if (frequencies[0] == 1.0)
      is_present = false;
  }

  /**
   * Given the kstring_t, cut the prefix and returns the suffix
   * containing only the samples. [used by fill_genotypes]
   **/
  char* get_samples(char *line) {
    int offset = 0;
    char* fmt_suffix;
    int fmt_size = 0;
    bool fmt_flag = false;
    while(!fmt_flag) {
      offset += strlen(line + offset) + 1;
      /**
       * From VCF format specification: "The first sub-field must always
       * be the genotype (GT) if it is present. There are no required
       * sub-fields."
       **/
      if(strncmp(line + offset, "GT\t", 3) == 0 || strncmp(line + offset, "GT:", 3) == 0) {
	fmt_flag = true;
	fmt_suffix = line + offset;
	char* fmt_end = strchr(fmt_suffix, '\t');
	if(fmt_end==NULL) {
	  cerr << "Error - No samples" << endl;
	  exit(1);
	}
	fmt_size = fmt_end - fmt_suffix;
      }
    }
    return line + offset + fmt_size + 1;
  }

  /**
   * Transform a char* representing a genotype to the struct GT. [used
   * by fill_genotypes]
   **/
  GT extract_genotype(char* gt) {
    if(!strcmp(gt, ".") or !strcmp(gt, "./.") or !strcmp(gt, ".|.")) {
      return GT();
    }
    bool phased = strchr(gt, '|') != NULL;
    int all1, all2;
    char *gt_ptr;
    char *all = strtok_r(gt,"/|", &gt_ptr);
    all1 = atoi(all);
    all = strtok_r(NULL,"/|", &gt_ptr);
    all2 = all==NULL ? all1 : atoi(all);

    assert(all1 <= alts.size() && all2 <= alts.size());

    return GT(all1, all2, phased);
  }

  /**
   * Builds the genotypes from the kstring_t and fills the
   * corresponding vector
   **/
  void fill_genotypes() {
    genotypes.resize(n_individuals);

    char *samples = get_samples(vcf_line.s);
    int i = 0;
    for(;;) {
      char *next_tab = strchr(samples, '\t');
      char *gt_ptr;
      char *gtc = strtok_r(samples,":\t", &gt_ptr);
      GT gt = extract_genotype(gtc);
      genotypes[i] = gt;
      if(next_tab == NULL) break;
      samples = next_tab+1;
    }
  }

  /**
   * Return true if variant v is k/2-near (on the right) to this
   * variant (|--- this -- v --->). The "sum_to_add" variable must be
   * used when "concatenating" more variants (reference may expand or
   * shrink due to indels).
   **/
  bool is_rknear_to(const Variant &v, const int &k, const int &sum_to_add = 0) const {
    return ref_pos + ref_size - min_size - 1 + sum_to_add + ceil((float)k / 2) >= v.ref_pos;
  }

  /**
   * Return true if this variant overlaps with variant v
   * i.e. they are incompatible
   **/
  bool overlaps_with(const Variant &v) const {
    return (ref_pos <= v.ref_pos) && (v.ref_pos < ref_pos + ref_size);
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

  void set_coverage(const int &i, const float &cov) {
    // maybe we can add some control here
    coverages[i] = cov;
  }

  void add_signature(const string &allele, const signature &S) {
    int allele_index = get_allele_index(allele);
    SIGNS[allele_index].push_back(S);
  }

  void set_genotype(const GT &gt) {
    // maybe we can add some control here
    computed_gt = gt;
  }
};

#endif
