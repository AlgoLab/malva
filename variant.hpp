#ifndef _VARIANT_HPP_
#define _VARIANT_HPP_

struct Variant {
  std::string seq_name;
  int ref_pos;                               // Variant position 0-based
  std::string ref_sub;                       // Reference base{s}
  int ref_size;                              // Len of reference base{s}
  std::vector<std::string> alts;             // List of alternatives
  int min_size;                              // Length of the shortest string (ref and alts)
  int max_size;                              // Length of the longest string (ref and alts)
  bool is_good = true;                       // false if no alternatives, i.e. only <CN>
  std::vector<std::pair<int,int>> genotypes; // full list of genotypes
  std::vector<bool> phasing;                 // true if genotype i-th is phased, false otherwise
  std::vector<int> positive_samples;         // indices of samples for which genotype is different from 00

  Variant() {}

  Variant(bcf_hdr_t *vcf_header, bcf1_t *vcf_record) {
    seq_name = bcf_hdr_id2name(vcf_header, vcf_record->rid);
    ref_pos = vcf_record->pos;
    ref_sub = vcf_record->d.allele[0];
    ref_size = ref_sub.size();
    for(int i = 1; i<vcf_record->n_allele; ++i) {
      char* curr_alt = vcf_record->d.allele[i];
      if(curr_alt[0] != '<')
        alts.push_back(std::string(curr_alt));
    }
    set_attributes();
    if(is_good)
      extract_genotypes(vcf_header, vcf_record);
  }

  /**
   * Set the is_good flag and the min/max size of the variant
   **/
  void set_attributes() {
    if(alts.size() == 0)
      is_good = false;
    else {
      min_size = (int)ref_sub.size();
      max_size = (int)ref_sub.size();
      for (uint i = 0; i < alts.size(); i++) {
        if ((int)alts[i].size() < min_size)
          min_size = alts[i].size();
        if ((int)alts[i].size() > max_size)
          max_size = alts[i].size();
      }
    }
  }
  
  void extract_genotypes(bcf_hdr_t *vcf_header, bcf1_t *vcf_record) {
    // number of individuals from header
    int n_individuals = bcf_hdr_nsamples(vcf_header);
    int32_t *gt_arr = NULL, ngt = 0;
    int ngt_ret_value = bcf_get_genotypes(vcf_header, vcf_record, &gt_arr, &ngt);
    /***
     * If the record contains GT fields, ngt == ngt_ret_value == #GT.
     * Otherwise, ngt_ret_value is <= 0 that means an error occurred.
     ***/
    if(ngt_ret_value <= 0) {
      // std::cout << "The record doesn't contain GT information" << std::endl;
      is_good = false;
      return;
    }

    int ploidy = ngt/n_individuals;
    for (int i=0; i<n_individuals; ++i) {
      int32_t *curr_gt = gt_arr + i*ploidy;
      // Here, I'm assuming ploidy = 2. Otherwise, loop til ploidy
      genotypes.push_back(std::make_pair(bcf_gt_allele(curr_gt[0]),
                                         bcf_gt_allele(curr_gt[1])));
      if(bcf_gt_allele(curr_gt[0]) > 0 or bcf_gt_allele(curr_gt[1]) > 0)
        positive_samples.push_back(i);
      if(bcf_gt_is_phased(curr_gt[1]))
        // this works, but I'm not 100% sure it's sufficient
        phasing.push_back(true);
      else
        phasing.push_back(false);
    }
    free(gt_arr);
  }

  /**
   * Return the i-th allele, 0 is the reference
   **/
  std::string get_allele(const int &i) const {
    if(i==0)
      return ref_sub;
    else
      return alts[i-1];
  }

  /**
   * Given an allele, return its position in the list (1-based since 0 is the reference)
   **/
  int get_allele_index(const std::string &a) const {
    int i = 1;
    for(const std::string &all : alts) {
      if(all.compare(a) == 0)
        return i;
      ++i;
    }
    return -1;
  }
  
  /**
   * Return true if (this) variant overlaps with var
   * i.e. they are incompatible
   **/
  bool overlaps_with(const Variant &var) const {
    return (ref_pos <= var.ref_pos) &&
      (var.ref_pos < ref_pos + ref_size);
  }

  /**
   * Return true if (this) variant is k/2-near (on the right) to var
   **/
  bool is_near_to(const Variant &var, const int &k, const int &sum_to_add) const {
    return ref_pos + ref_size - min_size - 1 + sum_to_add + ceil((float)k / 2) >= var.ref_pos;
  }

  /**
   * Return true if (this) variant is k/2-near (on the right) to var
   **/
  bool is_near_to(const Variant &var, const int &k) const {
    return is_near_to(var, k, 0);
  }
  
  /**
   * Return true if (this) variant is gt-compatible with var
   **/
  bool is_compatible_with(const Variant &var) {
    for(const uint i : positive_samples) {
      if(phasing[i]) {
        if((genotypes[i].first != 0 && var.genotypes[i].first != 0) ||
           (genotypes[i].second != 0 && var.genotypes[i].second != 0))
          return true;
      } else {
        if((genotypes[i].first != 0 || genotypes[i].second != 0) &&
           (var.genotypes[i].first != 0 || var.genotypes[i].second != 0))
          return true;
      }
    }
    return false;
  }
};

#endif
