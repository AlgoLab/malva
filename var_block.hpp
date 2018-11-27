#ifndef _VAR_BLOCK_HPP_
#define _VAR_BLOCK_HPP_

#include "variant.hpp"

typedef std::map<int, std::map<int, std::vector<std::vector<std::string>>>> VK_GROUP;

/**
 * Extend a container with another
 **/
template <typename T> void extend(T &V1, const T &V2) {
  V1.reserve(V1.size() + distance(V2.begin(), V2.end()));
  V1.insert(V1.end(), V2.begin(), V2.end());
}

class VB {
private: // attributes
  std::vector<Variant> variants;
  const char *reference;
  int k;
  int number_variants_out = 0;

public:
  VB(const char *_reference, const int &_k) {
    reference = _reference;
    k = _k;
  }
  ~VB() {}

  bool is_near_to_last(const Variant &v) {
    return are_near(&variants.back(), &v, k);
  }

  void add_variant(const Variant &v) {
    variants.push_back(v);
  }

  Variant get_variant(const int &i) const {
    return variants[i];
  }

  bool empty() {
    return variants.empty();
  }

  void clear() {
    variants.clear();
  }

  VK_GROUP extract_kmers() {
    VK_GROUP kmers;

    for(uint v_index = 0; v_index < variants.size(); ++v_index) {
      std::map<int, std::vector<std::vector<std::string>>> _kmers;

      Variant *v = &variants[v_index];

      std::vector<std::vector<int>> right_combs = get_combs_on_the_right(v_index);
      std::vector<std::vector<int>> left_combs = get_combs_on_the_left(v_index);

      std::vector<std::vector<int>> combs = combine_combs(left_combs, right_combs, v_index);

      for(const std::vector<int> &comb : combs) {
        std::vector<std::string> ref_subs = get_ref_subs(comb);
        std::set<std::vector<std::string>> alt_allele_combs = build_alternative_alleles_combs(comb, v_index);

        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // !!! the body of this for could be split in more methods !!!
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        for(const std::vector<std::string> aac : alt_allele_combs) {
	  std::vector<std::string> ksss; //kmers sequences
	  std::string mid_allele;

	  if(aac.size() == 1 && aac[0].size() >= (uint)k) {
	    mid_allele = aac[0];

	    std::string kmer = mid_allele.substr(0,k);
	    transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
	    ksss.push_back(kmer);

	    for(uint p = k; p < mid_allele.size(); ++p) {
	      char c = toupper(mid_allele[p]);
	      kmer.erase(0,1);
	      kmer += c;
	      ksss.push_back(kmer);
	    }
	  } else {
	    std::string kmer = "";
	    int mid_pos_in_kmer = 0;
	    for(uint j = 0; j<aac.size(); ++j) {
	      std::string rs;
	      if(j>=ref_subs.size())
		rs = "";
	      else
		rs = ref_subs[j];

	      // store the position of the mid allele and the mid allele inside the kmer to use later for extending the kmer
	      if(comb[j] == (int)v_index) {
		mid_pos_in_kmer = kmer.size();
		mid_allele = aac[j];
	      }
	      kmer += aac[j] + rs;
	    }

	    // get how much we must extend or cut
	    int first_part_size = mid_pos_in_kmer + mid_allele.size()/2;
	    int second_part_size = kmer.size() - first_part_size;
	    int missing_prefix = k/2 - first_part_size;
	    int missing_suffix = ceil((float)k/2) - second_part_size;

	    // extending/cutting on the left
	    if(missing_prefix >= 0) {
	      Variant *first_var_in_comb = &variants[comb.front()];
	      std::string prefix (reference + first_var_in_comb->ref_pos - missing_prefix, missing_prefix);
	      kmer = prefix + kmer;
	    } else
	      kmer.erase(0,abs(missing_prefix));

	    // extending/cutting on the right
	    if(missing_suffix >= 0) {
	      Variant *last_var_in_comb = &variants[comb.back()];
	      std::string suffix (reference + last_var_in_comb->ref_pos + last_var_in_comb->ref_size, missing_suffix);
	      kmer += suffix;
	    } else
	      kmer.erase(kmer.size() - abs(missing_suffix), abs(missing_suffix));
	    transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);

	    ksss.push_back(kmer);
	  }

          // add ksss (to _kmers)
          int allele_index = v->get_allele_index(mid_allele);
          if(_kmers.find(allele_index) != _kmers.end()) {
            _kmers[allele_index].push_back(ksss);
          } else {
            std::vector<std::vector<std::string>> tmp_kmers;
            tmp_kmers.push_back(ksss);
            _kmers[allele_index] = tmp_kmers;
          }
        }
        kmers[v_index] = _kmers;
      }
    }
    return kmers;
  }

  void output_variants(std::vector<GT> genotypes) {
    for(uint i=0; i<variants.size(); ++i) {
      Variant *v = &variants[i];
      std::cout << v->seq_name << '\t' << v->ref_pos + 1 << '\t'
                << v->idx << '\t' << v->ref_sub << '\t';
      uint varc = 0;
      for(const std::string &alt : v->alts) {
        std::cout << alt;
        ++varc;
        if(varc != v->alts.size())
          std::cout << ',';
      }
      std::cout << "\t" << v->quality
                << "\t" << v->filter
                << "\t" << v->info
                << "\tGT:GQ\t" << genotypes[i].first << ":" << (int)round(genotypes[i].second*100) << "\n";
    }
  }

private: // methods
  const VB &operator=(const VB &other) { return *this; }
  const VB &operator=(const VB &&other) { return *this; }

  //- General methods for managing variants --------------------------
  //------------------------------------------------------------------
  /**
   * Return true if variant v1 overlaps with variant v2
   * i.e. they are incompatible
   **/
  bool are_overlapping(const Variant *v1, const Variant *v2) const {
    return (v1->ref_pos <= v2->ref_pos) &&
      (v2->ref_pos < v1->ref_pos + v1->ref_size);
  }

  /**
   * Return true if variant v1 is k/2-near (on the right) to variant v2
   **/
  bool are_near(const Variant *v1, const Variant *v2, const int &k, const int &sum_to_add = 0) const {
    return v1->ref_pos + v1->ref_size - v1->min_size - 1 + sum_to_add + ceil((float)k / 2) >= v2->ref_pos;
  }

  /**
   * Return true if variant v1 is gt-compatible with variant v2
   **/
  bool are_gt_compatible(const Variant *v1, const Variant *v2) {
    for(const uint i : v1->positive_samples) {
      if(v1->phasing[i]) {
        if((v1->genotypes[i].first != 0 && v2->genotypes[i].first != 0) ||
           (v1->genotypes[i].second != 0 && v2->genotypes[i].second != 0))
          return true;
      } else {
        if((v1->genotypes[i].first != 0 || v1->genotypes[i].second != 0) &&
           (v2->genotypes[i].first != 0 || v2->genotypes[i].second != 0))
          return true;
      }
    }
    return false;
  }
  //------------------------------------------------------------------

  //- Methods for kmers generation  ----------------------------------
  //------------------------------------------------------------------
  /**
   * Given a mid variant, builds all the possible combinations on the right.
   * In building them, we check for:
   *   - k/2-proximity
   *   - overlapping variants
   *   - gt-compatibility
   **/
  std::vector<std::vector<int>> get_combs_on_the_right(const int &i) {
    Variant *mid_v = &variants[i];
    std::vector<std::vector<int>> right_combs; //possible combinations
    std::vector<int> right_sums; //total sum of the combinations (for checking (k/2)-proximity)
    bool halt_flag = false;
    for(uint j=i+1; j<variants.size() && !halt_flag; ++j) {
      Variant *curr_v = &variants[j];

      if(are_overlapping(mid_v, curr_v))
        continue;

      if(right_combs.empty()) { //first var to be added
        if(are_near(mid_v, curr_v, k) and are_gt_compatible(mid_v, curr_v)) {
          std::vector<int> new_comb (1, (int)j);
          right_combs.push_back(new_comb);
          right_sums.push_back(curr_v->ref_size - curr_v->min_size);
        }
      } else {
        if(are_gt_compatible(mid_v, curr_v)) {
          // add the var to all the compatible combinations
          bool added_flag = false;
          for(uint c = 0; c<right_combs.size(); ++c) {
            std::vector<int> comb = right_combs[c];
            int sum = right_sums[c];

            Variant *last_v_in_comb = &variants[comb.back()];

            if(!are_overlapping(last_v_in_comb, curr_v)) {
              added_flag = true;
              if(are_near(mid_v, curr_v, k, sum)) {
                right_combs[c].push_back(j);
                right_sums[c] += curr_v->ref_size - curr_v->min_size;
              }
            }
          }
          if(!added_flag) {
            //if the var has not been added to any combination
            std::vector<std::vector<int>> new_right_combs;
            std::vector<int> new_right_sums;
            for(uint c = 0; c<right_combs.size(); ++c) {
              // shorten the combinations and try to add the var
              std::vector<int> new_comb = right_combs[c];
              int new_sum = right_sums[c];

              Variant *last_v_in_comb = &variants[new_comb.back()];

              while(are_overlapping(last_v_in_comb, curr_v) and !new_comb.empty()) {
                new_comb.pop_back();
                new_sum -= last_v_in_comb->ref_size - last_v_in_comb->min_size;
                last_v_in_comb = &variants[new_comb.back()];
              }
              new_comb.push_back(j);
              if(are_near(mid_v, curr_v, k, new_sum)) {
                added_flag = true;
                new_right_combs.push_back(new_comb);
                new_right_sums.push_back(new_sum + curr_v->ref_size - curr_v->min_size);
              }
            }
            extend(right_combs, new_right_combs);
            extend(right_sums, new_right_sums);

            // if the var has not been added to any combination (neither to shortened ones),
            // we can halt the loop: no successive variants will be added (too far away)
            if(!added_flag)
              halt_flag = true;
          }
        }
      }
    }

    return right_combs;
  }

  /**
   * Given a mid variant, builds all the possible combinations on the left.
   * In building them, we check for:
   *   - k/2-proximity
   *   - overlapping variants
   *   - gt-compatibility
   **/
  std::vector<std::vector<int>> get_combs_on_the_left(const int &i) {
    Variant *mid_v = &variants[i];
    std::vector<std::vector<int>> left_combs; //possible combinations
    std::vector<int> left_sums; //total sum of the combinations (for checking (k/2)-proximity)

    bool halt_flag = false;
    for(int j=i-1; j>=0 && !halt_flag; --j) {
      Variant *curr_v = &variants[j];

      if(are_overlapping(curr_v, mid_v))
        continue;

      if(left_combs.empty()) { //first var to be added
        if(are_near(curr_v, mid_v, k) and are_gt_compatible(mid_v, curr_v)) {
          std::vector<int> new_comb (1, (int)j);
          left_combs.push_back(new_comb);
          left_sums.push_back(curr_v->ref_size - curr_v->min_size);
        }
      } else {
        if(are_gt_compatible(mid_v, curr_v)) {
          // add the var to all the compatible combinations
          bool added_flag = false;
          for(uint c = 0; c<left_combs.size(); ++c) {
            std::vector<int> comb = left_combs[c];
            int sum = left_sums[c];

            Variant *last_v_in_comb = &variants[comb.back()];

            if(!are_overlapping(curr_v, last_v_in_comb)) {
              added_flag = true;
              if(are_near(curr_v, mid_v, k, sum)) {
                left_combs[c].push_back(j);
                left_sums[c] += curr_v->ref_size - curr_v->min_size;
              }
            }
          }
          if(!added_flag) {
            //if the var has not been added to any combination
            std::vector<std::vector<int>> new_left_combs;
            std::vector<int> new_left_sums;
            for(uint c = 0; c<left_combs.size(); ++c) {
              // shorten the combinations and try to add the var
              std::vector<int> new_comb = left_combs[c];
              int new_sum = left_sums[c];

              Variant *last_v_in_comb = &variants[new_comb.back()];

              while(are_overlapping(curr_v, last_v_in_comb) and !new_comb.empty()) {
                new_comb.pop_back();
                new_sum -= last_v_in_comb->ref_size - last_v_in_comb->min_size;
                last_v_in_comb = &variants[new_comb.back()];
              }
              new_comb.push_back(j);
              if(are_near(curr_v, mid_v, k, new_sum)) {
                added_flag = true;
                new_left_combs.push_back(new_comb);
                new_left_sums.push_back(new_sum + curr_v->ref_size - curr_v->min_size);
              }
            }
            extend(left_combs, new_left_combs);
            extend(left_sums, new_left_sums);

            // if the var has not been added to any combination (neither to shortened ones),
            // we can halt the loop: no successive variants will be added (too far away)
            if(!added_flag)
              halt_flag = true;
          }
        }
      }
    }

    return left_combs;
  }

  /**
   * Combine all left_combs with all right_combs, and placing i-th variant between them
   **/
  std::vector<std::vector<int>> combine_combs(std::vector<std::vector<int>> &left_combs,
                                              const std::vector<std::vector<int>> &right_combs,
                                              const int &i) {
    std::vector<std::vector<int>> full_combs;
    if(left_combs.empty() && right_combs.empty()) {
      std::vector<int> comb;
      comb.push_back(i);
      full_combs.push_back(comb);
    } else if(left_combs.empty()) {
      std::vector<int> lcomb;
      lcomb.push_back(i);
      std::vector<int> comb = lcomb;
      for(const std::vector<int> &rcomb : right_combs) {
        extend(comb, rcomb);
        full_combs.push_back(comb);
        comb = lcomb;
      }
    } else {
      for(std::vector<int> &lcomb : left_combs) {
        std::reverse(lcomb.begin(), lcomb.end());
        lcomb.push_back(i);
        std::vector<int> comb = lcomb;
        if(right_combs.empty()) {
          full_combs.push_back(comb);
        }
        for(const std::vector<int> &rcomb : right_combs) {
          extend(comb, rcomb);
          full_combs.push_back(comb);
          comb = lcomb;
        }
      }
    }
    return full_combs;
  }

  /**
   * Return the reference substring between considered variants
   **/
  std::vector<std::string> get_ref_subs(const std::vector<int> &comb) {
    std::vector<std::string> ref_subs;

    int last_end = -1;
    for(const int &index : comb) {
      Variant *v = &variants[index];
      if(last_end == -1) {
        last_end = v->ref_pos + v->ref_size;
        continue;
      }
      std::string ref_sub (reference + last_end, v->ref_pos - last_end);
      ref_subs.push_back(ref_sub);
      last_end = v->ref_pos + v->ref_size;
    }
    return ref_subs;
  }

  /**
   * Builds and returns all the possible combination of alternative alleles,
   * with respect to GTs.
   * !!! For now, I'm assuming phased GT !!!
   **/
  std::set<std::vector<std::string>> build_alternative_alleles_combs(const std::vector<int> &comb,
                                                                     const int &central_index) {
    // A set to avoid duplicate elements
    std::set<std::vector<std::string>> aacs;
    Variant *central_v = &variants[central_index];
    // For each individual having this variant
    for(const int &gt_i : central_v->positive_samples) {
      std::vector<std::string> aac;
      // If GT is not 0|0, then build the vector of alleles
      if(central_v->genotypes[gt_i].first != 0) {
        for(const int &j : comb)
          aac.push_back(variants[j].get_allele(variants[j].genotypes[gt_i].first));
        aacs.insert(aac);
      }
      aac.clear();
      if(central_v->genotypes[gt_i].second != 0) {
        for(const int &j : comb)
          aac.push_back(variants[j].get_allele(variants[j].genotypes[gt_i].second));
        aacs.insert(aac);
      }
    }
    return aacs;
  }
};

#endif
