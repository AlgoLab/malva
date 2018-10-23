#ifndef _VAR_BLOCK_HPP_
#define _VAR_BLOCK_HPP_

#include "variant.hpp"

typedef std::map<int, std::map<int, std::vector<std::string>>> VK_GROUP;

/**
 * Extend a container with another
 **/
template <typename T> void extend(T &V1, const T &V2) {
  V1.reserve(V1.size() + distance(V2.begin(), V2.end()));
  V1.insert(V1.end(), V2.begin(), V2.end());
}

// /**
//  * Print the content of a container of basic elements.
//  **/
// template <typename T> void print(std::vector<T> A, std::string s) {
//   std::cout << "-- " << s << " ";
//   for(const auto &a : A) {
//     std::cout << a << " ";
//   }
//   std::cout << std::endl;
// }

class VB {
private: // attributes
  std::vector<Variant> variants;
  const char *reference;
  int k;

public:
  VB(const char *_reference, int _k) {
    reference = _reference;
    k = _k;
  }
  ~VB() {}

  bool is_near_to_last(const Variant &v) {
    return variants.back().is_near_to(v, k);
  }
  
  void add_variant(const Variant &v) {
    variants.push_back(v);
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
      std::map<int, std::vector<std::string>> _kmers;

      Variant v = variants[v_index];

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
          std::string kmer = "";
          int mid_pos_in_kmer = 0;
          std::string mid_allele;
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
            Variant first_var_in_comb = variants[comb.front()];
            std::string prefix (reference + first_var_in_comb.ref_pos - missing_prefix, missing_prefix);
            kmer = prefix + kmer;
          } else
            kmer.erase(0,abs(missing_prefix));

          // extending/cutting on the right
          if(missing_suffix >= 0) {
            Variant last_var_in_comb = variants[comb.back()];
            std::string suffix (reference + last_var_in_comb.ref_pos + last_var_in_comb.ref_size, missing_suffix);
            kmer += suffix;
          } else
            kmer.erase(kmer.size() - abs(missing_suffix), abs(missing_suffix));

          // std::cout << "- " << kmer << " (" << kmer.size() << ")" << std::endl;
          //std::cout << missing_prefix << " - (" << first_part_size << "|" << second_part_size << ") - " << missing_suffix << std::endl;

          // add k-mer (to _kmers)
          int allele_index = v.get_allele_index(mid_allele);
          if(_kmers.find(allele_index) != _kmers.end()) {
            _kmers[allele_index].push_back(kmer);
          } else {
            std::vector<std::string> tmp_kmers;
            tmp_kmers.push_back(kmer);
            _kmers[allele_index] = tmp_kmers;
          }
        }
        kmers[v_index] = _kmers;
      }
    }
    return kmers;
  }
  
private: // methods
  const VB &operator=(const VB &other) { return *this; }
  const VB &operator=(const VB &&other) { return *this; }
  /**
   * Given a mid variant, builds all the possible combinations on the right.
   * In building them, we check for:
   *   - k/2-proximity
   *   - overlapping variants
   *   - gt-compatibility
   **/
  std::vector<std::vector<int>> get_combs_on_the_right(const int &i) {
    Variant mid_v = variants[i];
    std::vector<std::vector<int>> right_combs; //possible combinations
    std::vector<int> right_sums; //total sum of the combinations (for checking (k/2)-proximity)
    bool halt_flag = false;
    for(uint j=i+1; j<variants.size() && !halt_flag; ++j) {
      Variant curr_v = variants[j];

      if(mid_v.overlaps_with(curr_v))
        continue;

      if(right_combs.empty()) { //first var to be added
        if(mid_v.is_near_to(curr_v, k) and mid_v.is_compatible_with(curr_v)) {
          std::vector<int> new_comb (1, (int)j);
          right_combs.push_back(new_comb);
          right_sums.push_back(curr_v.ref_size - curr_v.min_size);
        }
      } else {
        if(mid_v.is_compatible_with(curr_v)) {
          // add the var to all the compatible combinations
          bool added_flag = false;
          for(uint c = 0; c<right_combs.size(); ++c) {
            std::vector<int> comb = right_combs[c];
            int sum = right_sums[c];

            Variant last_v_in_comb = variants[comb.back()];
            
            if(!last_v_in_comb.overlaps_with(curr_v)) {
              added_flag = true;
              if(mid_v.is_near_to(curr_v, k, sum)) {
                right_combs[c].push_back(j);
                right_sums[c] += curr_v.ref_size - curr_v.min_size;
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

              Variant last_v_in_comb = variants[new_comb.back()];
            
              while(last_v_in_comb.overlaps_with(curr_v) and !new_comb.empty()) {
                new_comb.pop_back();
                new_sum -= last_v_in_comb.ref_size - last_v_in_comb.min_size;
                last_v_in_comb = variants[new_comb.back()];
              }
              new_comb.push_back(j);
              if(mid_v.is_near_to(curr_v, k, new_sum)) {
                added_flag = true;
                new_right_combs.push_back(new_comb);
                new_right_sums.push_back(new_sum + curr_v.ref_size - curr_v.min_size);
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
    Variant mid_v = variants[i];
    std::vector<std::vector<int>> left_combs; //possible combinations
    std::vector<int> left_sums; //total sum of the combinations (for checking (k/2)-proximity)

    bool halt_flag = false;
    for(int j=i-1; j>=0 && !halt_flag; --j) {
      Variant curr_v = variants[j];

      if(curr_v.overlaps_with(mid_v))
        continue;

      if(left_combs.empty()) { //first var to be added
        if(curr_v.is_near_to(mid_v, k) and mid_v.is_compatible_with(curr_v)) {
          std::vector<int> new_comb (1, (int)j);
          left_combs.push_back(new_comb);
          left_sums.push_back(curr_v.ref_size - curr_v.min_size);
        }
      } else {
        if(mid_v.is_compatible_with(curr_v)) {
          // add the var to all the compatible combinations
          bool added_flag = false;
          for(uint c = 0; c<left_combs.size(); ++c) {
            std::vector<int> comb = left_combs[c];
            int sum = left_sums[c];

            Variant last_v_in_comb = variants[comb.back()];
            
            if(!curr_v.overlaps_with(last_v_in_comb)) {
              added_flag = true;
              if(curr_v.is_near_to(mid_v, k, sum)) {
                left_combs[c].push_back(j);
                left_sums[c] += curr_v.ref_size - curr_v.min_size;
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

              Variant last_v_in_comb = variants[new_comb.back()];

              while(curr_v.overlaps_with(last_v_in_comb) and !new_comb.empty()) {
                new_comb.pop_back();
                new_sum -= last_v_in_comb.ref_size - last_v_in_comb.min_size;
                last_v_in_comb = variants[new_comb.back()];
              }
              new_comb.push_back(j);
              if(curr_v.is_near_to(mid_v, k, new_sum)) {
                added_flag = true;
                new_left_combs.push_back(new_comb);
                new_left_sums.push_back(new_sum + curr_v.ref_size - curr_v.min_size);
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
      Variant v = variants[index];
      if(last_end == -1) {
        last_end =v.ref_pos + v.ref_size;
        continue;
      }
      std::string ref_sub (reference + last_end, v.ref_pos - last_end);
      ref_subs.push_back(ref_sub);
      last_end = v.ref_pos + v.ref_size;
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
    Variant central_v = variants[central_index];
    // For each individuals
    for(uint gt_i=0; gt_i<central_v.genotypes.size(); ++gt_i) {
      std::vector<std::string> aac;
      // If GT is not 0|0, then build the vector of alleles
      if(central_v.genotypes[gt_i].first != 0) {
        for(const int &j : comb)
          aac.push_back(variants[j].get_allele(variants[j].genotypes[gt_i].first));
        aacs.insert(aac);
      }
      aac.clear();
      if(central_v.genotypes[gt_i].second != 0) {
        for(const int &j : comb)
          aac.push_back(variants[j].get_allele(variants[j].genotypes[gt_i].second));
        aacs.insert(aac);
      }
    }
    return aacs;
  }
};

#endif
