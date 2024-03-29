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

#ifndef _VAR_BLOCK_HPP_
#define _VAR_BLOCK_HPP_

#include "variant.hpp"
#include <map>
#include <unordered_set>
#include <string_view>
#include "xxhash.h"

using namespace std;

typedef map<int, map<int, vector<vector<string>>>> VK_GROUP;

/**
 * Custom hash function for unordered_set in build_alleles_combs
 **/
struct VectorHash
{
  size_t operator()(const vector<string_view> &v) const
  {
    size_t hash = 0;
    for (string_view i : v)
    {
      hash ^= XXH3_64bits(i.data(), i.size()) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
    }
    return hash;
  }
};

/**
 * Extend a container with another
 **/
template <typename T>
void extend(T &V1, const T &V2)
{
  V1.reserve(V1.size() + distance(V2.begin(), V2.end()));
  V1.insert(V1.end(), V2.begin(), V2.end());
}

class VB
{
private: // attributes
  vector<Variant> variants;
  int k;
  float error_rate;
  int number_variants_out = 0;

public:
  VB(const int &_k, const float _error_rate)
  {
    k = _k;
    error_rate = _error_rate;
  }
  ~VB() {}

  bool is_near_to_last(const Variant &v)
  {
    return are_near(variants.back(), v, k);
  }

  void add_variant(const Variant &v) { variants.push_back(v); }

  void set_variant_coverage(const int &v, const int &i, const float &cov)
  {
    variants[v].set_coverage(i, cov);
  }

  Variant get_variant(const int &i) const { return variants[i]; }

  bool empty() { return variants.empty(); }

  void clear() { variants.clear(); }

  VK_GROUP extract_kmers(const string &reference, const bool haploid)
  {
    VK_GROUP kmers;
    for (uint v_index = 0; v_index < variants.size(); ++v_index)
    {
      map<int, vector<vector<string>>> _kmers;

      Variant *v = &variants[v_index];

      if (!v->is_present || v->ref_pos < k || v->ref_pos > (int)reference.size() - k)
      {
        // Here we could use k/2 but with k the probability to not
        // exceed the reference length due to insertions is
        // higher. Anyway, it shouldn't be common to have variants
        // near to reference ends
        kmers[v_index] = _kmers;
        continue;
      }

      vector<vector<int>> right_combs =
          get_combs_on_the_right(v_index);
      vector<vector<int>> left_combs = get_combs_on_the_left(v_index);

      vector<vector<int>> combs =
          combine_combs(left_combs, right_combs, v_index);

      for (const vector<int> &comb : combs)
      {
        vector<string> ref_subs = get_ref_subs(comb, reference);
        unordered_set<vector<string_view>, VectorHash> alt_allele_combs = build_alleles_combs(comb, v_index, haploid);

        for (const vector<string_view> &aac : alt_allele_combs)
        {
          vector<string> ksss; // kmers sequences
          string_view mid_allele;
          if (aac.size() == 1 && aac[0].size() >= (uint)k)
          {
            mid_allele = aac[0];

            string kmer(mid_allele, 0, k);
            ksss.emplace_back(kmer);

            for (uint p = k; p < mid_allele.size(); ++p)
            {
              char c = mid_allele[p];
              kmer.erase(0, 1);
              kmer += c;
              ksss.emplace_back(kmer);
            }
          }
          else
          {
            string kmer = "";
            int mid_pos_in_kmer = 0;
            for (uint j = 0; j < aac.size(); ++j)
            {
              string rs;
              if (j >= ref_subs.size())
                rs = "";
              else
                rs = ref_subs[j];

              // store the position of the mid allele and the mid allele inside
              // the kmer to use later for extending the kmer
              if (comb[j] == (int)v_index)
              {
                mid_pos_in_kmer = kmer.size();
                mid_allele = aac[j];
              }
              kmer += aac[j];
              kmer += rs;
            }

            // get how much we must extend or cut
            int first_part_size = mid_pos_in_kmer + mid_allele.size() / 2;
            int second_part_size = kmer.size() - first_part_size;
            int missing_prefix = k / 2 - first_part_size;
            int missing_suffix = ceil((float)k / 2) - second_part_size;

            // extending/cutting on the left
            if (missing_prefix >= 0)
            {
              Variant *first_var_in_comb = &variants[comb.front()];
              string prefix(reference,
                            first_var_in_comb->ref_pos - missing_prefix,
                            missing_prefix);
              kmer = prefix + kmer;
            }
            else
              kmer.erase(0, abs(missing_prefix));

            // extending/cutting on the right
            if (missing_suffix >= 0)
            {
              Variant *last_var_in_comb = &variants[comb.back()];
              string suffix(reference,
                            last_var_in_comb->ref_pos + last_var_in_comb->ref_size,
                            missing_suffix);
              kmer += suffix;
            }
            else
              kmer.erase(kmer.size() - abs(missing_suffix),
                         abs(missing_suffix));

            ksss.emplace_back(kmer);
          }

          // add ksss (to _kmers)
          int allele_index = v->get_allele_index(mid_allele);
          if (_kmers.find(allele_index) != _kmers.end())
          {
            _kmers[allele_index].push_back(ksss);
          }
          else
          {
            vector<vector<string>> tmp_kmers;
            tmp_kmers.emplace_back(ksss);
            _kmers[allele_index] = tmp_kmers;
          }
        }
        kmers[v_index] = _kmers;
      }
    }
    return kmers;
  }

  /**
   * Method to compute and store the genotype of each variant of the block.
   **/
  void genotype(const int &max_cov, const bool haploid)
  {
    string best_geno = "0/0";
    if (haploid)
      best_geno = "0";

    for (uint i = 0; i < variants.size(); ++i)
    {
      Variant *v = &variants[i];

      // If some allele is too covered, assign best_geno (ie 0 or 0/0)
      // with quality 0 and continue
      bool continue_flag = false;
      for (const int &cov : v->coverages)
      {
        if (cov > max_cov)
        {
          GT gt = {best_geno, 0};
          v->add_genotype(gt);
          continue_flag = true;
          continue;
        }
      }
      if (continue_flag)
        continue;

      // The variant wasn't present in any sample: we have only the
      // coverage of the reference allele
      if (v->coverages.size() == 1)
      {
        GT gt = {best_geno, 1};
        v->add_genotype(gt);
        continue;
      }

      // No allele of the variant is covered
      uint total_sum = accumulate(v->coverages.begin(), v->coverages.end(), 0);
      if (total_sum == 0)
      {
        GT gt = {best_geno, 0};
        v->add_genotype(gt);
        continue;
      }

      if (haploid)
      {
        for (uint g1 = 0; g1 < v->coverages.size(); ++g1)
        {
          uint truth = v->coverages[g1];
          uint error = total_sum - truth;

          double log_prior = 2 * log(v->frequencies[g1]);
          double log_posterior = log_binomial(truth + error, truth) +
                                 truth * log(1 - error_rate) +
                                 error * log(error_rate / (v->coverages.size() - 1));

          double log_prob = log_prior + log_posterior;
          double prob = 0;
          if (!isinf(log_prob))
            prob = exp(log_prob);
          GT gen = {to_string(g1), prob};
          v->add_genotype(gen);
        }
      }
      else
      {
        for (uint g1 = 0; g1 < v->coverages.size(); ++g1)
        {
          for (uint g2 = g1; g2 < v->coverages.size(); ++g2)
          {
            double log_prior;
            double log_posterior;
            if (g1 == g2)
            {
              log_prior = 2 * log(v->frequencies[g1]);
              uint truth = v->coverages[g1];
              uint error = total_sum - truth;
              log_posterior = log_binomial(truth + error, truth) +
                              truth * log(1 - error_rate) +
                              error * log(error_rate / (v->coverages.size() - 1));
            }
            else
            {
              log_prior = log(2 * v->frequencies[g1] * v->frequencies[g2]);
              uint truth1 = v->coverages[g1];
              uint truth2 = v->coverages[g2];
              uint error = total_sum - truth1 - truth2;
              log_posterior = log_binomial(truth1 + truth2 + error, truth1 + truth2) +
                              log_binomial(truth1 + truth2, truth1) +
                              truth1 * log((1 - error_rate) / 2) +
                              truth2 * log((1 - error_rate) / 2);
              if (v->coverages.size() > 2)
                // otherwise /0 -> nan
                log_posterior += error * log(error_rate / (v->coverages.size() - 2));
            }

            double log_prob = log_prior + log_posterior;
            double prob = 0;
            if (!isinf(log_prob))
              prob = exp(log_prob);
            GT gen = {to_string(g1) + "/" + to_string(g2), prob};
            v->add_genotype(gen);
          }
        }
      }
    }
  }

  /**
   * Method to output the variants of the block in VCF format.
   * ! Clean this method! Filter is set to "PASS" -
   * see variant.hpp !
   **/
  void output_variants(const bool haploid, const bool verbose)
  {
    for (uint i = 0; i < variants.size(); ++i)
    {
      Variant *v = &variants[i];
      cout << v->seq_name << '\t' << v->ref_pos + 1 << '\t' << v->idx
           << '\t' << v->ref_sub << '\t';
      uint varc = 0;
      for (const string &alt : v->alts)
      {
        cout << alt;
        ++varc;
        if (varc != v->alts.size())
          cout << ',';
      }
      cout << "\t";
      if (isnan(v->quality))
        cout << ".";
      else
        cout << v->quality;
      string info = ".";
      if (verbose)
      {
        // Adds coverages to v->info (here I'm assuming v->info is '.')
        info = "COVS=";
        for (const auto &cov : v->coverages)
          info += to_string((int)cov) + ",";
        info.pop_back();
      }
      // Adds gts to v->info
      string best_geno = haploid ? "0" : "0/0";
      double best_qual = 0.0;
      double total_qual = 0.0;
      for (const auto gt : v->computed_gts)
      {
        total_qual += gt.second;
      }
      string geno = "";
      double qual = 0.0;
      if (verbose)
        info += ";GTS=";
      for (const auto gt : v->computed_gts)
      {
        geno = gt.first;
        qual = gt.second / total_qual;
        if (qual > best_qual)
        {
          best_geno = geno;
          best_qual = qual;
        }
        if (verbose)
          info += geno + ":" + to_string(qual) + ",";
      }
      if (verbose)
        info.pop_back();
      cout << "\t" << v->filter << "\t" << info
           << "\tGT:GQ\t" << best_geno << ":"
           << (int)round(best_qual * 100) << "\n";
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
  bool are_overlapping(const Variant &v1, const Variant &v2) const
  {
    return (v1.ref_pos <= v2.ref_pos) &&
           (v2.ref_pos < v1.ref_pos + v1.ref_size);
  }

  /**
   * Return true if variant v1 is k/2-near (on the right) to variant v2
   **/
  bool are_near(const Variant &v1, const Variant &v2, const int &k,
                const int &sum_to_add = 0) const
  {
    return v1.ref_pos + v1.ref_size - v1.min_size - 1 + sum_to_add +
               ceil((float)k / 2) >=
           v2.ref_pos;
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
  vector<vector<int>> get_combs_on_the_right(const int &i)
  {
    Variant *mid_v = &variants[i];
    vector<vector<int>> right_combs; // possible combinations
    vector<int> right_sums;          // total sum of the combinations (for checking
                                     // (k/2)-proximity)
    bool halt_flag = false;
    for (uint j = i + 1; j < variants.size() && !halt_flag; ++j)
    {
      Variant *curr_v = &variants[j];

      if (!curr_v->is_present)
        continue;

      if (are_overlapping(*mid_v, *curr_v))
        continue;

      if (right_combs.empty())
      { // first var to be added
        if (are_near(*mid_v, *curr_v, k))
        {
          vector<int> new_comb(1, (int)j);
          right_combs.emplace_back(new_comb);
          right_sums.emplace_back(curr_v->ref_size - curr_v->min_size);
        }
      }
      else
      {
        // add the var to all the compatible combinations
        bool added_flag = false;
        for (uint c = 0; c < right_combs.size(); ++c)
        {
          vector<int> comb = right_combs[c];
          int sum = right_sums[c];

          Variant *last_v_in_comb = &variants[comb.back()];

          if (!are_overlapping(*last_v_in_comb, *curr_v))
          {
            added_flag = true;
            if (are_near(*mid_v, *curr_v, k, sum))
            {
              right_combs[c].emplace_back(j);
              right_sums[c] += curr_v->ref_size - curr_v->min_size;
            }
          }
        }
        if (!added_flag)
        {
          // if the var has not been added to any combination
          vector<vector<int>> new_right_combs;
          vector<int> new_right_sums;
          for (uint c = 0; c < right_combs.size(); ++c)
          {
            // shorten the combinations and try to add the var
            vector<int> new_comb = right_combs[c];
            int new_sum = right_sums[c];

            Variant *last_v_in_comb = &variants[new_comb.back()];

            while (are_overlapping(*last_v_in_comb, *curr_v) and
                   !new_comb.empty())
            {
              new_comb.pop_back();
              new_sum -= last_v_in_comb->ref_size - last_v_in_comb->min_size;
              last_v_in_comb = &variants[new_comb.back()];
            }
            new_comb.emplace_back(j);
            if (are_near(*mid_v, *curr_v, k, new_sum))
            {
              added_flag = true;
              new_right_combs.emplace_back(new_comb);
              new_right_sums.emplace_back(new_sum + curr_v->ref_size -
                                          curr_v->min_size);
            }
          }
          extend(right_combs, new_right_combs);
          extend(right_sums, new_right_sums);

          // if the var has not been added to any combination (neither to
          // shortened ones), we can halt the loop: no successive variants
          // will be added (too far away)
          if (!added_flag)
            halt_flag = true;
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
  vector<vector<int>> get_combs_on_the_left(const int &i)
  {
    Variant *mid_v = &variants[i];
    vector<vector<int>> left_combs; // possible combinations
    vector<int> left_sums;          // total sum of the combinations (for checking
                                    // (k/2)-proximity)

    bool halt_flag = false;
    for (int j = i - 1; j >= 0 && !halt_flag; --j)
    {
      Variant *curr_v = &variants[j];

      if (!curr_v->is_present)
        continue;

      if (are_overlapping(*curr_v, *mid_v))
        continue;

      if (left_combs.empty())
      { // first var to be added
        if (are_near(*curr_v, *mid_v, k))
        {
          vector<int> new_comb(1, (int)j);
          left_combs.emplace_back(new_comb);
          left_sums.emplace_back(curr_v->ref_size - curr_v->min_size);
        }
      }
      else
      {
        // add the var to all the compatible combinations
        bool added_flag = false;
        for (uint c = 0; c < left_combs.size(); ++c)
        {
          vector<int> comb = left_combs[c];
          int sum = left_sums[c];

          Variant *last_v_in_comb = &variants[comb.back()];

          if (!are_overlapping(*curr_v, *last_v_in_comb))
          {
            added_flag = true;
            if (are_near(*curr_v, *mid_v, k, sum))
            {
              left_combs[c].emplace_back(j);
              left_sums[c] += curr_v->ref_size - curr_v->min_size;
            }
          }
        }
        if (!added_flag)
        {
          // if the var has not been added to any combination
          vector<vector<int>> new_left_combs;
          vector<int> new_left_sums;
          for (uint c = 0; c < left_combs.size(); ++c)
          {
            // shorten the combinations and try to add the var
            vector<int> new_comb = left_combs[c];
            int new_sum = left_sums[c];

            Variant *last_v_in_comb = &variants[new_comb.back()];

            while (are_overlapping(*curr_v, *last_v_in_comb) &&
                   !new_comb.empty())
            {
              new_comb.pop_back();
              new_sum -= last_v_in_comb->ref_size - last_v_in_comb->min_size;
              last_v_in_comb = &variants[new_comb.back()];
            }
            new_comb.emplace_back(j);
            if (are_near(*curr_v, *mid_v, k, new_sum))
            {
              added_flag = true;
              new_left_combs.emplace_back(new_comb);
              new_left_sums.emplace_back(new_sum + curr_v->ref_size -
                                         curr_v->min_size);
            }
          }
          extend(left_combs, new_left_combs);
          extend(left_sums, new_left_sums);

          // if the var has not been added to any combination (neither to
          // shortened ones), we can halt the loop: no successive variants
          // will be added (too far away)
          if (!added_flag)
            halt_flag = true;
        }
      }
    }

    return left_combs;
  }

  /**
   * Combine all left_combs with all right_combs, and placing i-th variant
   *between them
   **/
  vector<vector<int>>
  combine_combs(vector<vector<int>> &left_combs,
                const vector<vector<int>> &right_combs,
                const int &i)
  {
    vector<vector<int>> full_combs;
    if (left_combs.empty() && right_combs.empty())
    {
      vector<int> comb;
      comb.emplace_back(i);
      full_combs.emplace_back(comb);
    }
    else if (left_combs.empty())
    {
      vector<int> lcomb;
      lcomb.emplace_back(i);
      vector<int> comb = lcomb;
      for (const vector<int> &rcomb : right_combs)
      {
        extend(comb, rcomb);
        full_combs.emplace_back(comb);
        comb = lcomb;
      }
    }
    else
    {
      for (vector<int> &lcomb : left_combs)
      {
        reverse(lcomb.begin(), lcomb.end());
        lcomb.emplace_back(i);
        vector<int> comb = lcomb;
        if (right_combs.empty())
        {
          full_combs.emplace_back(comb);
        }
        else
        {
          for (const vector<int> &rcomb : right_combs)
          {
            extend(comb, rcomb);
            full_combs.emplace_back(comb);
            comb = lcomb;
          }
        }
      }
    }
    return full_combs;
  }

  /**
   * Return the reference substring between considered variants
   **/
  vector<string> get_ref_subs(const vector<int> &comb, const string &reference)
  {
    vector<string> ref_subs;

    int last_end = -1;
    for (const int &index : comb)
    {
      Variant *v = &variants[index];
      if (last_end == -1)
      {
        last_end = v->ref_pos + v->ref_size;
        continue;
      }
      string ref_sub(reference,
                     last_end,
                     v->ref_pos - last_end);
      ref_subs.emplace_back(ref_sub);
      last_end = v->ref_pos + v->ref_size;
    }
    return ref_subs;
  }

  /**
   * Given two haplotypes, build all possible combinations from them
   * (used when gt information is unphased). Example: given 0/1, 1/3,
   * 1/1, I want: 0,1,1; 0,3,1; 1,1,1; 1,3,1 (with repetitions).
   **/
  vector<vector<string_view>> combine_haplotypes(const vector<string_view> &hap1,
                                                 const vector<string_view> &hap2)
  {
    int n = hap1.size();   // number of alleles in each haplotype
    int N = pow(2, n - 1); // number of possible haplotypes
    vector<vector<string_view>> HAPs(2 * N, vector<string_view>(n));

    for (int level = 0; level < n; ++level)
    {
      vector<string_view> alleles({hap1[level], hap2[level]});
      int rep = pow(2, n - 1 - level);

      for (int col = 0; col < N; ++col)
      {
        HAPs[col][level] = alleles[(col / rep) % 2];
        HAPs[col + N][level] = alleles[(col / rep + 1) % 2];
      }
    }
    return HAPs;
  }

  /**
   * Builds and returns all the possible combination of alleles (haplotypes),
   * with respect to GTs.
   **/
  unordered_set<vector<string_view>, VectorHash>
  build_alleles_combs(const vector<int> &comb,
                      const int central_index,
                      const bool haploid)
  {
    // A set to avoid duplicate elements
    unordered_set<vector<string_view>, VectorHash> aacs;
    Variant *central_v = &variants[central_index];
    // For each individual having this variant
    for (int gt_i = 0; gt_i < (int)central_v->genotypes.size(); ++gt_i)
    {
      if (haploid)
      {
        vector<string_view> hap1;
        hap1.reserve(comb.size());
        for (const int j : comb)
        {
          hap1.emplace_back(variants[j].get_allele(variants[j].genotypes[gt_i].first));
        }

        aacs.insert(hap1);
      }
      else
      {
        bool phased_combination = true;
        vector<string_view> hap1;
        hap1.reserve(comb.size());
        vector<string_view> hap2;
        hap2.reserve(comb.size());
        for (const int j : comb)
        {
          phased_combination &= variants[j].phasing[gt_i];
          hap1.emplace_back(variants[j].get_allele(variants[j].genotypes[gt_i].first));
          hap2.emplace_back(variants[j].get_allele(variants[j].genotypes[gt_i].second));
        }

        if (phased_combination)
        {
          aacs.insert(hap1);
          aacs.insert(hap2);
        }
        else
        {
          vector<vector<string_view>> all_haplotypes = combine_haplotypes(hap1, hap2);
          for (const vector<string_view> &hap : all_haplotypes)
          {
            aacs.insert(hap);
          }
        }
      }
    }
    return aacs;
  }

  /**
   * Log of the binomial coefficient computed using Stirling
   * approximation (see https://en.wikipedia.org/wiki/Stirling's_approximation)
   **/
  double log_binomial(const int &n, const int &k)
  {
    if (n == 0 || n == k || k == 0)
      return 0;
    return n * log(n) - k * log(k) - (n - k) * log(n - k);
  }
};

#endif