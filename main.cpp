#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <getopt.h>
#include <omp.h>
#include <zlib.h>

#include "bloom_filter.hpp"
#include "hts_log.h"
#include "kmc_api/kmc_file.h"
#include "kseq.h"
#include "vcf.h"

KSEQ_INIT(gzFile, gzread)

using namespace std;

auto start_t = chrono::high_resolution_clock::now();

void pelapsed(const string &s = "") {
  auto now_t = chrono::high_resolution_clock::now();
  cerr << "[malva/" << s << "] Time elapsed "
       << chrono::duration_cast<chrono::milliseconds>(now_t - start_t).count()
       << endl;
}

typedef map<int, vector<pair<int, vector<string>>>>
    VK_GROUP; // VARIANT - KMERS GROUP

static const char *USAGE_MESSAGE =
    "Usage: malva [-k KMER-SIZE] [-c MIN-COV] [-x MAX-COV] <reference.fa> "
    "<variants.vcf> "
    "<kmc_output_prefix>\n"
    "\n"
    "Top notch description of this tool\n"
    "\n"
    "      -h, --help                        display this help and exit\n"
    "      -l, --loose                       loose mode. If set then the "
    "reference kmers will not be removed (default: false)\n"
    "      -t, --threads                     number of threads (default: 1)\n"
    "\n"
    "      -k, --kmer-size                   size of the kmers to index "
    "(default:15)\n"
    "      -c, --min-coverage                minimum coverage for validated "
    "variants (default:30)\n"
    "      -x, --max-coverage                maximum coverage for validated "
    "variants.\n"
    "                                        If this value is greter than 250 "
    "then no upper bound is set (default: unlimited)\n"
    "      -u, --max-combinations            maximum (even) number of neighbouring "
    "variants considered for each signature (default:4)"
    "\n";

namespace opt {
static bool strictMode = true;
static size_t klen = 15;
static size_t ref_klen = 31;
static size_t vcfBlockSize = 10000;
static size_t bfSize = ((size_t)0b1 << 34);
static size_t maxCombinations = 4;
static size_t minCoverage = 30;
static size_t maxCoverage = 251;
static size_t nThreads = 1;
static string fastaPath;
static string vcfPath;
static string sequencePath;
} // namespace opt

static const char *shortopts = "k:r:u:c:x:b:t:lh";

static const struct option longopts[] = {
    {"kmer-size", required_argument, NULL, 'k'},
    {"ref-kmer", required_argument, NULL, 'r'},
    {"max-combinations", required_argument, NULL, 'u'},
    {"min-coverage", required_argument, NULL, 'c'},
    {"max-coverage", required_argument, NULL, 'x'},
    {"bf-size", required_argument, NULL, 'b'},
    {"loose", no_argument, NULL, 'l'},
    {"threads", no_argument, NULL, 't'},
    {"help", no_argument, NULL, 'h'},
    {NULL, 0, NULL, 0}};

void parseOptions(int argc, char **argv) {
  bool die = false;
  for (char c;
       (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'b':
      // Let's consider this as MB
      arg >> opt::bfSize;
      opt::bfSize = opt::bfSize * (0b1 << 20);
      break;
    case 'c':
      arg >> opt::minCoverage;
      break;
    case 'x':
      arg >> opt::maxCoverage;
      break;
    case 'u':
      arg >> opt::maxCombinations;
      break;
    case 'k':
      arg >> opt::klen;
      break;
    case 'r':
      arg >> opt::ref_klen;
      break;
    case 'l':
      opt::strictMode = false;
      break;
    case 't':
      arg >> opt::nThreads;
      break;
    case '?':
      die = true;
      break;
    case 'h':
      cout << USAGE_MESSAGE;
      exit(EXIT_SUCCESS);
    }
  }

  if (argc - optind < 3) {
    cerr << "malva : missing arguments\n";
    die = true;
  } else if (argc - optind > 3) {
    cerr << "malva : too many arguments\n";
    die = true;
  }

  if (die) {
    cerr << "\n" << USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }

  opt::fastaPath = argv[optind++];
  opt::vcfPath = argv[optind++];
  opt::sequencePath = argv[optind++];
}

struct Variant {
  string seqName;
  int refPos;          // Variant position 0-based
  string refSub;       // Reference base{s}
  vector<string> alts; // List of alternatives
  int minLength;       // Length of the shortest string at that position
  bool goodVar = true; // false if no variants, i.e. only <CN>

  Variant() {}

  Variant(const string &_seqName, const int _refPos, const char *_refSub,
          char **_alts, const size_t _altsNum)
      : seqName(_seqName), refPos(_refPos), refSub(_refSub) {
    for (size_t i = 0; i < _altsNum; ++i) {
      if (_alts[i][0] != '<')
        alts.push_back(string(_alts[i]));
    }
    setMinSize();
  }

  vector<string> getAllStrings() const {
    vector<string> allStrings = alts;
    allStrings.push_back(refSub);
    return allStrings;
  }

  void setMinSize() {
    if (alts.size() == 0)
      goodVar = false;
    else {
      minLength = alts[0].length();
      for (size_t i = 1; i < alts.size(); i++) {
        if ((int)alts[i].length() < minLength)
          minLength = alts[i].length();
      }
      minLength = min(minLength, (int)refSub.length());
    }
  }

  void print() const {
    cout << refPos << " " << refSub << " ";
    for (const string alt : alts)
      cout << alt << "-";
    cout << endl;
  }
};

// Extract kmer/s from single variant
vector<pair<int, vector<string>>> extractKmers(const char *reference,
                                               const Variant &var) {
  vector<pair<int, vector<string>>> KMERs;

  int altID = 0;
  for (auto &alt : var.alts) {
    vector<string> kmers;

    int altSize = alt.length();
    if (altSize <= (int)opt::klen) {
      int prefStart = var.refPos - opt::klen / 2 + altSize / 2;
      int borderLen = (opt::klen - altSize) / 2;
      string pref(reference + prefStart, borderLen);

      // we take a char after the alternative and not before
      if (altSize % 2 == 0)
        ++borderLen;
      int suffStart = var.refPos + var.refSub.length();
      string suff(reference + suffStart, borderLen);
      string kmer = pref + alt + suff;
      kmers.push_back(kmer);
    } else { // the alternative is longer than k
      int i = 0;
      while (i < altSize - (int)opt::klen + 1) {
        string kmer = alt.substr(i, opt::klen);
        kmers.push_back(kmer);
        ++i;
      }
    }
    KMERs.push_back(make_pair(altID, kmers));
    ++altID;
  }
  return KMERs;
}

/**
   Given two variants, return true if they overlap
 **/
bool overlapping(const Variant &var1, const Variant &var2) {
  return (var1.refPos <= var2.refPos) &&
         (var2.refPos < var1.refPos + (int)var1.refSub.size());
}

/**
   Given a variations block and a center i, we extract the [leftIndex,
 rightIndex]
   window containing variants "near" the center.
   We do not check for overlapping variants: !!we can compute bigger windows!!,
   anyway we cut the kmers.
 **/
pair<int, int> getWindowIndices(vector<Variant> &varBlock, const int &i) {
  bool leftWithin = true;
  bool rightWithin = true;
  int leftIndex = i - 1;
  int rightIndex = i + 1;

  // on the left we do not consider the size of the center, we start from its
  // begin
  int leftSum = 0;
  // on the rigth we consider the size of the center, since we start from its
  // begin
  int rightSum = varBlock[i].refSub.size() - varBlock[i].minLength;

  while (leftWithin || rightWithin) {
    if (leftWithin) {
      if (leftIndex >= 0) {
        leftSum +=
            varBlock[leftIndex].refSub.size() - varBlock[leftIndex].minLength;
        if (varBlock[i].refPos - varBlock[leftIndex].refPos <=
            (int)opt::klen / 2 + leftSum) {
          --leftIndex;
        } else {
          leftWithin = false;
          ++leftIndex;
        }
      } else {
        leftIndex = 0;
        leftWithin = false;
      }
    }

    if (rightWithin) {
      if (rightIndex < (int)varBlock.size()) {
        rightSum +=
            varBlock[rightIndex].refSub.size() - varBlock[rightIndex].minLength;
        if (varBlock[rightIndex].refPos - varBlock[i].refPos <
            ceil((float)opt::klen / 2) + rightSum) {
          ++rightIndex;
        } else {
          rightWithin = false;
          --rightIndex;
        }
      } else {
        rightIndex = varBlock.size() - 1;
        rightWithin = false;
      }
    }
  }
  return make_pair(leftIndex, rightIndex);
}

/**
   Extend a vector with another
 **/
template <typename T> void extend(T &V1, const T &V2) {
  V1.reserve(V1.size() + distance(V2.begin(), V2.end()));
  V1.insert(V1.end(), V2.begin(), V2.end());
}

/**
   Compute all the Q-combinations of a list, where Q is opt::maxComination
   --- adapted from rosettacode
   If rightFlag is true, the first element of the list is fixed; the last
   otherwise (ie the center must be fixed).
**/
void computeQCombinations(vector<int> L, vector<vector<int>> &qCombinations) {
  int N = L.size()-1; //first or last are fixed
  std::string bitmask(opt::maxCombinations/2, 1); // Q leading 1's
  bitmask.resize(N, 0); // N-Q trailing 0's
  do {
    vector<int> comb ({L[0]});
    for (int i = 0; i < N; ++i) {
      if (bitmask[i])
        comb.push_back(L[i+1]);
    }
    qCombinations.push_back(comb);
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
}

/**
   Compute all the possible sequences of variants which can extend the center on
 the right
   (here we check for overlapping variants).
 **/
vector<vector<int>>
getPossibleCombinationOnTheRight(const vector<Variant> &varBlock,
                                 const int &i,
                                 const int &rightIndex) {
  vector<vector<int>> rightValidVariants(1, {i});
  int currI = i + 1;
  /**
  # Combining the center only with the last maxCombinations of the window
  if(opt::maxCombinations>0) {
    currI = rightIndex-opt::maxCombinations+1;
    if(currI<i+1)
      currI = i + 1;
  }
  **/
  while (currI <= rightIndex) {
    bool alreadyAddedFlag = false;
    for (auto &validVarsSequence : rightValidVariants) {
      int lastVar = validVarsSequence.back();
      if (!overlapping(varBlock[lastVar], varBlock[currI])) {
        validVarsSequence.push_back(currI);
        alreadyAddedFlag = true;
      }
    }
    if (!alreadyAddedFlag) {
      if (!overlapping(varBlock[i], varBlock[currI])) {
        vector<int> validVarsSequence({i, currI});
        rightValidVariants.push_back(validVarsSequence);
      }
    }
    ++currI;
  }
  vector<vector<int>> rightVariants;
  for(const auto &vars : rightValidVariants) {
    if(vars.size() <= opt::maxCombinations/2+1)
      rightVariants.push_back(vars);
    else {
      //Getting all Q-combinations, with Q=opt::maxCombinations, such
      //that the first element is fixed (the center)
      vector<vector<int>> qCombinations;
      computeQCombinations(vars, qCombinations);
      extend(rightVariants, qCombinations);
    }
  }
  return rightVariants;
}

/**
   Compute all the possible sequences of variants which can extend the center on
 the left
   (here we check for overlapping variants).
 **/
vector<vector<int>>
getPossibleCombinationOnTheLeft(const vector<Variant> &varBlock, const int &i,
                                const int &leftIndex) {
  vector<vector<int>> leftValidVariants(1, {i});

  int currI = i - 1;
  /**
  # Combining the center only with the first maxCombinations of the window
  if(opt::maxCombin>0) {
    currI = leftIndex+opt::maxCombin-1;
    if(currI>i-1)
      currI = i - 1;
  }
  **/
  while (leftIndex <= currI) {
    bool alreadyAddedFlag = false;
    for (auto &validVarsSequence : leftValidVariants) {
      int lastVar = validVarsSequence.back();
      if (!overlapping(varBlock[currI], varBlock[lastVar])) {
        validVarsSequence.push_back(currI);
        alreadyAddedFlag = true;
      }
    }
    if (!alreadyAddedFlag) {
      if (!overlapping(varBlock[currI], varBlock[i])) {
        vector<int> validVarsSequence({i, currI});
        leftValidVariants.push_back(validVarsSequence);
      }
    }
    --currI;
  }
  vector<vector<int>> leftVariants;
  for(const auto &vars : leftValidVariants) {
    if(vars.size() <= opt::maxCombinations/2+1)
      leftVariants.push_back(vars);
    else {
      //Getting all Q-combinations, with Q=opt::maxCombinations, such
      //that the last element is fixed (the center)
      vector<vector<int>> qCombinations;
      computeQCombinations(vars, qCombinations);
      extend(leftVariants, qCombinations);
    }
  }
  return leftVariants;
}

/**
   Extend on the right the kmers of a vector to reach length k/2
 **/
void extendOnTheRight(vector<string> &kmers, const char *reference,
                      const string &alternative, const Variant &var) {
  int firstStartAltSize = alternative.size() / 2;
  for (string &kmer : kmers) {
    if (kmer.size() - firstStartAltSize < ceil((float)opt::klen / 2)) {
      int startPos = var.refPos + var.refSub.size();
      int len = opt::klen;
      string postBorder(reference + startPos, len);

      kmer = kmer + postBorder;
      kmer = kmer.substr(0, ceil((float)opt::klen / 2) + firstStartAltSize);
    }
  }
}

/**
   Compute the kmers extending the considered alternative of the center on the
 right
 **/
vector<string>
computeKmersOnRightPortion(const char *reference,
                           const vector<Variant> &varBlock,
                           const vector<vector<int>> &rightValidVariants,
                           const string &alternative, const int &centerI) {
  vector<string> kmers;
  int firstStartAltSize = alternative.size() / 2;
  for (auto &rightValidSeq : rightValidVariants) {
    vector<string> varExtKmers;
    varExtKmers.push_back(alternative);
    int lastI = centerI;
    for (auto &currI : rightValidSeq) {
      if (currI != centerI) {
        vector<string> tmpKmers;
        vector<string> currAlts = varBlock[currI].getAllStrings();

        int startPos = varBlock[lastI].refPos + varBlock[lastI].refSub.size();
        int len = varBlock[currI].refPos -
                  (varBlock[lastI].refPos + varBlock[lastI].refSub.size());
        string refInsideVars(reference + startPos, len);

        for (const auto &kmer : varExtKmers) {
          for (const auto &currAlt : currAlts) {
            string tmp = kmer + refInsideVars + currAlt;
            if (tmp.size() - firstStartAltSize > ceil((float)opt::klen / 2)) {
              tmp =
                  tmp.substr(0, ceil((float)opt::klen / 2) + firstStartAltSize);
            }
            tmpKmers.push_back(tmp);
          }
        }
        varExtKmers = tmpKmers;
      }
      lastI = currI;
    }

    Variant lastVar = varBlock[rightValidSeq.back()];
    extendOnTheRight(varExtKmers, reference, alternative, lastVar);

    extend(kmers, varExtKmers);
  }
  return kmers;
}

/**
   Extend on the left the kmers of a vector to reach length k
 **/
void reachKOnTheLeft(vector<string> &kmers, const char *reference,
                     const Variant &var) {
  for (auto &kmer : kmers) {
    if (kmer.size() < opt::klen) {
      int startPos = max(0, var.refPos - (int)opt::klen);
      int len = var.refPos - startPos;
      string preBorder(reference + startPos, len);
      kmer = preBorder + kmer;
      kmer = kmer.substr(kmer.size() - opt::klen, opt::klen);
    }
  }
}

/**
   Compute the kmers extending the considered (half)kmer on the left
 **/
vector<string> extendOnTheLeft(const char *reference,
                               const vector<Variant> &varBlock,
                               const vector<int> &leftValidSeq,
                               const string &halfKmer, const int &centerI) {
  int lastI = centerI;
  vector<string> varExtKmers;
  varExtKmers.push_back(halfKmer);
  for (auto &currI : leftValidSeq) {
    if (currI != centerI) {
      vector<string> tmpKmers;
      vector<string> currStrings = varBlock[currI].getAllStrings();

      int startPos = varBlock[currI].refPos + varBlock[currI].refSub.size();
      int len = varBlock[lastI].refPos -
                (varBlock[currI].refPos + varBlock[currI].refSub.size());
      string refInsideVars(reference + startPos, len);

      for (const auto &kmer : varExtKmers) {
        for (const auto &currAlt : currStrings) {
          string tmp = currAlt + refInsideVars + kmer;
          if (tmp.size() > opt::klen) {
            tmp = tmp.substr(tmp.size() - opt::klen, opt::klen);
          }
          tmpKmers.push_back(tmp);
        }
      }
      varExtKmers = tmpKmers;
    }
    lastI = currI;
  }

  Variant firstVar = varBlock[leftValidSeq.back()];
  reachKOnTheLeft(varExtKmers, reference, firstVar);

  return varExtKmers;
}

VK_GROUP
extractKmersFromBlock(const char *reference, vector<Variant> &varBlock,
                      const size_t fromPos, const size_t toPos) {

  VK_GROUP KMERs;

  omp_lock_t writelock;
  omp_init_lock(&writelock);

#pragma omp parallel for shared(reference) shared(varBlock) shared(KMERs)
  for (size_t i = fromPos; i < toPos; ++i) {
    vector<pair<int, vector<string>>> kmers;

    // Computing window centered on the considered center (i)
    pair<int, int> windowIndices = getWindowIndices(varBlock, i);
    int leftIndex = windowIndices.first;
    int rightIndex = windowIndices.second;

    if (leftIndex == rightIndex || opt::maxCombinations == 0) { // Only the center in the window
      vector<pair<int, vector<string>>> altID2currKmers =
          extractKmers(reference, varBlock[i]);
      extend(kmers, altID2currKmers);
    } else {
      // All possible combination to extend the center on the right and on the
      // left
      vector<vector<int>> rightValidVariants =
          getPossibleCombinationOnTheRight(varBlock, i, rightIndex);
      vector<vector<int>> leftValidVariants =
          getPossibleCombinationOnTheLeft(varBlock, i, leftIndex);

      // If no combination (due to overlaps) -> only the center
      if (rightValidVariants.size() == 1 && rightValidVariants[0].size() == 1 &&
          leftValidVariants.size() == 1 && leftValidVariants[0].size() == 1) {
        vector<pair<int, vector<string>>> altID2currKmers =
            extractKmers(reference, varBlock[i]);
        extend(kmers, altID2currKmers);
      } else {
        // Extend the center on the right
        vector<pair<int, vector<string>>> altID2currKmers;
        int altID = 0;
        // FIXME: create method getVariants
        vector<string> startingStrings = varBlock[i].getAllStrings();
        startingStrings.pop_back();
        for (const auto &startAlt : startingStrings) {
          vector<string> currKmers;
          vector<string> halfKmers = computeKmersOnRightPortion(
              reference, varBlock, rightValidVariants, startAlt, i);
          // Extend the (half) kmers on the left
          for (auto &halfKmer : halfKmers) {
            for (auto &leftValidSeq : leftValidVariants) {
              vector<string> varExtKmers = extendOnTheLeft(
                  reference, varBlock, leftValidSeq, halfKmer, i);
              extend(currKmers, varExtKmers);
            }
          }
          for (const auto &kmer : currKmers) {
            vector<string> TMP({kmer});
            altID2currKmers.push_back(make_pair(altID, TMP));
          }
          ++altID;
        }
        extend(kmers, altID2currKmers);
      }
    }
    omp_set_lock(&writelock);
    KMERs.insert(make_pair(i, kmers));
    omp_unset_lock(&writelock);
  }

  return KMERs;
}

int sumKmers(const VK_GROUP &KMERs) {
  int N = 0;
  for (VK_GROUP::const_iterator it = KMERs.begin(); it != KMERs.end(); ++it) {
    for (const auto &p : it->second) {
      for (const auto k : p.second) {
        ++N;
      }
    }
  }
  return N;
}

void printKmers(const VK_GROUP &KMERs) {
  for (VK_GROUP::const_iterator it = KMERs.begin(); it != KMERs.end(); ++it) {
    int vID = it->first;
    cout << vID << endl;
    for (const auto &p : it->second) {
      int altID = p.first;
      cout << "\t" << altID << endl;
      for (const auto k : p.second) {
        cout << "\t\t" << k << endl;
        // cout << vID << " " << altID << " " << k << endl;
      }
    }
  }
}

void addKmersToBF(const VK_GROUP &vk, BF &bf) {
  for (const auto &variant_line : vk) {
    for (const auto &alt_kmers : variant_line.second) {
      for (const auto &kmer : alt_kmers.second) {
        bf.add_key(kmer);
      }
    }
  }
}

static size_t number_variants_out = 0;

// TODO: print in VCF
void outputGoodVariants(const bcf_hdr_t *hdr,
                        const map<int, vector<int>> &good_variants,
                        const vector<Variant> &vars) {
  for (const auto &elem : good_variants) {
    int var_idx = elem.first;
    cout << vars[var_idx].seqName << '\t' << vars[var_idx].refPos + 1 << '\t'
         << number_variants_out++ << '\t' << vars[var_idx].refSub << '\t';
    size_t varc = 0;
    for (const int alt_id : elem.second) {
      cout << vars[var_idx].alts[alt_id];
      ++varc;
      if (varc != elem.second.size())
        cout << ',';
    }
    cout << "\t100\tPASS\tVT=SNP\n";
  }
}

// FIXME: do we need this method?
bool contains(vector<int> &V, const int &i) {
  for (const auto &v : V) {
    if (v == i)
      return true;
  }
  return false;
}

int main(int argc, char *argv[]) {
  hts_set_log_level(HTS_LOG_OFF);

  parseOptions(argc, argv);

  omp_set_dynamic(0);
  omp_set_num_threads(opt::nThreads);

  gzFile fastain = gzopen(opt::fastaPath.c_str(), "r");
  htsFile *vcf = bcf_open(opt::vcfPath.c_str(), "r");
  kseq_t *reference = kseq_init(fastain);
  kseq_read(reference);
  bcf_hdr_t *vcf_header = bcf_hdr_read(vcf);

  kstring_t htxt = {0, 0, 0};
  bcf_hdr_format(vcf_header, 0, &htxt);
  cout << htxt.s;
  free(htxt.s);

  // Reading VCF file
  vector<Variant> vars;
  size_t begin_p = 0; // first kmer to analyze
  size_t end_p = 0;   // last var to analyze is var[end_p - 1]
  BF bf(opt::bfSize);

  // Step 1: add vcf kmers to bloom filter
  bcf1_t *vcf_record = bcf_init();
  while (bcf_read(vcf, vcf_header, vcf_record) == 0) {
    bcf_unpack(vcf_record, BCF_UN_STR);
    Variant v(bcf_hdr_id2name(vcf_header, vcf_record->rid), vcf_record->pos,
              vcf_record->d.allele[0], vcf_record->d.allele + 1,
              vcf_record->n_allele - 1);
    if (v.alts.size() > 0) {
      vars.push_back(v);
      ++end_p;
      if (vars.size() >= opt::vcfBlockSize) {
        // Fill in variants that are "close" to vars[complete_kmers_till - 1]
        while (bcf_read(vcf, vcf_header, vcf_record) == 0 &&
               vcf_record->pos - vars[end_p - 1].refPos <=
                   ceil(static_cast<float>(opt::klen) / 2)) {
          bcf_unpack(vcf_record, BCF_UN_STR);
          Variant v(bcf_hdr_id2name(vcf_header, vcf_record->rid),
                    vcf_record->pos, vcf_record->d.allele[0],
                    vcf_record->d.allele + 1, vcf_record->n_allele - 1);
          if (v.alts.size() > 0)
            vars.push_back(v);
        }
        // FIXME: I think we are skipping one variant everytime we get here
        VK_GROUP variants_kmers =
            extractKmersFromBlock(reference->seq.s, vars, begin_p, end_p);
        addKmersToBF(variants_kmers, bf);
        size_t new_begin = end_p - 1;
        while (new_begin != 0 &&
               vars[end_p - 1].refPos - vars[new_begin].refPos <=
                   ceil(static_cast<float>(opt::klen) / 2)) {
          --new_begin;
        }
        end_p -= new_begin;
        begin_p = end_p;
        vars.erase(vars.begin(), vars.begin() + new_begin);
      }
    }
  }
  if (vars.size() > 0) {
    VK_GROUP variants_kmers =
        extractKmersFromBlock(reference->seq.s, vars, begin_p, vars.size());
    addKmersToBF(variants_kmers, bf);
  }
  bcf_hdr_destroy(vcf_header);
  bcf_destroy(vcf_record);
  bcf_close(vcf);

  BF refbf(opt::bfSize);

  if (opt::strictMode) {
    for (size_t p = (opt::ref_klen - opt::klen) / 2;
         p < reference->seq.l - opt::ref_klen; ++p) {
      string refk(reference->seq.s + p, opt::klen);
      if (bf.test_key(refk)) {
        string context(reference->seq.s + p - ((opt::ref_klen - opt::klen) / 2),
                       opt::ref_klen);
        refbf.add_key(context);
      }
    }
  }

  kseq_destroy(reference);
  gzclose(fastain);

  bf.switch_mode();

  pelapsed("BF creation complete");

  // Step 2: test variants present in sequence fasta
  {
    CKMCFile kmer_db;
    uint32 klen, mode, min_counter, pref_len, sign_len, min_c, counter;
    uint64 tot_kmers, max_c;
    if (!kmer_db.OpenForListing(opt::sequencePath)) {
      cerr << "ERROR: cannot open " << opt::sequencePath << endl;
      return 1;
    }
    kmer_db.Info(klen, mode, min_counter, pref_len, sign_len, min_c, max_c,
                 tot_kmers);
    CKmerAPI kmer_obj(klen);
    string context;
    while (kmer_db.ReadNextKmer(kmer_obj, counter)) {
      kmer_obj.to_string(context);
      if (!refbf.test_key(context)) {
        string kmer(context.c_str() + ((opt::ref_klen - opt::klen) / 2),
                    opt::klen);
        if (bf.test_key(kmer))
          bf.increment(kmer, counter);
      }
    }
  }

  pelapsed("BF weights created");

  // Step 3: check if variants in vcf are covered enough
  fastain = gzopen(opt::fastaPath.c_str(), "r");
  vcf = bcf_open(opt::vcfPath.c_str(), "r");
  reference = kseq_init(fastain);
  kseq_read(reference);
  vcf_header = bcf_hdr_read(vcf);

  // Reading VCF file
  vars.clear();

  vcf_record = bcf_init();

  begin_p = 0;
  end_p = 0;
  while (bcf_read(vcf, vcf_header, vcf_record) == 0) {
    bcf_unpack(vcf_record, BCF_UN_STR);
    Variant v(bcf_hdr_id2name(vcf_header, vcf_record->rid), vcf_record->pos,
              vcf_record->d.allele[0], vcf_record->d.allele + 1,
              vcf_record->n_allele - 1);
    if (v.alts.size() > 0) {
      vars.push_back(v);
      ++end_p;
      if (vars.size() >= opt::vcfBlockSize) {
        // Fill in variants that are "close" to vars[complete_kmers_till - 1]
        while (bcf_read(vcf, vcf_header, vcf_record) == 0 &&
               vcf_record->pos - vars[end_p - 1].refPos <=
                   ceil(static_cast<float>(opt::klen) / 2)) {
          bcf_unpack(vcf_record, BCF_UN_STR);
          Variant v(bcf_hdr_id2name(vcf_header, vcf_record->rid),
                    vcf_record->pos, vcf_record->d.allele[0],
                    vcf_record->d.allele + 1, vcf_record->n_allele - 1);
          if (v.alts.size() > 0)
            vars.push_back(v);
        }
        // FIXME: I think we are skipping one variant everytime we get here
        VK_GROUP variants_kmers =
            extractKmersFromBlock(reference->seq.s, vars, begin_p, end_p);
        map<int, vector<int>> good_variants;
        for (const auto &var_line : variants_kmers) {
          const auto var_idx = var_line.first;
          for (const auto &alt_pair : var_line.second) {
            const auto alt_idx = alt_pair.first;
            size_t kpt = 0; // kmers passing threshold
            for (const auto &kmer : alt_pair.second) {
              if (bf.get_count(kmer) >= opt::minCoverage &&
                  bf.get_count(kmer) <= opt::maxCoverage)
                ++kpt;
            }
            if (kpt == alt_pair.second.size()) {
              if (good_variants.find(var_idx) == good_variants.end()) {
                good_variants[var_idx] = vector<int>();
              }
              // FIXME: is this check right?
              if (!contains(good_variants[var_idx], alt_idx)) {
                good_variants[var_idx].push_back(alt_idx);
              }
            }
          }
        }
        outputGoodVariants(vcf_header, good_variants, vars);
        size_t new_begin = end_p - 1;
        while (new_begin != 0 &&
               vars[end_p - 1].refPos - vars[new_begin].refPos <=
                   ceil(static_cast<float>(opt::klen) / 2)) {
          --new_begin;
        }
        end_p -= new_begin;
        begin_p = end_p;
        vars.erase(vars.begin(), vars.begin() + new_begin);
      }
    }
  }
  if (vars.size() > 0) {
    VK_GROUP variants_kmers =
        extractKmersFromBlock(reference->seq.s, vars, begin_p, end_p);
    map<int, vector<int>> good_variants;
    for (const auto &var_line : variants_kmers) {
      const auto var_idx = var_line.first;
      for (const auto &alt_pair : var_line.second) {
        const auto alt_idx = alt_pair.first;
        size_t kpt = 0; // kmers passing threshold
        for (const auto &kmer : alt_pair.second) {
          if (bf.get_count(kmer) >= opt::minCoverage &&
              bf.get_count(kmer) <= opt::maxCoverage)
            ++kpt;
        }
        if (kpt == alt_pair.second.size()) {
          if (good_variants.find(var_idx) == good_variants.end()) {
            good_variants[var_idx] = vector<int>();
          }
          // FIXME: is this check right?
          if (!contains(good_variants[var_idx], alt_idx)) {
            good_variants[var_idx].push_back(alt_idx);
          }
        }
      }
    }
    outputGoodVariants(vcf_header, good_variants, vars);
  }

  bcf_hdr_destroy(vcf_header);
  bcf_destroy(vcf_record);
  bcf_close(vcf);
  kseq_destroy(reference);
  gzclose(fastain);

  cout.flush();

  pelapsed("Execution completed");

  return 0;
}
