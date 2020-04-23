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

#ifndef _MALVAARGUMENT_PARSER_HPP_
#define _MALVAARGUMENT_PARSER_HPP_

#include <sstream>

#include <getopt.h>

using namespace std;

static const char *USAGE_MESSAGE =
  "Usage: malva-geno [-k KMER-SIZE] [-r REF-KMER-SIZE] [-c MAX-COV] "
  "<reference.fa> <variants.vcf> <kmc_output_prefix>\n"
  "\n"
  "Top notch description of this tool\n"
  "\n"
  "      -h, --help                        display this help and exit\n"
  "      -k, --kmer-size                   size of the kmers to index (default:35)\n"
  "      -r, --ref-kmer-size               size of the reference kmers to index (default:43)\n"
  "      -e, --error-rate                  expected sample error rate (default:0.001)\n"
  "      -s, --samples                     file containing the list of (VCF) samples to consider (default:-, i.e. all samples)\n"
  "      -f, --freq-key                    a priori frequency key in the INFO column of the input VCF (default:AF)\n"
  "      -c, --max-coverage                maximum coverage for variant alleles (default:200)\n"
  "      -b, --bf-size                     bloom filter size in GB (default:4)\n"
  "      -p, --strip-chr                   strip \"chr\" from sequence names (default:false)\n"
  "      -u, --uniform                     use uniform a priori probabilities (default:false)\n"
  "      -v, --verbose                     output COVS and GTS in INFO column (default: false)\n"
  "      -1, --haploid                     run MALVA in haploid mode (default: false)\n"
  // "      -t, --threads                     number of threads (default: 1)\n"
  "\n";

namespace opt {
static uint k = 35;
static uint ref_k = 43;
static float error_rate = 0.001;
static string samples = "-";
static string freq_key = "AF";
static uint max_coverage = 200;
static uint64_t bf_size = ((uint64_t)0b1 << 35);
static bool strip_chr = false;
static bool uniform = false;
static bool verbose = false;
static bool haploid = false;
// static size_t nThreads = 1;
static string fasta_path;
static string vcf_path;
static string kmc_sample_path;
}

static const char *shortopts = "k:r:e:s:f:c:b:hpuv1";

static const struct option longopts[] = {
    {"kmer-size", required_argument, NULL, 'k'},
    {"ref-kmer-size", required_argument, NULL, 'r'},
    {"error-rate", required_argument, NULL, 'e'},
    {"freq-key", required_argument, NULL, 'f'},
    {"samples", required_argument, NULL, 's'},
    {"max-coverage", required_argument, NULL, 'c'},
    {"bf-size", required_argument, NULL, 'b'},
    {"strip-chr", required_argument, NULL, 'p'},
    {"uniform", required_argument, NULL, 'u'},
    {"verbose", no_argument, NULL, 'v'},
    {"haplod", no_argument, NULL, '1'},
    // {"threads", no_argument, NULL, 't'},
    {"help", no_argument, NULL, 'h'},
    {NULL, 0, NULL, 0}};

void parse_arguments(int argc, char **argv) {
  bool die = false;
  for (char c;
       (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'p':
      opt::strip_chr = true;
      break;
    case 'u':
      opt::uniform = true;
      break;
    case 'k':
      arg >> opt::k;
      break;
    case 'r':
      arg >> opt::ref_k;
      break;
    case 'e':
      arg >> opt::error_rate;
      break;
    case 's':
      arg >> opt::samples;
      break;
    case 'f':
      arg >> opt::freq_key;
      break;
    case 'c':
      arg >> opt::max_coverage;
      break;
    case 'b':
      // Let's consider this as GB
      arg >> opt::bf_size;
      opt::bf_size = opt::bf_size * ((uint64_t)0b1 << 33);
      break;
    case 'v':
      opt::verbose = true;
      break;
    case '1':
      opt::haploid = true;
      break;
    // case 't':
    //   arg >> opt::nThreads;
    //   break;
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

  opt::fasta_path = argv[optind++];
  opt::vcf_path = argv[optind++];
  opt::kmc_sample_path = argv[optind++];
}

#endif
