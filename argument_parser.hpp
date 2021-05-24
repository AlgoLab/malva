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

#include <iostream>
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

struct OPT {
 uint k = 35;
 uint ref_k = 43;
 float error_rate = 0.001;
 string samples = "-";
 string freq_key = "AF";
 uint max_coverage = 200;
 uint64_t bf_size = ((uint64_t)0b1 << 35);
 bool strip_chr = false;
 bool uniform = false;
 bool verbose = false;
 bool haploid = false;
//  size_t nThreads = 1;
 string fasta_path;
 string vcf_path;
 string kmc_sample_path;
};

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

void parse_arguments(int argc, char **argv, struct OPT& opt);

#endif
