#ifndef _MALVAARGUMENT_PARSER_HPP_
#define _MALVAARGUMENT_PARSER_HPP_

#include <sstream>

#include <getopt.h>

static const char *USAGE_MESSAGE =
  "Usage: malva [-a] [-k KMER-SIZE] [-r REF-KMER-SIZE] [-c MIN-COV] <reference.fa> <variants.vcf> <kmc_output_prefix>\n";

namespace opt {
  static bool strict_mode = true;
  static bool all_variants = false;
  static uint k = 31;
  static uint ref_k = 43;
  static uint bf_size = ((uint)0b1 << 31);
  static uint min_coverage = 30;
  // static uint max_coverage = 251;
  // static size_t nThreads = 1;
  static std::string fasta_path;
  static std::string vcf_path;
  static std::string kmc_sample_path;
}

static const char *shortopts = "lak:r:b:c:h";

static const struct option longopts[] = {
  {"loose", no_argument, NULL, 'l'},
  {"all-variant", no_argument, NULL, 'a'},
  {"kmer-size", required_argument, NULL, 'k'},
  {"ref-kmer", required_argument, NULL, 'r'},
  {"bf-size", required_argument, NULL, 'b'},
  {"min-coverage", required_argument, NULL, 'c'},
  // {"max-coverage", required_argument, NULL, 'x'},
  // {"threads", no_argument, NULL, 't'},
  {"help", no_argument, NULL, 'h'},
  {NULL, 0, NULL, 0}};

void parse_arguments(int argc, char **argv) {
  bool die = false;
  for (char c;
       (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'l':
      opt::strict_mode = false;
      break;
    case 'a':
      opt::all_variants = true;
      break;
    case 'k':
      arg >> opt::k;
      break;
    case 'r':
      arg >> opt::ref_k;
      break;
    case 'b':
      // Let's consider this as MB
      arg >> opt::bf_size;
      opt::bf_size = opt::bf_size * (0b1 << 20);
      break;
    case 'c':
      arg >> opt::min_coverage;
      break;
    // case 'x':
    //   arg >> opt::max_coverage;
    //   break;
    // case 't':
    //   arg >> opt::nThreads;
    //   break;
    case '?':
      die = true;
      break;
    case 'h':
      std::cout << USAGE_MESSAGE;
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
    std::cerr << "\n" << USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }

  opt::fasta_path = argv[optind++];
  opt::vcf_path = argv[optind++];
  opt::kmc_sample_path = argv[optind++];
}

#endif
