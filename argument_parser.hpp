#ifndef _MALVAARGUMENT_PARSER_HPP_
#define _MALVAARGUMENT_PARSER_HPP_

#include <sstream>

#include <getopt.h>

static const char *USAGE_MESSAGE =
  "Usage: malva [-k KMER-SIZE] <reference.fa> <variants.vcf>";

namespace opt {
  static uint k = 31;
  static uint bf_size = ((uint)0b1 << 31);
  // static size_t minCoverage = 30;
  // static size_t maxCoverage = 251;
  // static size_t nThreads = 1;
  static std::string fasta_path;
  static std::string vcf_path;
  // static string sample_path;
}

static const char *shortopts = "k:r:u:c:x:b:t:lh";

static const struct option longopts[] = {
  {"kmer-size", required_argument, NULL, 'k'},
  // {"min-coverage", required_argument, NULL, 'c'},
  // {"max-coverage", required_argument, NULL, 'x'},
  // {"bf-size", required_argument, NULL, 'b'},
  // {"threads", no_argument, NULL, 't'},
  {"help", no_argument, NULL, 'h'},
  {NULL, 0, NULL, 0}};

void parse_arguments(int argc, char **argv) {
  bool die = false;
  for (char c;
       (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'b':
      // Let's consider this as MB
      arg >> opt::bf_size;
      opt::bf_size = opt::bf_size * (0b1 << 20);
      break;
      // case 'c':
      //   arg >> opt::minCoverage;
      //   break;
      // case 'x':
      //   arg >> opt::maxCoverage;
      //   break;
    case 'k':
      arg >> opt::k;
      break;
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
  /**
     if (argc - optind < 3) {
     cerr << "malva : missing arguments\n";
     die = true;
     } else if (argc - optind > 3) {
     cerr << "malva : too many arguments\n";
     die = true;
     }
  **/
  if (die) {
    std::cerr << "\n" << USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }

  opt::fasta_path = argv[optind++];
  opt::vcf_path = argv[optind++];
  // opt::sequencePath = argv[optind++];
}

#endif
