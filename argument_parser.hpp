#ifndef _MALVAARGUMENT_PARSER_HPP_
#define _MALVAARGUMENT_PARSER_HPP_

#include <sstream>

#include <getopt.h>

static const char *USAGE_MESSAGE =
    "Usage: malva [-k KMER-SIZE] [-r REF-KMER-SIZE] [-c MIN-COV] "
    "<reference.fa> <variants.vcf> <kmc_output_prefix>\n";

namespace opt {
static bool strict_mode = true;
static uint k = 31;
static uint ref_k = 43;
static uint read_len = 150;
static float error_rate = 0.001;
static uint64_t bf_size = ((uint64_t)0b1 << 33);
static uint min_coverage = 2;
static uint sample_coverage = 30;
static std::string pop = "EUR";
// static size_t nThreads = 1;
static std::string fasta_path;
static std::string vcf_path;
static std::string kmc_sample_path;
} // namespace opt

static const char *shortopts = "lk:r:b:c:e:n:s:p:h";

static const struct option longopts[] = {
    {"loose", no_argument, NULL, 'l'},
    {"kmer-size", required_argument, NULL, 'k'},
    {"ref-kmer", required_argument, NULL, 'r'},
    {"bf-size", required_argument, NULL, 'b'},
    {"min-coverage", required_argument, NULL, 'c'},
    {"error-rate", required_argument, NULL, 'e'},
    {"read-length", required_argument, NULL, 'n'},
    {"sample-coverage", required_argument, NULL, 's'},
    {"population", required_argument, NULL, 'p'},
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
    case 'e':
      arg >> opt::error_rate;
      break;
    case 'n':
      arg >> opt::read_len;
      break;
    case 's':
      arg >> opt::sample_coverage;
      break;
    case 'p':
      arg >> opt::pop;
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
