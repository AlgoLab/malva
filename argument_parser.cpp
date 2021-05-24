//
// Created by Marco on 24/05/2021.
//

#include "argument_parser.hpp"

using namespace std;

void parse_arguments(int argc, char **argv, struct OPT& opt) {
    bool die = false;
    for (char c;
         (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'p':
                opt.strip_chr = true;
                break;
            case 'u':
                opt.uniform = true;
                break;
            case 'k':
                arg >> opt.k;
                break;
            case 'r':
                arg >> opt.ref_k;
                break;
            case 'e':
                arg >> opt.error_rate;
                break;
            case 's':
                arg >> opt.samples;
                break;
            case 'f':
                arg >> opt.freq_key;
                break;
            case 'c':
                arg >> opt.max_coverage;
                break;
            case 'b':
                // Let's consider this as GB
                arg >> opt.bf_size;
                opt.bf_size = opt.bf_size * ((uint64_t)0b1 << 33);
                break;
            case 'v':
                opt.verbose = true;
                break;
            case '1':
                opt.haploid = true;
                break;
                // case 't':
                //   arg >> opt.nThreads;
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

    opt.fasta_path = argv[optind++];
    opt.vcf_path = argv[optind++];
    opt.kmc_sample_path = argv[optind++];
}