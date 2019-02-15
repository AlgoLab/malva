---
---
# MALVA: genotyping by Mapping-free ALternate-allele detection of known VAriants

MALVA is an alignment-free genotyping tool.
Given a set of reads in FASTA/Q format and a set of variants in VCF format, it outputs a VCF file including genotype calls for each variant.

## Dependencies

MALVA requires the following libraries and tools:

* [sdsl-lite v2.1.1](https://github.com/simongog/sdsl-lite/tree/v2.1.1)
* [KMC v3.1.0](https://github.com/refresh-bio/KMC/tree/v3.1.0)
* [htslib v1.9](https://github.com/samtools/htslib/tree/1.9)

The libraries are provided as submodules of this repository.

## Download and Installation

First clone the repository using the `--recursive` flag.

```bash
git clone --recursive https://github.com/AlgoLab/malva.git
```

Compile the required libraries if needed.

```bash
cd malva/sdsl-lite/build
./build.sh
cd ../../KMC
make
cd ../htslib
make
```

Compile MALVA using make in the root directory of the git repository.

```bash
cd <path-to-malva-local-repo>
make
```

## Usage

Use KMC3 to count the k-mers in your dataset and provide its output to MALVA.

```
Usage: malva [-k KMER-SIZE] [-r REF-KMER-SIZE] [-c MAX-COV] <reference.fa> <variants.vcf> <kmc_output_prefix>

Parameters

      -h, --help                        display this help and exit
      -l, --loose                       loose mode. If set then the reference kmers will not be removed (default: false)
      -k, --kmer-size                   size of the kmers to index (default:35)
      -r, --ref-kmer-size               size of the reference kmers to index (default:43)
      -n, --read-len                    length of input reads (default:150)
      -e, --error-rate                  expected sample error rate (default:0.001)
      -p, --population                  population to consider while reading input VCF (default:EUR)
      -c, --max-coverage                maximum coverage for variant alleles (default:200)
      -b, --bf-size                     bloom filter size in GB (default:4)
      -v, --verbose                     output a detailed VCF (more information in the INFO column)
```

## Authors

* [Marco Previtali](https://algolab.eu/people/previtali/)
* [Luca Denti](https://algolab.eu/people/luca-denti/)
* [Giulia Bernardini](https://algolab.eu/people/giulia-bernardini)
* [Paola Bonizzoni](https://algolab.eu/people/bonizzoni/)
* [Alexander Schönhuth](https://homepages.cwi.nl/~as/)

For inquiries on this software please contact either MP or LD.
