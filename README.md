# MALVA: Mapping-free ALternate-allele detection of known VAriants

Detect known alterante alleles directly from a sample of reads.

## Dependencies

MALVA requires the following libraries and tools:

* [sdsl-lite v2.1.1](https://github.com/simongog/sdsl-lite/tree/v2.1.1)
* [KMC v3](https://github.com/refresh-bio/KMC/tree/v3.0.0)
* [htslib v1.8](https://github.com/samtools/htslib/tree/1.8)

This repository comes with them as submodules so you don't need to clone them separately.

## Download and Installation

To download and compile the code run the following commands.

First clone the repository and cd into it.

```bash
git clone --recursive https://github.com/AlgoLab/malva.git
cd malva
```

If you have KMC3, sdsl-lite, and htslib already installed you can skip the following commands.

**Please Note**: KMC seems to not compile using `g++-7` under Ubuntu; please switch to `g++-6` in this step to avoid compilation errors.

```bash
cd sdsl-lite/build
./build.sh
cd ../../KMC
make
cd ../htslib
make
cd ..
```

You can now compile `malva` from the root of you local copy of the repository simply by running make.

```bash
cd <path-to-malva-local-repo>
make
```

## Usage

You can run the tool using the `MALVA` bash script as follows.

```
Usage: MALVA [-k ksize] [-m min_coverage] [-x max_coverage] [-l] <reference> <variants> <sample>
Preprocesses the VCF file and calls malva and outputs a VCF containing the
known variants in <sample> to STDOUT.

Ideally this script will be depracated soon.

Arguments:
     -h              print this help and exit
     -k              kmers length (default:15)
     -r              reference kmers length (default:31)
     -m              minimum coverage for filtering variants (default:30)
     -x              maximum coverage for filtering variants (default: unlimited)
     -l              loose mode (default: false)
     -t              number of threads used by KMC (default: 1)
     -u              maximum (even) number of neighbouring variants considered
                     for each signature (default: 4)

Positional arguments:
     <reference>     reference file in FASTA/FASTQ format
     <variants>      variants file in VCF format
     <sample>        sample file in FASTA/FASTQ format
```

If you already counted the kmers in your sample using KMC3 you can directly use `malva` as follows.

```
Usage: malva [-k KMER-SIZE] [-c MIN-COV] [-x MAX-COV] <reference.fa> <variants.vcf> <kmc_output_prefix>

Top notch description of this tool

      -h, --help                        display this help and exit
      -l, --loose                       loose mode. If set then the reference kmers will not be removed (default: false)
      -t, --threads                     number of threads (default: 1)

      -k, --kmer-size                   size of the kmers to index (default:15)
      -c, --min-coverage                minimum coverage for validated variants (default: 30)
      -x, --max-coverage                maximum coverage for validated variants.
                                        If this value is greter than 250 then no upper bound is set (default: unlimited)
      -u, --max-combinations            maximum (even) number of neighbouring variants considered
                                        for each signature (default: 4)
```

## Authors

* [Marco Previtali](https://algolab.eu/people/previtali/)
* [Luca Denti](https://algolab.eu/people/luca-denti/)
* [Giulia Bernardini](https://algolab.eu/people/giulia-bernardini)
* [Paola Bonizzoni](https://algolab.eu/people/bonizzoni/)

For inquiries on this software please contact either MP or LD.
