# MALVA: genotyping by Mapping-free ALternate-allele detection of known VAriants

Alignment-free genotyping of a set of known variants (in VCF format) directly from a sample of reads.

## Dependencies

MALVA requires the following libraries and tools:

* [sdsl-lite v2.1.1](https://github.com/simongog/sdsl-lite/tree/v2.1.1)
* [KMC v3.1.0](https://github.com/refresh-bio/KMC/tree/v3.1.0)
* [htslib v1.9](https://github.com/samtools/htslib/tree/1.9)

This repository comes with them as submodules so you don't need to clone them separately.

## Download and Installation

To download and compile the code run the following commands.

First clone the repository and cd into it.

```bash
git clone --recursive https://github.com/AlgoLab/malva.git
cd malva
```

If you have KMC3, sdsl-lite, and htslib already installed you can skip the following commands.

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

```
./malva <reference> <variants> <kmc_output_prefix>
```

## Authors

* [Marco Previtali](https://algolab.eu/people/previtali/)
* [Luca Denti](https://algolab.eu/people/luca-denti/)
* [Giulia Bernardini](https://algolab.eu/people/giulia-bernardini)
* [Paola Bonizzoni](https://algolab.eu/people/bonizzoni/)

For inquiries on this software please contact either MP or LD.
