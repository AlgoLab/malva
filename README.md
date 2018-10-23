# MALVA: Mapping-free ALternate-allele detection of known VAriants

#### !!! This is the development branch !!!

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

```
./malva [-k ksize] [-m min_coverage] [-x max_coverage] <reference> <variants> <kmc_output_prefix>
```

## Authors

* [Marco Previtali](https://algolab.eu/people/previtali/)
* [Luca Denti](https://algolab.eu/people/luca-denti/)
* [Giulia Bernardini](https://algolab.eu/people/giulia-bernardini)
* [Paola Bonizzoni](https://algolab.eu/people/bonizzoni/)

For inquiries on this software please contact either MP or LD.
