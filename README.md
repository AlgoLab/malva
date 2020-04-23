[![Anaconda-Server Badge](https://anaconda.org/bioconda/malva/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/malva/badges/platforms.svg)](https://anaconda.org/bioconda/malva)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/malva/badges/license.svg)](https://anaconda.org/bioconda/malva)

# MALVA: genotyping by Mapping-free ALternate-allele detection of known VAriants

Alignment-free genotyping of a set of known variants (in VCF format) directly from a sample of reads.

## Install

`MALVA` is available on bioconda.
```bash
conda create -n malvatest malva
```
will create an environment named `malvatest` that includes `MALVA` and its dependencies.

## Install from source code

### Dependencies

MALVA requires the following libraries and tools:

* [sdsl-lite v2.1.1](https://github.com/simongog/sdsl-lite/tree/v2.1.1)
* [KMC >= v3.1.1](https://github.com/refresh-bio/KMC/tree/v3.1.1)
* [htslib >= v1.10.2](https://github.com/samtools/htslib/tree/1.10.2)

This repository comes with them as submodules so you don't need to clone them separately.

### Download and installation

To download and compile the code run the following commands.

First clone the repository and cd into it.

```bash
git clone --recursive https://github.com/AlgoLab/malva.git
cd malva
```

If you have KMC3, sdsl-lite, and htslib already installed you can skip the following commands.

```bash
cd sdsl-lite
./install.sh ..
cd ../KMC
make
cd ../htslib
make
cd ..
```

You can now compile MALVA from the root of you local copy of the repository simply by running make.

```bash
cd <path-to-malva-local-repo>
make
```

## Usage
```
Usage: malva-geno [-k KMER-SIZE] [-r REF-KMER-SIZE] [-c MAX-COV] <reference> <variants> <kmc_output_prefix>

Arguments:
    -h, --help                        display this help and exit
    -k, --kmer-size                   size of the kmers to index (default:35)
    -r, --ref-kmer-size               size of the reference kmers to index (default:43)
    -e, --error-rate                  expected sample error rate (default:0.001)
    -s, --samples                     file containing the list of (VCF) samples to consider (default:-, i.e. all samples)
    -f, --freq-key                    a priori frequency key in the INFO column of the input VCF (default:AF)
    -c, --max-coverage                maximum coverage for variant alleles (default:200)
    -b, --bf-size                     bloom filter size in GB (default:4)
    -p, --strip-chr                   strip "chr" from sequence names (default:false)
    -u, --uniform                     use uniform a priori probabilities (default:false)
    -v, --verbose                     output COVS and GTS in INFO column (default: false)
    -1, --haploid                     run MALVA in haploid mode (default: false)

Positional arguments:
    <reference>                       reference file in FASTA format (may be gzipped)
    <variants>                        variants file in VCF format (may be gzipped)
    <kmc_output_prefix>               prefix of KMC output
```

The file needed by malva whose prefix is `<kmc_output_prefix>` can be computed with KMC as follows:
```
cd <path-to-malva-local-repo>
./KMC/bin/KMC -k<REF-KMER-SIZE> <sample> <kmc_output_prefix> <kmc_tmp_dir>
```

Anyway, we provide a bash script that you can use to run the full pipeline `KMC+malva-geno`:
```
Usage: MALVA [-k KMER-SIZE] [-r REF-KMER-SIZE] [-c MAX-COV] <reference> <variants> <sample>

Arguments:
     -h              print this help and exit
     -k              size of the kmers to index (default:35)
     -r              size of the reference kmers to index (default:43)
     -e              expected sample error rate (default:0.001)
     -s              file containing the list of (VCF) samples to consider (default:-, i.e. all samples)
     -f              a priori frequency key in the INFO column of the input VCF (default:AF)
     -c              maximum coverage for variant alleles (default:200)
     -b              bloom filter size in GB (default:4)
     -m              max amount of RAM in GB - KMC parameter (default:4)
     -p              strip "chr" from sequence names (dafault:false)
     -u              use uniform a priori probabilities (default:false)
     -v              output COVS and GTS in INFO column (default: false)
     -1              run MALVA in haploid mode (default: false)

Positional arguments:
    <reference>     reference file in FASTA format (can be gzipped)
    <variants>      variants file in VCF format (can be gzipped)
    <sample>        sample file in FASTA/FASTQ format (can be gzipped)
```

## Example
After you compiled `malva`, you can test it on the example data provided:
```
cd example
tar xvfz data.tar.gz
../MALVA -k 35 -r 43 -b 1 -f EUR_AF chr20.fa chr20.vcf chr20.sample.fa > chr20.genotyped.vcf
```

The last command is equivalent to run:
```
mkdir -p kmc_tmp
../KMC/bin/kmc -m4 -k43 -fm chr20.sample.fa kmc.out kmc_tmp
../malva-geno -k 35 -r 43 -b 1 -f EUR_AF chr20.fa chr20.vcf kmc.out > chr20.genotyped.vcf
```

This should take less than 1 minute to complete. You can also verify
the correcteness of the output VCF `chr20.genotyped.vcf` by comparing
it with [chr20.malva.vcf](https://github.com/AlgoLab/malva/blob/master/example/chr20.malva.vcf).

## Haploid mode - Example
To run MALVA in haploid mode just use the `-1` argument.
```
cd example
tar xvfz haploid.tar.gz
../MALVA -1 -k 35 -r 43 -b 1 -f AF haploid.fa haploid.vcf haploid.fq > haploid.genotyped.vcf
```

### Note
- The tool has been tested only on 64bit Linux system.

## Authors

* [Marco Previtali](https://algolab.eu/people/previtali/)
* [Luca Denti](https://algolab.eu/people/luca-denti/)
* [Giulia Bernardini](https://algolab.eu/people/giulia-bernardini)
* [Paola Bonizzoni](https://algolab.eu/people/bonizzoni/)
* [Alexander Schönhuth](https://homepages.cwi.nl/~as/)

For inquiries on this software please contact either MP or LD.

## License
MALVA is distributed under the GPL-3.0-or-later license.
