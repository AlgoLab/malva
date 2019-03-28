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
    -p, --population                  population to consider while reading input VCF (default:EUR)
    -c, --max-coverage                maximum coverage for variant alleles (default:200)
    -b, --bf-size                     bloom filter size in GB (default:4)

Positional arguments:
    <reference>                       reference file in FASTA format (may be gzipped)
    <variants>                        variants file in VCF format (may be gzipped)
    <kmc_output_prefix>               prefix of KMC output
```

The file needed by malva whose prefix is `<kmc_output_prefix>` can be computed with KMC as follows:
```
cd <path-to-malva-local-repo>
./KMC/bin/KMC -k <REF-KMER-SIZE> <sample> <kmc_output_prefix>
```

Anyway, we provide a bash script that you can use to run the full pipeline `KMC+malva-geno`:
```
Usage: MALVA [-k KMER-SIZE] [-r REF-KMER-SIZE] [-c MAX-COV] <reference> <variants> <sample>

Arguments:
     -h              print this help and exit
     -k              size of the kmers to index (default:35)
     -r              size of the reference kmers to index (default:43)
     -e              expected sample error rate (default:0.001)
     -p              population to consider while reading input VCF (default:EUR)
     -c              maximum coverage for variant alleles (default:200)
     -b              bloom filter size in GB (default:4)
     -m              max amount of RAM in GB - KMC parameter (default:4)

Positional arguments:
    <reference>     reference file in FASTA format (can be gzipped)
    <variants>      variants file in VCF format (can be gzipped)
    <sample>        sample file in FASTA/FASTQ format (can be gzipped)
```

## Example
After you compile `malva`, you can test it on the example data provided (note that we set the Bloom filter size to 1GB):
```
cd example
tar xvfz data.tar.gz
../MALVA chr20.fa chr20.vcf chr20.sample.fa -b 1 > chr20.genotyped.vcf
```

This should take less than 1 minute to complete. You can also verify
the correcteness of the output VCF `chr20.genotyped.vcf` by comparing
it with [chr20.malva.vcf](https://github.com/AlgoLab/malva/blob/master/example/chr20.malva.vcf).

### Note
- The tool has been tested only on 64bit Linux system.
- The current release is optimized for working with the VCF files provided by the 1000GenomesProject (it uses the value encoded as `*_AF` in the `INFO` field to compute the _a priori_ frequency of each allele). If your VCF file uses a different format, let us know

## Authors

* [Marco Previtali](https://algolab.eu/people/previtali/)
* [Luca Denti](https://algolab.eu/people/luca-denti/)
* [Giulia Bernardini](https://algolab.eu/people/giulia-bernardini)
* [Paola Bonizzoni](https://algolab.eu/people/bonizzoni/)
* [Alexander Schönhuth](https://homepages.cwi.nl/~as/)

For inquiries on this software please contact either MP or LD.

## License
MALVA is distributed under the GPL-3.0-or-later license.
