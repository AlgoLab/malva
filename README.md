[![Anaconda-Server Badge](https://anaconda.org/bioconda/malva/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/malva/badges/platforms.svg)](https://anaconda.org/bioconda/malva)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/malva/badges/license.svg)](https://anaconda.org/bioconda/malva)

# MALVA: genotyping by Mapping-free ALternate-allele detection of known VAriants

Alignment-free genotyping of a set of known variants (in VCF format) directly from a sample of reads.

## Install

`MALVA` is available on bioconda.
```shell
$ conda create -n malvatest -c bioconda malva
```
will create an environment named `malvatest` that includes `MALVA` and its dependencies.

## Install from source code

### Dependencies

To manually compile MALVA you'll need a C++17-compliant compiler, [CMake](https://cmake.org), and the following libraries installed in your system:
* [sdsl-lite v2.1.1](https://github.com/simongog/sdsl-lite) under [GPLv3 License](https://github.com/simongog/sdsl-lite/blob/master/COPYING)
* [KMC >= v2.3](https://github.com/refresh-bio/KMC) under GPLv3 License
* [htslib >= v1.10.2](https://github.com/samtools/htslib) under [MIT/Expat License](https://github.com/samtools/htslib/blob/develop/LICENSE)
* [zstd >= 1.5.1](https://github.com/facebook/zstd) under [BSD License](https://github.com/facebook/zstd/blob/dev/LICENSE) and [GPLv2](https://github.com/facebook/zstd/blob/dev/COPYING)
* [zstdstream](https://github.com/TorchCraft/TorchCraftAI/blob/main/common/zstdstream.cpp) (included in this repository) under [MIT License](https://github.com/TorchCraft/TorchCraftAI/blob/main/LICENSE)
* [xxHash](https://github.com/Cyan4973/xxHash) (included in this repository) under BSD 2-Clause License
* [zlib](https://www.zlib.net) under [zlib license](https://www.zlib.net/zlib_license.html)

Use your favorite system-wide package manager to install them before compiling MALVA.

Alternatively, you can also use conda (and bioconda) to install dependencies. For example (please adapt to you system setup):

``` shell
conda create -n malvadeps -c conda-forge -c bioconda htslib kmc sdsl-lite cmake zstd-static cxx-compiler
```

Notice that these dependencies are needed only if you want to compile MALVA from sources,
since otherwise it is already available on Bioconda in binary form (see [above](#install)).


### Download and installation

To download and compile the code run the following commands.

```shell
git clone https://github.com/AlgoLab/malva.git
cd malva
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

If the compilation is successful, the `malva-geno` binary will be copied to the `${PROJECT_ROOT}/bin` directory.

## Usage
```
Usage: malva-geno <subcommand> [-k KMER-SIZE] [-r REF-KMER-SIZE] [-c MAX-COV] <reference> <variants> <kmc_output_prefix>

Arguments:
	<subcommand>                      either index to create the reference index or call to call call the genotypes.

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
kmc -k<REF-KMER-SIZE> <sample> <kmc_output_prefix> <kmc_tmp_dir>
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
```shell
cd example
tar xvfz data.tar.gz
../MALVA -k 35 -r 43 -b 1 -f EUR_AF chr20.fa chr20.vcf chr20.sample.fa > chr20.genotyped.vcf
```

The last command is equivalent to run:
```shell
mkdir -p kmc_tmp
kmc -m4 -k43 -fm chr20.sample.fa kmc.out kmc_tmp
../bin/malva-geno index -k 35 -r 43 -b 1 -f EUR_AF chr20.fa chr20.vcf kmc.out
../bin/malva-geno call -k 35 -r 43 -b 1 -f EUR_AF chr20.fa chr20.vcf kmc.out > chr20.genotyped.vcf
```

This should take less than 1 minute to complete. You can also verify the correcteness of the output VCF `chr20.genotyped.vcf` by comparing it with `chr20.malva.vcf`.
For a quick comparison look [MALVA-TEST](https://github.com/Bunco3/malva/tree/path_to_version_2_0/malva_test): Workflow for testing MALVA Output.

## Haploid mode - Example
To run MALVA in haploid mode just use the `-1` argument.
```
cd example
tar xvfz haploid.tar.gz
../MALVA -1 -k 35 -r 43 -b 1 -f AF haploid.fa haploid.vcf haploid.fq > haploid.genotyped.vcf
```
This should take less than 1 minute to complete. You can also verify the correcteness of the output VCF `haploid.genotyped.vcf` by comparing it with `haploid.malva.vcf`.
For a quick comparison look [MALVA-TEST](https://github.com/Bunco3/malva/tree/path_to_version_2_0/malva_test): Workflow for testing MALVA Output.

### Note
- The tool has been tested only on 64bit Linux system.

## Authors

* [Marco Previtali](https://algolab.eu/people/previtali/)
* [Luca Denti](https://algolab.eu/people/luca-denti/)
* [Giulia Bernardini](https://algolab.eu/people/giulia-bernardini)
* [Paola Bonizzoni](https://algolab.eu/people/bonizzoni/)
* [Alexander Schönhuth](https://homepages.cwi.nl/~as/)
* Marco Burgio (code refactoring and `malva-test` utility)

For inquiries on this software please contact either MP or LD.

## License
MALVA is distributed under the GPL-3.0-or-later license.
