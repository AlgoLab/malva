Let EXP_FOLD be the folder where you want to run the experiments.

1. Download data:
```
cd ${EXP_FOLD}

# Reference
wget ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz

# VCFs
bash ${MALVA_LOCAL_REPO}/exps/genome-wide/get_vcf.sh

# Sample
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam.bai

# GATK resources
bash ${MALVA_LOCAL_REPO}/exps/genome-wide/get_gatk_resources.sh
```
