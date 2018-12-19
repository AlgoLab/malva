#!/bin/bash

base_url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
panel="integrated_call_male_samples_v3.20130502.ALL.panel"
chr1_vcf="ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
vcfs=("ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz")
pop="EUR"
out_vcf="1kgen.ph3.EUR.vcf"

function fetch {
    wget ${base_url}/$1
}

fetch ${panel}
tmpf="samples.tmp"
grep -P "\t${pop}\t" ${panel} | cut -f 1 > ${tmpf}

echo ${chr1_vcf}
fetch ${chr1_vcf}
fetch ${chr1_vcf}.tbi
bcftools view -S ${tmpf} ${chr1_vcf} > ${out_vcf}
rm ${chr1_vcf}
rm ${chr1_vcf}.tbi

for vcf in ${vcfs[*]}
do
    echo ${vcf}
    fetch ${vcf}
    fetch ${vcf}.tbi
    bcftools view -S ${tmpf} ${in_vcf} >> ${out_vcf}
    rm ${vcf}
    rm ${vcf}.tbi
done
rm ${tmpf}
