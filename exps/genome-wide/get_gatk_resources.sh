#!/bin/bash

function get_vcf {
  fname=$1
  base_url=ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19
  wget ${base_url}/${fname}.gz
  gunzip ${fname}.gz
  sed -i 's/=chr/=/g' ${fname}
  sed -i 's/^chr//g' ${fname}
}

for fname in dbsnp_138.hg19.vcf hapmap_3.3.hg19.sites.vcf 1000G_omni2.5.hg19.sites.vcf 1000G_phase1.snps.high_confidence.hg19.sites.vcf Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
do
  get_vcf ${fname}
done
