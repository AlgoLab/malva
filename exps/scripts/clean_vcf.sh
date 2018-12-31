#!/bin/bash

vcf=$1
#out_vcf="$(dirname ${vcf})/$(basename ${vcf} .vcf).clean.vcf"

grep "^#" ${vcf} #> ${out_vcf}
grep -v "^#" ${vcf} | awk '{$8="."; print $0}' | sed 's/ /\t/g' #>> ${out_vcf}
