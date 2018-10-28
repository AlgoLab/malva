#!/bin/bash

vcf=$1
sample=$2
out=$3

bcftools view -s ${sample} ${vcf} > ${out}.vcf
grep -v "0|0" ${out}.vcf > ${out}.truth.vcf
