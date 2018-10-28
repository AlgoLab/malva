#!/bin/bash

vcf=$1
panel=$2
spop=$3

grep -P "\t${spop}\t" ${panel} | cut -f 1 > sample.${spop}.tmp
bcftools view -S sample.${spop}.tmp ${vcf}
rm sample.${spop}.tmp
