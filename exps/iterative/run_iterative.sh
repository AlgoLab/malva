#!/bin/bash

ref=$1
vcf=$2
kmc_fold=$3
truth=$4
out_prefix=$5

./malva -s 30 -n 150 -e 0.001 -k 31 -r 33 -c 2 ${ref} ${vcf} ${kmc_fold}/k33.original > ${out_prefix}.first_iteration.vcf

python3 build_second_vcf.py ${vcf} ${out_prefix}.first_iteration.vcf > ${out_prefix}.first_iteration_no00.vcf  2> /dev/null

./malva -s 30 -n 150 -e 0.001 -k 35 -r 57 -c 2 ${ref} ${out_prefix}.first_iteration_no00.vcf ${kmc_fold}/k57.original > ${out_prefix}.second_iteration.vcf

python3 combine_iterative.py ${out_prefix}.first_iteration.vcf ${out_prefix}.second_iteration.vcf > ${out_prefix}.vcf 2> /dev/null

python3 exps/scripts/vcf_analysis.py -p ${truth} ${out_prefix}.vcf 2> /dev/null
