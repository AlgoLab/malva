#!/bin/bash

in_fold=$1
truth=$2

echo "k,K,sample,vcf,TP,FP,FN,P,R,Time(s),RAM(MB)"
for vcf in $(ls ${in_fold}/*.vcf)
do
    k=$(basename ${vcf} | cut -f 1 -d'.' | sed 's/k//g')
    K=$(basename ${vcf} | cut -f 2 -d'.' | sed 's/K//g')
    sample=$(basename ${vcf} | cut -f 3 -d'.')
    vcf_type=$(basename ${vcf} | cut -f 4 -d'.')

    time_f=$(echo $vcf | sed 's/.vcf$/.time/g')
    time=$(grep "wall clock" ${time_f} | cut -f8 -d' ' | awk -F':' '{if (NF == 2) {print $1 * 60 + $2} else {print $1 * 60 * 60 + $2 * 60 + $3}}' | cut -f 1 -d'.')
    RAM=$(( $(grep "Maximum resident set size" ${time_f} | cut -f6 -d' ') / 1000))
    
    echo $k,$K,$sample,$vcf_type,$(python3 compare_vcfs.py $truth $vcf),$time,$RAM
done
