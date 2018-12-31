#!/bin/bash

ref=$1
truth=$2
vcf=$3
out=$4

export HGREF=$ref
/home/prj_vcfgeno/hap.py-install/bin/hap.py -r ${ref} ${truth} ${vcf} -o ${out}
