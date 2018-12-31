#!/bin/bash

fold=$1

echo "k,K,P_SNP,R_SNP,P_INDEL,R_INDEL,Time_sec,RAM_GB"
for csv in $(ls ${fold}/happy/*.summary.csv)
do
    fname=$(basename ${csv} .summary.csv)
    k=$(echo ${fname} | cut -f 1 -d'.')
    k=${k:1:2}
    K=$(echo ${fname} | cut -f 2 -d'.')
    K=${K:1:2}
    echo -n $k,$K,
    log=${fold}/malva/malva.${fname}.time
    for vartype in "SNP" "INDEL"
    do
        P=$(grep "${vartype},ALL" ${csv} | cut -f 12 -d',')
        P=${P:0:5}
        R=$(grep "${vartype},ALL" ${csv} | cut -f 11 -d',')
        R=${R:0:5}
        echo -n $P,$R,
    done

    RAM=$(grep "Maximum" ${log} | cut -f 6 -d' ')
    RAM=${RAM:0:2}
    time=$(grep "wall" ${log} | cut -f8 -d' ' | awk -F':' '{if (NF == 2) {print $1 * 60 + $2} else {print $1 * 60 * 60 + $2 * 60 + $3}}' | cut -f 1 -d'.')
    echo ${time},${RAM}
done
