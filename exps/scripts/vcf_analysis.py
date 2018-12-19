import sys, argparse

from pysam import VariantFile

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns; sns.set()

n_elem = 4
idx_map = { '0/0' : 0,
            '0/1' : 1,
            '1/1' : 2,
            '*/*' : 3}

def main():
    parser = argparse.ArgumentParser(description='Script to compare two VCFs based on ref position and genotyped allele.')
    parser.add_argument('truth_path', help='Path to first VCF (used as truth)')
    parser.add_argument('vcf_path', help='Path to second VCF')
    parser.add_argument('-p', dest='print_flag', action='store_true', default=False, help='Set this to print FN/FPs)')
    args = parser.parse_args()

    truth = VariantFile(args.truth_path, 'r')
    output = VariantFile(args.vcf_path, 'r')

    truth_gts = {}
    total = 0
    for record in truth:
        vidx = (record.id, record.pos, "-".join(record.alts))
        gt = ""
        for (type_name, content) in record.samples.items()[0][1].items():
            if type_name == 'GT':
                gt = str(min(content[0],content[1])) + "/" + str(max(content[0],content[1]))
        is_good = True
        for alt in record.alts:
            if alt[0] == '<':
                is_good = False
                break
        if is_good:
            if gt != "0/0":
                total += 1
                truth_gts[vidx] = gt

    table = [[0 for j in range(0,n_elem)] for i in range(0,n_elem)]

    tp = 0
    fp = 0
    out_vcf = open("{}.FP.vcf".format(args.vcf_path[:-4]), 'w')
    for record in output:
        vidx = (record.id, record.pos, "-".join(record.alts))
        gt = ""
        add_flag = False
        truth = ""
        for (type_name, content) in record.samples.items()[0][1].items():
            if type_name == 'GT':
                gt = str(min(content[0],content[1])) + "/" + str(max(content[0],content[1]))
        if vidx in truth_gts:
            if truth_gts[vidx] == gt:
                tp += 1
            else:
                fp += 1
                truth = truth_gts[vidx]
                add_flag = True
        else:
            if gt != "0/0":
                truth = "0/0"
                fp+=1
                add_flag = True
        if add_flag:
            if args.print_flag:
                out_vcf.write("{}\t{}".format(truth, record))
            if truth not in idx_map:
                truth = '*/*'
            if gt not in idx_map:
                gt = '*/*'
            table[idx_map[truth]][idx_map[gt]] += 1

    out_vcf.close()

    fn = total - tp
    P = round(100*tp/(tp+fp),3) if tp+fp != 0 else 0
    R = round(100*tp/(tp+fn),3) if tp+fn != 0 else 0

    out_csv = open("{}.csv".format(args.vcf_path[:-4]), 'w')
    out_csv.write("TP,FP,FN,P,R\n{},{},{},{},{}\n".format(tp, fp, fn, P, R))
    out_csv.close()
    #print("{:>7},{:>7},{:>7},{:>7},{:>7}".format("TP", "FP", "FN", "P", "R"))
    #print("{:>7},{:>7},{:>7},{:>7.3f},{:>7.3f}".format(tp, fp, fn, P, R))

    ax = sns.heatmap(table,
                     annot=True,
                     fmt='d',
                     linewidths=.5,
                     xticklabels=[lab for lab in idx_map],
                     yticklabels=[lab for lab in idx_map])
    ax.invert_yaxis()
    plt.xlabel("given")
    plt.ylabel("truth")
    plt.savefig("{}.pdf".format(args.vcf_path[:-4]))

if __name__ == "__main__":
    main()
