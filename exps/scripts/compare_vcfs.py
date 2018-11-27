import sys, argparse

from pysam import VariantFile

def main():
    parser = argparse.ArgumentParser(description='Script to compare two VCFs based on ref position and genotyped allele.')
    parser.add_argument('--all', dest='all_mode', help='Compare all lines, not only the one with genotype different from 0|0',
                        required=False, action='store_true')
    parser.add_argument('truth', help='Path to first VCF (used as truth)')
    parser.add_argument('vcf', help='Path to second VCF')
    args = parser.parse_args()

    truth = VariantFile(args.truth, 'r')
    output = VariantFile(args.vcf, 'r')

    truth_alts = set()
    total = 0
    for record in truth:
        gt = ""
        for (type_name, content) in record.samples.items()[0][1].items():
            if type_name == 'GT':
                gt = str(content[0]) + "/" + str(content[1])
        for alt in record.alts:
            if alt[0] != '<':
                truth_alts.add((record.pos, record.ref, gt))

    tp = 0
    fp = 0
    for record in output:
        gt = ""
        for (type_name, content) in record.samples.items()[0][1].items():
            if type_name == 'GT':
                gt = str(content[0]) + "/" + str(content[1])
        if (record.pos, record.ref, gt) in truth_alts:
            tp+=1
        else:
            fp+=1
    fn = len(truth_alts) - tp

    P = round(100*tp/(tp+fp),3) if tp+fp != 0 else 0
    R = round(100*tp/(tp+fn),3) if tp+fn != 0 else 0

    print("{:>7},{:>7},{:>7},{:>7},{:>7}".format("TP", "FP", "FN", "P", "R"))
    print("{:>7},{:>7},{:>7},{:>7.3f},{:>7.3f}".format(tp, fp, fn, P, R))

if __name__ == "__main__":
    main()
