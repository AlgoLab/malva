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
        alts = []
        for (type_name, content) in record.samples.items()[0][1].items():
            if type_name == 'GT':
                ap1 = int(content[0])
                if ap1 != 0:
                    alts.append(record.alts[ap1-1])
                ap2 = int(content[1])
                if ap2 != 0:
                    alts.append(record.alts[ap2-1])
        for alt in alts:
            if alt[0] != '<':
                truth_alts.add((record.pos, alt))
        if args.all_mode and len(alts) == 0:
            truth_alts.add((record.pos, '@'))

    truth_alts = set(truth_alts)

    tp = 0
    fp = 0
    for record in output:
        alts = []
        for (type_name, content) in record.samples.items()[0][1].items():
            if type_name == 'GT':
                if content == (0,0):
                    if args.all_mode:
                        alts = ["@"]
                else:
                    alts = record.alts
        for alt in alts:
            if (record.pos, alt) in truth_alts:
                tp+=1
            else:
                fp+=1
    fn = len(truth_alts) - tp

    # print("{:>7},{:>7},{:>7},{:>7},{:>7}".format("TP", "FP", "FN", "P", "R"))
    # print("{:>7},{:>7},{:>7},{:>7.3f},{:>7.3f}".format(tp, fp, fn, 100*tp/(tp+fp), 100*tp/(tp+fn)))
    P = 100*tp/(tp+fp) if tp+fp != 0 else 0
    R = 100*tp/(tp+fn) if tp+fn != 0 else 0
    print("{},{},{},{},{}".format(tp, fp, fn, round(P, 3), round(R, 3)), sep=',')

    return 0

if __name__ == "__main__":
    main()
