import sys, argparse

from pysam import VariantFile

# def is_snp(record):
#     return len(record.ref) == 1 and len(record.alts) == 1 and len(record.alts[0]) == 1

def main():
    parser = argparse.ArgumentParser(description='Script to compare two VCFs based on ref position and genotyped allele.')
    # parser.add_argument('-o', dest='only_gt', help='Compare only variants with genotype different from 0|0',
    #                     required=False, action='store_true')
    parser.add_argument('truth', help='Path to first VCF (used as truth)')
    parser.add_argument('vcf', help='Path to second VCF')
    parser.add_argument('-p', dest='print', action='store_true', default=False, help='Set this to print FN/FPs)')
    parser.add_argument('-c', dest='check', default='@', help='Path to third VCF (used to check second FN/FPs)')
    args = parser.parse_args()

    truth = VariantFile(args.truth, 'r')
    output = VariantFile(args.vcf, 'r')

    truth_gts = {}
    for record in truth:
        # if not is_snp_for_vargeno(record):
        #     continue
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
            truth_gts[(record.id, "-".join(record.alts))] = gt

    tp = 0
    fp = 0
    FPs = {}
    for record in output:
        # if not is_snp(record):
        #     continue
        gt = ""
        for (type_name, content) in record.samples.items()[0][1].items():
            if type_name == 'GT':
                gt = str(min(content[0],content[1])) + "/" + str(max(content[0],content[1]))
        if truth_gts[(record.id, "-".join(record.alts))] == gt:
            tp+=1
        else:
            FPs[(record.id, "-".join(record.alts))] = gt
            fp+=1
    fn = len(truth_gts) - tp

    P = round(100*tp/(tp+fp),3) if tp+fp != 0 else 0
    R = round(100*tp/(tp+fn),3) if tp+fn != 0 else 0

    print("{:>7},{:>7},{:>7},{:>7},{:>7}".format("TP", "FP", "FN", "P", "R"), file=sys.stderr)
    print("{:>7},{:>7},{:>7},{:>7.3f},{:>7.3f}".format(tp, fp, fn, P, R), file=sys.stderr)

    if args.print:
        truth = VariantFile(args.truth, 'r')
        print(truth.header, end="")
        for record in truth:
            gt = ""
            for (type_name, content) in record.samples.items()[0][1].items():
                if type_name == 'GT':
                    gt = str(min(content[0],content[1])) + "/" + str(max(content[0],content[1]))
            if (record.id, "-".join(record.alts)) in FPs:
                print(record, end="")

    if args.check != '@':
        check = VariantFile(args.check, 'r')
        check_gts = {}
        for record in check:
            gt = ""
            for (type_name, content) in record.samples.items()[0][1].items():
                if type_name == 'GT':
                    gt = str(min(content[0],content[1])) + "/" + str(max(content[0],content[1]))
            check_gts[(record.id, "-".join(record.alts))] = gt
        total_fps = 0
        not_in_third = 0
        correct_in_third = 0
        same_in_third = 0
        different_in_third = 0
        for (key,gt) in FPs.items():
            total_fps += 1
            if key in check_gts:
                if check_gts[key] == truth_gts[key]:
                    correct_in_third += 1
                else:
                    if check_gts[key] == gt:
                        same_in_third += 1
                    else:
                        different_in_third += 1
            else:
                not_in_third += 1

        print("", file=sys.stderr)
        print("{:>7},{:>7},{:>7},{:>7},{:>7}".format("FPs", "no3rd", "correct3rd", "same3rd", "diff3rd"), file=sys.stderr)
        print("{:>7},{:>7},{:>7},{:>7},{:>7}".format(total_fps, not_in_third, correct_in_third, same_in_third, different_in_third), file=sys.stderr)

if __name__ == "__main__":
    main()
