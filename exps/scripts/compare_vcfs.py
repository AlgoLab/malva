import sys

from pysam import VariantFile

def main():
    truth = VariantFile(sys.argv[1], 'r')
    output = VariantFile(sys.argv[2], 'r')

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
            truth_alts.add((record.pos, alt))

    truth_alts = set(truth_alts)

    tp = 0
    fp = 0
    for record in output:
        for alt in record.alts:
            if (record.pos, alt) in truth_alts:
                tp+=1
            else:
                fp+=1
    fn = len(truth_alts) - tp

    # print("{:>7},{:>7},{:>7},{:>7},{:>7}".format("TP", "FP", "FN", "P", "R"))
    # print("{:>7},{:>7},{:>7},{:>7.3f},{:>7.3f}".format(tp, fp, fn, 100*tp/(tp+fp), 100*tp/(tp+fn)))
    print("{},{},{},{},{}".format(tp, fp, fn, round(100*tp/(tp+fp), 3), round(100*tp/(tp+fn), 3)), sep=',')

    return 0

if __name__ == "__main__":
    main()
