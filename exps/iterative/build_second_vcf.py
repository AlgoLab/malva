import sys

from pysam import VariantFile

def main():
    in_vcf_path = sys.argv[1]
    malva_vcf_path = sys.argv[2]

    idx_to_keep = set()
    for record in VariantFile(malva_vcf_path, 'r'):
        idx = (record.id, record.pos, "-".join(record.alts))
        gt = ""
        for (type_name, content) in record.samples.items()[0][1].items():
            if type_name == 'GT':
                gt = str(min(content[0],content[1])) + "/" + str(max(content[0],content[1]))
        if gt != "0/0":
            idx_to_keep.add(idx)

    in_vcf = VariantFile(in_vcf_path, 'r')
    print(in_vcf.header, end='')
    for record in in_vcf:
        idx = (record.id, record.pos, "-".join(record.alts))
        if idx in idx_to_keep:
            print(record, end='')

if __name__ == '__main__':
    main()
