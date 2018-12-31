import sys

from pysam import VariantFile

def main():
    vcf_path = sys.argv[1]
    vcf = VariantFile(vcf_path, 'r', drop_samples = True)

    vcf.header.add_line("##INFO=<ID=CAF,Number=.,Type=String,Description=\"An ordered, comma delimited list of allele frequencies, starting with the reference allele followed by alternate alleles as ordered in the ALT column.\">")

    for record in vcf.header.records:
        if record.key == "FORMAT":
            record.remove()

    print('\n'.join(str(vcf.header).split('\n')[:-2]))
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

    for record in vcf:
        afs = record.info.get('EUR_AF')
        rf = 1.0 - sum(afs)
        record.info.__setitem__('CAF', '{},{}'.format(rf, ','.join([str(af) for af in afs])))
        print(record, end='')

if __name__ == "__main__":
    main()
