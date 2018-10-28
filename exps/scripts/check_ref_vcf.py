import sys

from Bio import SeqIO
from pysam import VariantFile

def main():
    fa_path = sys.argv[1]
    vcf_path = sys.argv[2]

    ref = SeqIO.read(fa_path, "fasta")
    
    vcf = VariantFile(vcf_path)

    total = 0
    wrong = 0
    for var in vcf.fetch():
        total+=1
        ref_all = var.ref
        pos = var.pos - 1 #VCF is 1-based

        ref_ss = ref[pos:pos+len(ref_all)].seq
        if ref_all != ref_ss:
            wrong+=1
    print("{}/{} are wrong".format(wrong, total))

if __name__ == '__main__':
    main()
