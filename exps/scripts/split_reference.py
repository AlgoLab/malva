import os, sys

import gzip

from Bio import SeqIO

def main():
    fa_path = sys.argv[1]

    path = os.path.dirname(os.path.realpath(fa_path))

    print("Indexing input fa...")
    fa = SeqIO.index(fa_path, "fasta")

    chroms = [str(c) for c in list(range(1,23))] + ['X', 'Y']

    for c in chroms:
        print(c)
        out = "{}/chr{}.fa".format(path, c)
        ref = fa[c]
        SeqIO.write(ref, out, "fasta")

        #with open(out, 'rb') as f_in, gzip.open('{}.gz'.format(out), 'wb') as f_out:
        #    f_out.writelines(f_in)
        #os.remove(out)

if __name__ == '__main__':
    main()
