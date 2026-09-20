#!/usr/bin/python
#Usage:
#    python 02_snp_geno.py <parents_vcf> <f2_snp_file> <output_geno_file>
#    - parents_vcf: Parents SNP VariantSet file containing genotype information for both parents (col 10 and 11).
#    - f2_snp_file: F2 single-sample SNP file containing chr (col 1), pos (col 2), and GT (col 10).
#    - output_geno_file: Output file containing mapped genotype classifications ('heter', 'A157', or 'E454') for each SNP position.

import sys
import gzip

dict_par = {}
with gzip.open(sys.argv[1], 'rt') as fin1:
    for line in fin1:
        if not line.startswith('#'):
            cols = line.split()
            if len(cols) >= 11:
                dict_par[(cols[0], cols[1])] = [cols[9].split(':')[0], cols[10].split(':')[0]]

with gzip.open(sys.argv[2], 'rt') as fin2, open(sys.argv[3], 'w') as fout:
    for line in fin2:
        if not line.startswith('#'):
            cols = line.split()
            if len(cols) >= 10:
                loc = (cols[0], cols[1])
                if loc in dict_par:
                    child_gt = cols[9].split(':')[0]
                    p1_gt, p2_gt = dict_par[loc]
                    if child_gt in {'0/1', '1/0', '0|1', '1|0'}:
                        fout.write(f"{cols[0]}\t{cols[1]}\theter\n")
                    elif child_gt == p1_gt:
                        fout.write(f"{cols[0]}\t{cols[1]}\tA157\n")
                    elif child_gt == p2_gt:
                        fout.write(f"{cols[0]}\t{cols[1]}\tE454\n")