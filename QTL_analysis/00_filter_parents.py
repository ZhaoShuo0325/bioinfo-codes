#!/usr/bin/env python3
import gzip
import sys

vcf_in = sys.argv[1]
vcf_out = sys.argv[2]
n_tot = n_keep = 0

with gzip.open(vcf_in, "rt") as fin, open(vcf_out, "w") as fout:
    for line in fin:
        if line.startswith("#"):
            if line.startswith("#CHROM"):
                h = line.strip().split("\t")
                iA, iE, iY = h.index("A157"), h.index("E454"), h.index("YS1")
            fout.write(line)
            continue

        n_tot += 1
        row = line.strip().split("\t")
        gA, gE, gY = (
            row[iA][:3].replace("|", "/"),
            row[iE][:3].replace("|", "/"),
            row[iY][:3].replace("|", "/"),
        )

        if gY in ("0/1", "1/0") and (
            (gA, gE) == ("0/0", "1/1") or (gA, gE) == ("1/1", "0/0")
        ):
            fout.write(line)
            n_keep += 1

print(f"Total: {n_tot}, Kept: {n_keep}")