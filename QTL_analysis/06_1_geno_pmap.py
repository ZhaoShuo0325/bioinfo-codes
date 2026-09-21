'''
Author: Shuo Zhao && 18904530325@163.com
Date: 2026-09-21 12:38:19
LastEditors: Shuo Zhao && 18904530325@163.com
LastEditTime: 2026-09-21 16:58:00
FilePath: /02_work/04_QTL/06_1_geno_pmap.py
Description: 

'''
#!/usr/bin/python
# Usage:
#    python 06_1_geno_pmap.py <input_F2_matrix> <sample_list> <output_geno> <output_pmap>
#    - input_matrix: Unified bin genotype matrix file (from previous step).
#    - sample_list: Path to the file containing the list of sample IDs.
#    - output_geno: Output genotype CSV file formatted for R/qtl2 (geno.csv).
#    - output_pmap: Output physical map CSV file for R/qtl2 (pmap.csv).
#
import sys
import pandas as pd

input_matrix = sys.argv[1]
sample_list_path = sys.argv[2]
output_geno = sys.argv[3]
output_pmap = sys.argv[4]

with open(sample_list_path, "r") as f:
  samples = [line.strip() for line in f if line.strip()]

df = pd.read_csv(
    input_matrix,
    sep=r"\s+",
    header=None,
    names=["chr", "marker", "start", "end"] + samples,
)

df["marker"] = df["chr"].astype(str) + "_" + df["marker"].astype(str)
df["pos"] = (df["start"] + df["end"]) / 2.0 / 1e6
df[["marker", "chr", "pos"]].to_csv(output_pmap, index=False)

geno_transposed = (
    df[samples]
    .replace({"A157": "A", "E454": "B", "heter": "H", "-": pd.NA})
    .T
)
geno_transposed.columns = df["marker"].values
geno_transposed.index.name = "id"
geno_transposed.to_csv(output_geno)