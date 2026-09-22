'''
Author: Shuo Zhao && 18904530325@163.com
Date: 2026-09-21 14:08:22
LastEditors: Shuo Zhao && 18904530325@163.com
LastEditTime: 2026-09-22 10:43:49
FilePath: /Code_Notes/QTL_analysis/06_2_prep_pheno.py
Description: 

'''
#!/usr/bin/python
# Usage:
#    python 06_2_prep_pheno.py <input_pheno_file> <sample_list> <output_pheno_csv>
#    - input_pheno_file: Path to the input phenotype file (Excel format).
#    - sample_list: Path to the file containing the list of sample IDs.
#    - output_pheno_csv: Path to the output phenotype CSV file.
#
import sys
import pandas as pd

input_file = sys.argv[1]
sample_list_path = sys.argv[2]
output_report = sys.argv[3]

df = pd.read_excel(input_file)
df[df.columns[0]] = df[df.columns[0]].astype(str).str.replace('×', 'x').str.strip()
id_col, valid = df.columns[0], [df.columns[0]]

for c in df.columns[1:]:
  s = pd.to_numeric(df[c], errors="coerce")
  if df[c].notna().sum() and df[c].notna().sum() == s.notna().sum():
    df[c], valid = s, valid + [c]
  else:
    print(f"exclude: {c}")

df[valid].set_index(id_col).reindex(
    [l.strip() for l in open(sample_list_path) if l.strip()]
).reset_index().to_csv(output_report, index=False, na_rep="NA")