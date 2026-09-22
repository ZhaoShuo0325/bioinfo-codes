'''
Author: Shuo Zhao && 18904530325@163.com
Date: 2026-09-21 12:36:53
LastEditors: Shuo Zhao && 18904530325@163.com
LastEditTime: 2026-09-22 10:43:37
FilePath: /Code_Notes/QTL_analysis/04_heter_geno.py
Description: 

'''
#!/usr/bin/python
# Usage:
#    python 04_heter_geno.py <03-input_detail_file> <output_report_file> <output_geno_file>
#    - 03-input_detail_file: Detailed bin genotype file
#    - output_report_file: Output summary statistics report for this sample.
#    - output_geno_file: Output file containing only the final bin genotype.

import sys

input_file = sys.argv[1]
output_report = sys.argv[2]
output_geno = sys.argv[3]

total_bins = 0
a157_bins = 0
e454_bins = 0
heter_bins = 0

with open(input_file, 'rt') as fin, open(output_geno, 'wt') as fout_geno:
  for line in fin:
    cols = line.split()
    if cols:
      total_bins += 1
      geno_call = cols[7]

      fout_geno.write(f'{geno_call}\n')

      if geno_call == 'A157':
        a157_bins += 1
      elif geno_call == 'E454':
        e454_bins += 1
      elif geno_call == 'heter':
        heter_bins += 1

a157_ratio = a157_bins / total_bins
e454_ratio = e454_bins / total_bins
heter_ratio = heter_bins / total_bins

with open(output_report, 'wt') as fout_rep:
  fout_rep.write(f'Total_Bins\t{total_bins}\n')
  fout_rep.write(f'A157_Bins\t{a157_bins}\t{a157_ratio:.4f}\n')
  fout_rep.write(f'E454_Bins\t{e454_bins}\t{e454_ratio:.4f}\n')
  fout_rep.write(f'Heter_Bins\t{heter_bins}\t{heter_ratio:.4f}\n')