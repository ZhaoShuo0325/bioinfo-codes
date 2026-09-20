#!/usr/bin/python
# Usage:
#    python 03_bin_geno.py <01-win_file> <02-snp_geno_file> <output_detail_file> <output_geno_file>
#    - bin_breaks_file: File containing chromosome bin start and end positions (col 1 and col 2).
#    - snp_geno_file: Single-sample chromosome SNP genotype file containing chr (col 1), pos (col 2), and genotype state (col 3).
#    - output_detail_file: Detailed output file containing bin name, coordinates, SNP counts for each class, and final bin genotype call.
#    - output_geno_file: Clean output file containing only the final bin genotype sequence for genetic mapping.

import sys

dict_breaks = {}
with open(sys.argv[1], 'rt') as fin1:
  for line in fin1:
    cols = line.split()
    if cols:
      c_name = cols[0]
      if c_name not in dict_breaks:
        dict_breaks[c_name] = {'start': [], 'end': []}
      dict_breaks[c_name]['start'].append(int(cols[1]))
      dict_breaks[c_name]['end'].append(int(cols[2]))

# dict_chrom_pos[chr_name][(chr_name, pos)] = geno
dict_chrom_pos = {}
with open(sys.argv[2], 'rt') as fin2:
  for line in fin2:
    cols = line.split()
    if cols:
      c_name, pos, geno = cols[0], int(cols[1]), cols[2]
      if c_name not in dict_chrom_pos:
        dict_chrom_pos[c_name] = {}
      dict_chrom_pos[c_name][(c_name, pos)] = geno

fout1_path = sys.argv[3]
fout2_path = sys.argv[4]

with open(fout1_path, 'w') as fout1, open(fout2_path, 'w') as fout2:
  for target_chr in dict_breaks:
    if target_chr not in dict_chrom_pos:
      continue

    list_breaks = dict_breaks[target_chr]['start']
    list_end = dict_breaks[target_chr]['end']
    dict_pos = dict_chrom_pos[target_chr]

    for k, (start, end) in enumerate(zip(list_breaks, list_end)):
      va, vb, vh, num = 0, 0, 0, 0

      for (c, p), geno in dict_pos.items():
        if start <= p < end:
          if geno == 'A157':
            va += 1
            num += 1
          elif geno == 'E454':
            vb += 1
            num += 1
          elif geno == 'heter':
            vh += 1
            num += 1

      bin_name = f'bin{k+1}'
      if num == 0:
        fout1.write(
            f'{target_chr}\t{bin_name}\t{start}\t{end}\t{va}\t{vb}\t{vh}\t-\n'
        )
        fout2.write(f'{target_chr}\t{bin_name}\t-\n')
      else:
        ratio = float(va - vb) / num
        if ratio >= 0.5:
          geno_call = 'A157'
        elif ratio <= -0.5:
          geno_call = 'E454'
        else:
          geno_call = 'heter'

        fout1.write(
            f'{target_chr}\t{bin_name}\t{start}\t{end}\t{va}\t{vb}\t{vh}\t{geno_call}\n'
        )
        fout2.write(
            f'{target_chr}\t{bin_name}\t{geno_call}\n'
        )