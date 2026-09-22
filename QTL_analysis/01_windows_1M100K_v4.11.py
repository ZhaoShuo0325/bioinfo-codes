'''
Author: Shuo Zhao && 18904530325@163.com
Date: 2026-09-19 16:38:41
LastEditors: Shuo Zhao && 18904530325@163.com
LastEditTime: 2026-09-22 10:42:48
FilePath: /Code_Notes/QTL_analysis/01_windows_1M100K_v4.11.py
Description: 

'''
#!/usr/bin/python
#Usage:
#    python 01_windows_1M100K_v4.11.py <input_parents_file> <input_snp_file> <output_win_file>
#    - input_parents_file: Parents SNP VariantSet file
#    - input_snp_file: Preprocessed table file containing chr (col 1), pos (col 2), and GT (col 10).
#    - output_win_file: Output file for sliding window heterozygosity (1Mb window, 100kb step).
import sys
import gzip

valid_set = set()
opener = gzip.open if sys.argv[1].endswith('.gz') else open
with opener(sys.argv[1], 'rt') as f:
    for line in f:
        if not line.startswith('#'):
            cols = line.split()
            if len(cols) >= 2:
                valid_set.add((cols[0], cols[1]))

dict_all = {}
with gzip.open(sys.argv[2], 'rt') as f:
    for line in f:
        if not line.startswith('#'):
            cols = line.split()
            if len(cols) >= 10 and (cols[0], cols[1]) in valid_set:
                dict_all.setdefault(cols[0], {})[cols[1]] = cols[9].split(':')[0]

chr_lens = {
    'chr01': 88922896, 'chr02': 46362011, 'chr03': 61274481,
    'chr04': 69827751, 'chr05': 55371468, 'chr06': 59724636,
    'chr07': 58503165, 'chr08': 60055318, 'chr09': 68052101,
    'chr10': 61888340, 'chr11': 47470937, 'chr12': 60415087
}
het_types = {'0/1', '1/0', '0|1', '1|0'}
hom_types = {'0/0', '1/1', '0|0', '1|1'}

with open(sys.argv[3], 'w') as fout:
    for h, length in chr_lens.items():
        if h not in dict_all: continue
        dict_ch = dict_all[h]
        
        for start in range(0, length, 100000):
            end = min(start + 1000000, length)
            n = heter = 0
            
            for k_pos, gt in dict_ch.items():
                if start <= int(k_pos) <= end:
                    if gt in het_types:
                        n += 1
                        heter += 1
                    elif gt in hom_types:
                        n += 1
                        
            if n > 0:
                fout.write(f"{h}\t{start}\t{end}\t{float(heter)/n:.2f}\n")