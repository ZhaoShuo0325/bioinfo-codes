#!/usr/bin/env python
'''
Author: Shuo Zhao && 18904530325@163.com
Date: 2026-10-15 08:41:32
LastEditors: Shuo Zhao && 18904530325@163.com
LastEditTime: 2026-10-19 09:39:40
FilePath: /Code_Notes/Genome_Assembly/split_fasta.py
Description: 

'''
import sys
import argparse

parser = argparse.ArgumentParser(description="Extract and split chromosomes into H1 and H2 from a FASTA file.")
parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
parser.add_argument("-o1", "--output-h1", required=True, help="Output FASTA file for H1 (odd chromosomes)")
parser.add_argument("-o2", "--output-h2", required=True, help="Output FASTA file for H2 (even chromosomes)")

args = parser.parse_args()

input_fasta = args.input
output_h1 = args.output_h1
output_h2 = args.output_h2

sequences = {}
current_header = None
current_lines = []

with open(input_fasta, "r") as f:
    for line in f:
        line_str = line.strip()
        if line_str.startswith(">"):
            if current_header:
                sequences[current_header] = current_lines
            current_header = line_str[1:]
            current_lines = []
        else:
            if current_header:
                current_lines.append(line_str)
if current_header:
    sequences[current_header] = current_lines

h1_file = open(output_h1, "w")
h2_file = open(output_h2, "w")

for i in range(1, 25):
    old_name = f"HiC_scaffold_{i}"
    if old_name in sequences:
        new_idx = (i + 1) // 2
        new_name = f"chr{new_idx:02d}"
        
        target_file = h1_file if (i % 2 != 0) else h2_file
        
        target_file.write(f">{new_name}\n")
        for seq_line in sequences[old_name]:
            target_file.write(seq_line + "\n")

h1_file.close()
h2_file.close()