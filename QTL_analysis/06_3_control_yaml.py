'''
Author: Shuo Zhao && 18904530325@163.com
Date: 2026-09-21 14:27:21
LastEditors: Shuo Zhao && 18904530325@163.com
LastEditTime: 2026-09-21 16:57:54
FilePath: /02_work/04_QTL/06_3_control_yaml.py
Description: 

'''
#!/usr/bin/python
# Usage:
#    python 06_3_control_yaml.py <06_1-geno_file> <06_2-pheno_file> <06_1-pmap_file> <output_yaml> <cross_type>
#    - 06_1-geno_file: Path to the genotype file.
#    - 06_2-pheno_file: Path to the phenotype file.
#    - 06_1-pmap_file: Path to the physical map file.
#    - output_yaml: Path to the output YAML control file.
#    - cross_type: Type of cross (e.g., f2).
#
import os
import sys


geno_file = sys.argv[1]
pheno_file = sys.argv[2]
pmap_file = sys.argv[3]
output_yaml = sys.argv[4]
cross_type = sys.argv[5]

geno_file = os.path.basename(geno_file)
pheno_file = os.path.basename(pheno_file)
pmap_file = os.path.basename(pmap_file)

yaml_content = f"""# R/qtl2 control file

crosstype: {cross_type}

geno:
  - {geno_file}

pheno: {pheno_file}

pmap:
  - {pmap_file}

sep: ","
na.strings: ["NA", "-", ""]
"""

with open(output_yaml, "w", encoding="utf-8") as f:
  f.write(yaml_content)