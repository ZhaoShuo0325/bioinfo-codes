#!/bin/bash
#SBATCH --partition=AMD_9654
#SBATCH --cpus-per-task=8
#SBATCH --job-name=geno_exp_vis
#SBATCH --output=%x.out
#SBATCH --error=%x.err

HOME="/public/home/zhaoshuo/work1"
FILE_LIST="$HOME/01_work/06_candidate/02_zz/04_genotype_exp/file_list.txt"
Rscript="/public/home/zhaoshuo/miniconda3/envs/R_env/bin/Rscript"
OUT_DIR="$HOME/01_work/06_candidate/02_zz/04_genotype_exp"

cd $OUT_DIR

while read -r filepath; do
    [ -z "$filepath" ] && continue
    "$Rscript" /public/home/zhaoshuo/work1/01_work/06_candidate/mt1186/00_scripts/geno_exp.R "$filepath"
done < "$FILE_LIST"