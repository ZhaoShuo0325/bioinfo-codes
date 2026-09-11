#!/bin/bash
#SBATCH --partition=AMD_9654
#SBATCH --cpus-per-task=1
#SBATCH --job-name=genotype_exp
#SBATCH --output=%x.out
#SBATCH --error=%x.err

HOME="/public/home/zhaoshuo/work1"
VCF="$HOME/01_work/06_candidate/02_zz/03_candidate/C07_44.0_45.0_candidate_SV_maf0.05.vcf.gz"
SAMPLE_LIST="$HOME/01_work/06_candidate/mt1186/02_featureCounts/zz88_list"
MAPPING_FILE="$HOME/01_work/06_candidate/02_zz/03_candidate/C07_44.0_45.0_gene_sv_mapping.txt"
TPM="$HOME/01_work/06_candidate/02_zz/02_featureCounts/zz88_TPM_filtered.txt"
OUT_DIR="$HOME/01_work/06_candidate/02_zz/04_genotype_exp"

tail -n +2 "$MAPPING_FILE" | while read -r gene sv; do
    # 1. 使用 -S 参数仅提取 SAMPLE_LIST 中指定样本的该 SV 基因型
    bcftools query -S "$SAMPLE_LIST" -i "ID=='$sv'" -f '[%SAMPLE\t%GT\n]' "$VCF" > /tmp/sv_geno_$$.txt

    # 2. 从 TPM 矩阵中提取该基因在各个样本中的表达量 (格式: 样本名 \t 表达量)
    awk -v target="$gene" '
        NR==1 {
            for (i=2; i<=NF; i++) samples[i] = $i
        }
        $1 == target {
            for (i=2; i<=NF; i++) print samples[i], $i
        }
    ' "$TPM" > /tmp/gene_exp_$$.txt

    # 3. 按样本名进行合并，输出：样本名、基因型、表达量
    join -1 1 -2 1 <(sort /tmp/sv_geno_$$.txt) <(sort /tmp/gene_exp_$$.txt) \
        | awk '{print $1, $2, $3}' \
        > "$OUT_DIR/${gene}_${sv}.txt"

    # 清理临时文件
    rm -f /tmp/sv_geno_$$.txt /tmp/gene_exp_$$.txt
done

#find "$PWD" -maxdepth 1 -name "*.txt" ! -name "file_list.txt" > file_list.txt