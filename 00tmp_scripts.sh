
###
 # @Author: Shuo Zhao && 18904530325@163.com
 # @Date: 2026-08-19 14:50:27
 # @LastEditors: Shuo Zhao && 18904530325@163.com
 # @LastEditTime: 2026-08-22 16:13:13
 # @FilePath: /Code_Notes/00tmp_scripts.sh
 # @Description: 
 # 
### 

#!/bin/bash
HOME="/public/home/zhaoshuo/work1"
GBZ="$HOME/01_work/02_graph/02_GiraffePrep/all_autoindex.giraffe.gbz"
DIST="$HOME/01_work/02_graph/02_GiraffePrep/all_autoindex.dist"
MIN="$HOME/01_work/02_graph/02_GiraffePrep/all_autoindex.shortread.withzip.min"
ZIP="$HOME/01_work/02_graph/02_GiraffePrep/all_autoindex.shortread.zipcodes"
FQ_DIR="$HOME/01_work/00_data/WGS_diploid" #re-sequencing
SAMPLE_LIST="$HOME/01_work/00_data/lq_diploid/list" #re-sequencing
OUT_DIR="$HOME/01_work/03_genotyping/01_giraffe"

for sample in $(cat $SAMPLE_LIST); do
    sed "s/50/64/g" work.sh | \
    sed "s/edta/giraffe_$sample/g" | \
    sed "s/%j/giraffe_$sample/g" > giraffe_$sample.sh

cat >> giraffe_$sample.sh << EOF
OUT="$OUT_DIR/${sample}"
mkdir -p \$OUT
F1=$FQ_DIR/${sample}/*_1.clean.fq.gz
F2=$FQ_DIR/${sample}/*_2.clean.fq.gz

vg giraffe -t 64 \
        -Z $GBZ \
        -d $DIST \
        -m $MIN \
        -z $ZIP \
        -f \$F1 \
        -f \$F2 \
        -o gam > \$OUT/${sample}_giraffe.gam
EOF
sbatch giraffe_$sample.sh
done

#!/bin/bash
HOME="/public/home/zhaoshuo/work1"
XG="$HOME/01_work/02_graph/02_GiraffePrep/all_autoindex.xg"
GBZ="$HOME/01_work/02_graph/02_GiraffePrep/all_autoindex.giraffe.gbz"
GAM_DIR="$HOME/01_work/03_genotyping/01_giraffe"
SAMPLE_LIST="$HOME/01_work/00_data/lq_diploid/list" #re-sequencing
OUT_DIR="$HOME/01_work/03_genotyping/02_pack"

for sample in $(cat $SAMPLE_LIST); do
    sed "s/50/32/g" work.sh | \
    sed "s/edta/pack_$sample/g" | \
    sed "s/%j/pack_$sample/g" > pack_$sample.sh

cat >> pack_$sample.sh << EOF
GAM="$GAM_DIR/${sample}/${sample}_giraffe.gam"
OUT="$OUT_DIR/${sample}"
mkdir -p \$OUT
vg pack -t 32 -Q 5 \
        -x $GBZ \
        -g \$GAM \
        -o "\$OUT/${sample}.pack"
EOF
sbatch pack_$sample.sh
done

#!/bin/bash
#SBATCH --partition=AMD_9A14
#SBATCH --cpus-per-task=128
#SBATCH --job-name=SVmerge
#SBATCH --output=%x.out
#SBATCH --error=%x.err

HOME="/public/home/zhaoshuo/work1"
REF="$HOME/data/reference/DM8.1_genome.ori.chr.fa"
SAMPLE_LIST="$HOME/01_work/00_data/lq_diploid/list" #re-sequencing
VCF_DIR="$HOME/01_work/03_genotyping/03_PPG1.2_call"

file_list="$VCF_DIR/vcf_list.txt"
> "$file_list"
while read -r sample; do
    [ -z "$sample" ] && continue
    echo "$VCF_DIR/${sample}/${sample}_q30.vcf" >> "$file_list"
done < "$SAMPLE_LIST"

SURVIVOR merge "$file_list" 1000 1 1 1 0 50 "$VCF_DIR/lq_PGG_1.2.vcf"
bgzip "$VCF_DIR/lq_PGG_1.2.vcf"
bcftools annotate -x INFO,^FORMAT/GT $VCF_DIR/lq_PGG_1.2.vcf.gz -Oz -o $VCF_DIR/lq_PGG_1.2.GT_only.vcf.gz
bcftools index -t "$VCF_DIR/lq_PGG_1.2.GT_only.vcf.gz"
bcftools sort "$VCF_DIR/lq_PGG_1.2.GT_only.vcf.gz" -Ov -o $VCF_DIR/lq_PGG_1.2_sorted.vcf
bgzip "$VCF_DIR/lq_PGG_1.2_sorted.vcf"
bcftools index -t "$VCF_DIR/lq_PGG_1.2_sorted.vcf.gz"

#!/bin/bash
#SBATCH --partition=AMD_9A14
#SBATCH --mem=64G
#SBATCH --cpus-per-task=64
#SBATCH --job-name=maf
#SBATCH --output=%x.out
#SBATCH --error=%x.err

HOME="/public/home/zhaoshuo/work1"
REF="$HOME/data/reference/DM8.1_genome.ori.chr.fa"
CHR_FILE="$HOME/01_work/00_data/01_116VCF/chrs.txt"
VCF="$HOME/01_work/03_genotyping/03_PPG1.2_call/lq_PGG_1.2_sorted.vcf.gz"
OUT_DIR="$HOME/01_work/03_genotyping/04_SV"

TMP_DIR="$HOME/01_work/03_genotyping/04_SV/tmp"
mkdir -p $TMP_DIR

export TMP_DIR VCF OUT_DIR
cat "$CHR_FILE" | xargs -P 64 -I {} bash -c '
    chr="{}"
    out_vcf="$TMP_DIR/${chr}.vcf.gz"

    bcftools view -r "$chr" -Oz -o "$out_vcf" "$VCF"
    bcftools index -t "$out_vcf"
'
cat "$CHR_FILE" | xargs -P 64 -I {} bash -c '
    chr="{}"
    out_vcf="$OUT_DIR/${chr}.maf005.MF05.vcf.gz"

    vcftools --gzvcf $TMP_DIR/${chr}.vcf.gz --maf 0.05 --max-missing 0.5 --recode --recode-INFO-all --stdout | bgzip > $out_vcf
    bcftools index -t "$out_vcf"
'

#!/bin/bash
#SBATCH --partition=AMD_9654
#SBATCH --cpus-per-task=64
#SBATCH --job-name=merge
#SBATCH --output=%x.out
#SBATCH --error=%x.err

HOME="/public/home/zhaoshuo/work1"
REF="$HOME/data/reference/DM8.1_genome.ori.chr.fa"
SAMPLE_LIST="$HOME/01_work/05_GWAS/phenotype_info/SGA_list" #re-sequencing
VCF_DIR="$HOME/01_work/03_genotyping/03_PPG1.2_call"
CHR_FILE="$HOME/01_work/00_data/01_116VCF/chrs.txt"
OUT_DIR="$HOME/01_work/05_GWAS/01_SV_GWAS/01_SGA/01_PPG1.2_VCF"

mkdir -p $OUT_DIR

file_list="$VCF_DIR/vcf_list.txt"
> "$file_list"
while read -r sample; do
    [ -z "$sample" ] && continue
    echo "$VCF_DIR/${sample}/${sample}_q30.vcf" >> "$file_list"
done < "$SAMPLE_LIST"

SURVIVOR merge "$file_list" 1000 1 1 1 0 50 "$OUT_DIR/SGA_PPG1.2_SV.vcf"
bgzip "$OUT_DIR/SGA_PPG1.2_SV.vcf"
bcftools annotate -x INFO,^FORMAT/GT $OUT_DIR/SGA_PPG1.2_SV.vcf.gz -Oz -o $OUT_DIR/SGA_PPG1.2_SV.GT_only.vcf.gz
bcftools sort "$OUT_DIR/SGA_PPG1.2_SV.GT_only.vcf.gz" -Oz -o $OUT_DIR/SGA_PPG1.2_SV_sorted.vcf.gz
bcftools index -t "$OUT_DIR/SGA_PPG1.2_SV_sorted.vcf.gz"
rm -rf $OUT_DIR/SGA_PPG1.2_SV.vcf.gz $OUT_DIR/SGA_PPG1.2_SV.GT_only.vcf.gz
vcftools --gzvcf $OUT_DIR/SGA_PPG1.2_SV_sorted.vcf.gz --maf 0.05 --recode --recode-INFO-all --stdout | bgzip > $OUT_DIR/SGA_PPG1.2_SV_maf005.vcf.gz
zcat $OUT_DIR/SGA_PPG1.2_SV_maf005.vcf.gz | awk '{if($0 ~ /^#/) print; else {gsub(/\.\/\./, "0/0"); print}}' | bgzip > $OUT_DIR/SGA_PPG1.2_SV_maf005_modified.vcf.gz
bcftools index -t $OUT_DIR/SGA_PPG1.2_SV_maf005_modified.vcf.gz

zcat SGA_PPG1.2_SV_sorted.vcf.gz | awk '{if($0 ~ /^#/) print; else {gsub(/\.\/\./, "0/0"); print}}' | bgzip > SGA_PPG1.2_SV_sorted_modified.vcf.gz

#!/bin/bash
#SBATCH --partition=AMD_9A14
#SBATCH --cpus-per-task=16
#SBATCH --job-name=plink
#SBATCH --output=%x.out
#SBATCH --error=%x.err

HOME="/public/home/zhaoshuo/work1"
CHR_FILE="$HOME/01_work/00_data/01_116VCF/chrs.txt"
VCF_DIR="$HOME/01_work/05_GWAS/01_SV_GWAS/01_SGA/01_PPG1.2_VCF"
PHENO_FILE="$HOME/01_work/05_GWAS/phenotype_info/SGA_pheno.txt"
OUT_DIR="$HOME/01_work/05_GWAS/01_SV_GWAS/01_SGA/02_PPG1.2_PLINK"

mkdir -p $OUT_DIR

# 转换
plink --vcf $VCF_DIR/SGA_PPG1.2_SV_sorted_modified_maf005.vcf.gz --make-bed --out $OUT_DIR/all --allow-extra-chr --set-missing-var-ids @:#
# LD 修剪
#plink --bfile $OUT_DIR/all --indep-pairwise 50kb 1 0.2 --out $OUT_DIR/all_pruned --allow-extra-chr
#plink --bfile $OUT_DIR/all --extract $OUT_DIR/all_pruned.prune.in --make-bed --out $OUT_DIR/all_pruned --allow-extra-chr
# 计算协变量
plink --bfile $OUT_DIR/all --pca 10 --allow-no-sex --out $OUT_DIR/pca_result
awk '{print 1, $3, $4, $5}' $OUT_DIR/pca_result.eigenvec > $OUT_DIR/gemma_cov.txt

vcftools --gzvcf GSA_SV_modified.vcf.gz --maf 0.05 --recode --recode-INFO-all --stdout | bgzip > GSA_SV_modified_maf005.vcf.gz

zcat GSA_SV_sorted.vcf.gz | awk '{if($0 ~ /^#/) print; else {gsub(/\.\/\./, "0/0"); print}}' | bgzip > GSA_SV_sorted_maf005_modified.vcf.gz