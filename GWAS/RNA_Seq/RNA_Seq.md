<!--
 * @Author: Shuo Zhao && 18904530325@163.com
 * @Date: 2026-08-27 11:18:14
 * @LastEditors: Shuo Zhao && 18904530325@163.com
 * @LastEditTime: 2026-08-28 09:20:46
 * @FilePath: /Code_Notes/RNA_Seq.md
 * @Description: 
 * 
-->

# 转录组分析 (RNA-Seq)

## 原始数据质控  
**详见 WGS-analysis** https://github.com/ZhaoShuo0325/bioinfo-codes/blob/main/WGS_analysis.md  
   - 使用 `fastp` 软件去接头序列
```bash
fastp -i r1.fq.gz -I r2.fq.gz -o r1.clean.fq.gz -O r2.clean.fq.gz -w 24 --detect_adapter_for_pe
```

## 比对到参考基因组
**HISAT2** 是二代测序 Reads 比对工具，用于将 RNA-seq 或 DNA 测序数据比对到参考基因组上。
1. 建立参考基因组索引  
```bash
# 只需运行一次
hisat2-build -p 8 DM8.1_genome.ori.chr.fa DM8.1_genome.ori.chr
```
2. 比对 Reads  
```bash
# 运行 HISAT2 比对，并排序建立索引
hisat2 -p 8 -x DM8.1_genome.ori.chr \
       -1 r1.clean.fq.gz -2 r2.clean.fq.gz \
       --summary-file align_summary.txt | \
samtools view -bS - | \
samtools sort -o ${sample}_sorted.bam
samtools index sorted.bam
```

## 表达量定量分析
**featureCounts** 通过统计比对到参考基因组上的 Reads 数量，转换成一个基因表达量矩阵。
1. 准备基因注释 gft 文件  
```bash
# 将原有的 gff3 文件进行格式转换
gffread DM8.1_gene.gff3 -T -o DM8.1_gene.gtf
```
2. 运行 featureCounts  
```bash
# featureCounts 可以输入多个 bam 文件，生成一个多样本的矩阵
featureCounts -T 16 -p \
              --countReadPairs \
              -t exon \
              -g gene_id \
              -a DM8.1_gene.gtf \
              -o gene_count_matrix.txt \
              *_sorted.bam
```
featureCounts 的输出文件样本名可能是 bam 文件路径名，需要进行转换
```R
df <- read.table("gene_count_matrix.txt", header = TRUE, row.names = 1, sep = "\t", comment.char = "#", check.names = FALSE)
# 在 featureCounts 输出中：前5列是 Chr, Start, End, Strand, Length，从第6列开始是各个样本的 BAM 路径
sample_cols <- 6:ncol(df)
colnames(df)[sample_cols] <- gsub(".*/([^/]+)_sorted\\.bam", "\\1", colnames(df)[sample_cols])
write.table(df, file = "gene_clean_matrix.txt", sep = "\t", quote = FALSE, col.names = NA)
```
3. 基因表达量标准化  
   - **FPKM (Fragments Per Kilobase Million)** 每百万比对片段中，比对到某基因、每千碱基长度的片段数  
```R
df <- read.table("gene_clean_matrix.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
# 提取长度和 Count 矩阵 第六列后为样本
lengths <- df$Length
counts <- as.matrix(df[, 6:ncol(df)])
# 计算 FPKM
fpkm <- t(t(counts / lengths) / colSums(counts, na.rm = TRUE)) * 1e9
write.table(fpkm, file = "gene_FPKM_matrix.txt", sep = "\t", quote = FALSE, col.names = NA)
```
   - **TPM (Transcripts Per Million)** 每百万转录本中，该基因所占转录本数量  
```R
df <- read.table("gene_clean_matrix.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
lengths <- df$Length
counts <- as.matrix(df[, 6:ncol(df)])
# 计算 TPM
rpk <- counts / (lengths / 1000)
tpm <- t(t(rpk) / colSums(rpk, na.rm = TRUE)) * 1e6
write.table(tpm, file = "gene_TPM_matrix.txt", sep = "\t", quote = FALSE, col.names = NA)
```
**过滤**：过滤所有样本中表达量均小于1的基因  
```bash
awk -F'\t' '
    NR==1 {print; next} 
    {
        max = 0
        for (i=2; i<=NF; i++) {
            if ($i+0 > max) max = $i+0
        }
        if (max >= 1) print
    }
' gene_TPM_matrix.txt > gene_TPM_filtered.txt
```

## 提取目标基因
1. 利用基因组注释文件，提取目标范围内全部基因
```bash
# 提取10号染色体，53.3-54.3Mb 范围内的基因
awk -F'\t' '
    $3=="gene" && ($1=="chr10" || $1=="10") && $4 <= 54300000 && $5 >= 53300000 {
        split($9, arr, ";");
        for (i in arr) {
            if (arr[i] ~ /^ID=/) {
                sub(/^ID=/, "", arr[i]);
                print arr[i]
            }
        }
    }
' DM8.1_gene.gff3 > C10_53.3_54.3_gene.txt
```
2. 利用过滤后表达量矩阵矩阵，提取有表达的候选基因
```bash
# 提取10号染色体，53.3-54.3Mb 范围内表达的基因
awk '
    NR==FNR {
        if (NR > 1) expr[$1] = 1
        next
    }
    $1 in expr {
        print $1
    }
' gene_TPM_filtered.txt C10_53.3_54.3_gene.txt > C10_53.3_54.3_candidate_gene.txt
# 输出这些基因的位置 bed 文件
awk -F'\t' '
    NR==FNR {
        genes[$1] = 1
        next
    }
    $3=="gene" {
        split($9, arr, ";")
        id = ""
        for (i in arr) {
            if (arr[i] ~ /^ID=/) {
                id = arr[i]
                sub(/^ID=/, "", id)
            }
        }
        if (id in genes) {
            chrom = $1
            start = $4 - 1
            end = $5
            strand = $7
            print chrom "\t" start "\t" end "\t" id "\t.\t" strand
        }
    }
' C10_53.3_54.3_candidate_gene.txt DM8.1_gene.gff3 > C10_53.3_54.3_candidate_gene.bed
```
3. 提取候选基因的表达量
```bash
awk '
    NR==FNR {
        genes[$1] = 1
        next
    }
    NR==1 {
        print
        next
    }
    $1 in genes {
        print
    }
' C10_53.3_54.3_candidate_gene.txt gene_TPM_filtered.txt > C10_53.3_54.3_candidate_expression.txt
```
4. 提取候选 SV
```bash
# 候选 SV：位于候选基因上下游 2000bp 内的 SV
# 生成新的 bed 文件
awk -F'\t' 'BEGIN {OFS="\t"} {
    start = $2 - 2000; 
    if (start < 0) start = 0; 
    end = $3 + 2000; 
    print $1, start, end
}' /public/home/zhaoshuo/work1/01_work/06_candidate/mt1186/03_candidate/C10_53.3_54.3_candidate_gene.bed > /public/home/zhaoshuo/work1/01_work/06_candidate/mt1186/03_candidate/C10_53.3_54.3_candidate_genes_slop2k.bed
# 提取 SV
bcftools view -R /public/home/zhaoshuo/work1/01_work/06_candidate/mt1186/03_candidate/C10_53.3_54.3_candidate_genes_slop2k.bed \
              -Oz -o /public/home/zhaoshuo/work1/01_work/06_candidate/mt1186/03_candidate/C10_53.3_54.3_candidate_SV.vcf.gz \
              /public/home/zhaoshuo/work1/01_work/05_GWAS/01_SV_GWAS/02_zz/01_PPG1.4_VCF/zz_PPG1.4_SV_sorted_modified.vcf.gz
bcftools index -t /public/home/zhaoshuo/work1/01_work/06_candidate/mt1186/03_candidate/C10_53.3_54.3_candidate_SV.vcf.gz
# 可以只保留 maf > 0.05 的 SV
vcftools --gzvcf C10_53.3_54.3_candidate_SV.vcf.gz --maf 0.05 --recode --recode-INFO-all --stdout | bgzip -c > C10_53.3_54.3_candidate_SV_maf0.05.vcf.gz
bcftools index -t C10_53.3_54.3_candidate_SV_maf0.05.vcf.gz
# SV 重命名
bcftools view C10_53.3_54.3_candidate_SV_maf0.05.vcf.gz | awk 'BEGIN{OFS="\t"} /^#/ {print; next} {$3="sv" ++i; print}' | bgzip -c > C10_53.3_54.3_candidate_SV_maf0.05.vcf.gz.tmp && mv C10_53.3_54.3_candidate_SV_maf0.05.vcf.gz.tmp C10_53.3_54.3_candidate_SV_maf0.05.vcf.gz && bcftools index -t C10_53.3_54.3_candidate_SV_maf0.05.vcf.gz
```
5. 候选基因与候选 SV 的对应关系
```bash
# 1. 将 SV VCF 转换为临时 BED 文件（处理 0-based 坐标）
bcftools query -f '%CHROM\t%POS\t%POS\t%ID\n' C10_53.3_54.3_candidate_SV_maf0.05.vcf.gz | awk 'BEGIN {OFS="\t"} {print $1, $2-1, $2, $4}' > /tmp/sv_temp.bed

# 2. 使用 bedtools window 寻找上下游 2000bp 内的对应关系，并输出 Gene 与 SV ID
bedtools window -a C10_53.3_54.3_candidate_gene.bed -b /tmp/sv_temp.bed -w 2000 | awk 'BEGIN {OFS="\t"; print "Gene_ID", "SV_ID"} {print $4, $10}' > C10_53.3_54.3_gene_sv_mapping.txt

# 3. 清理临时文件
rm /tmp/sv_temp.bed
```

## 目标基因表达与代谢物含量关联分析  

1. 生成目标基因表达列表
```bash
# 每个目标基因输出一个文件，第一列样本名，第二列 TPM
SAMPLE_LIST="zz88_list" # 样本文件
GENE_LIST="../01_corralation/GAME_gene.txt" # 目标基因文件
TPM_FILE="zz88_TPM_filtered.txt" # 表达水平文件
OUT_DIR="../01_corralation/01_GAME_mt" #输出目录

# 创建输出目录
mkdir -p "$OUT_DIR"

# 清理样本列表和基因列表的 \r，并去除基因名的版本号（如 .1）
tr -d '\r' < "$SAMPLE_LIST" > .clean_samples.txt
tr -d '\r' < "$GENE_LIST" | awk -F'.' '{print $1}' > .clean_genes.txt

while IFS= read -r gene; do
    [[ -z "$gene" ]] && continue
    
    awk -v target_gene="$gene" -v sample_file=".clean_samples.txt" '
    BEGIN {
        while ((getline s < sample_file) > 0) {
            if (s != "") samples[s] = 1
        }
        close(sample_file)
    }
    NR == 1 {
        for (i=1; i<=NF; i++) {
            header[i] = $i
            sub(/\r$/, "", header[i])
            if (header[i] in samples) {
                valid_cols[i] = 1
            }
        }
    }
    {
        g = $1
        sub(/\r$/, "", g)
        sub(/\..*$/, "", g)
    }
    g == target_gene {
        for (i=2; i<=NF; i++) {
            if (i in valid_cols) {
                print header[i] "\t" $i
            }
        }
    }
    ' "$TPM_FILE" > "$OUT_DIR/${gene}.txt"
    
done < .clean_genes.txt

# 清理临时文件
rm -f .clean_samples.txt .clean_genes.txt

# 一些被过滤掉的基因在TPM表中找不到，因此删除 0 字符文件
#find . -name "*.txt" -size 0 -delete
```

2. 生成目标代谢物含量总表
```bash
# 准备目标代谢物表型文件
while IFS= read -r mt; do
    mt=$(echo "$mt" | tr -d '\r')
    [[ -z "$mt" ]] && continue
    src="/public/home/zhaoshuo/work1/01_work/05_GWAS/03_mt_SV/01_all_pheno/${mt}_pheno.txt"
    if [ -f "$src" ]; then
        cp "$src" ./
    else
        echo "Warning: $src not found"
    fi
done < ../mt_SGA.txt
```
```python
import os
import pandas as pd

# 读取目标样本列表
with open("zz88_list") as f:
    target_samples = [line.strip() for line in f if line.strip()]

# 读取目标代谢物列表
with open("../mt_SGA.txt") as f:
    target_mts = [line.strip() for line in f if line.strip()]

# 读取原始Excel表型数据
df = pd.read_excel("../01_all_pheno/ZZ_metabolome_processed.xlsx")
start_col_index = 11  # 起始样本列索引

# 提取Excel中的所有样本列名（假设第12列开始是样本数据，表头或行名需对应；这里根据参考脚本处理转置后的列）
# 先将Excel转为：行为样本，列为代谢物的矩阵
df_subset = df.iloc[:, [0] + list(range(start_col_index, df.shape[1]))]
df_subset = df_subset.dropna(subset=[df_subset.columns[0]])
df_subset.iloc[:, 0] = df_subset.iloc[:, 0].astype(str).str.strip()

# 转置：行变成样本，列变成代谢物
transposed_df = df_subset.set_index(df_subset.columns[0]).T
transposed_df.index = transposed_df.index.astype(str).str.strip()

# 筛选目标代谢物（只保留在Excel中存在的代谢物）
valid_mts = [mt for mt in target_mts if mt in transposed_df.columns]

# 筛选目标样本（只保留在Excel中存在的样本）
valid_samples = [s for s in target_samples if s in transposed_df.index]

# 提取子矩阵并按目标样本列表的顺序排列
final_df = transposed_df.reindex(index=valid_samples, columns=valid_mts)

# 填充缺失值为 NA
final_df = final_df.fillna("NA")

# 保存结果
output_path = "combined_mt_pheno_values.txt"
final_df.to_csv(output_path, sep="\t", header=True, index=True, index_label="Sample")
print(f"提取完成，已保存至 {output_path}")
```