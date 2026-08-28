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

4. 样本-基因型-表达量 关联分析  
```bash
# 提取样本名
awk -F'\t' 'NR==1 {for(i=2;i<=NF;i++) print $i}' gene_TPM_matrix.txt > sample_list
```
```bash
bcftools query -r chr10:53842331 -S gene_list -f '[%SAMPLE\t%GT\n]' candidate_SV.vcf.gz > DM8C10G21370_expression.txt
```
```python
python3 -c '
import pandas as pd
file_path = "MYB210_TPM.txt"
df = pd.read_csv(file_path, sep="\t", header=None)
tpm = pd.read_csv("zz88_FPKM_matrix.txt", sep="\t", index_col=0)
df[2] = df[0].map(tpm.loc["DM8C10G21210"])
df.to_csv(file_path, sep="\t", header=False, index=False)
'
```