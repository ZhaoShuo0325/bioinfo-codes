<!--
 * @Author: Shuo Zhao && 18904530325@163.com
 * @Date: 2026-08-12 10:37:32
 * @LastEditors: Shuo Zhao && 18904530325@163.com
 * @LastEditTime: 2026-08-28 09:20:53
 * @FilePath: /Code_Notes/GWAS/GWAS.md
 * @Description: 
 * 
-->

# GWAS(Genome-Wide Association Study)
**GWAS（Genome-Wide Association Study，全基因组关联分析）** 是在全基因组范围内，寻找与特定表型显著关联的遗传变异位点，如SNP、InDel、SV等。  
**GWAS的理论基础——连锁不平衡（Linkage Disequilibrium, LD）** 如果某个变异位点紧邻或位于控制某个性状的基因上，那么带有这个特定等位基因的个体，其表型往往会表现出某种规律性的偏好，通过对群体中大量个体的全基因组变异进行扫描与分析，可以计算出每个位点与表型之间的统计学关联。  
**GWAS分析的核心流程** 包括：
   - **变异检测（Variant Calling）**：利用重测序数据获得群体样本的变异数据VCF文件。
   - **数据质控与过滤（Quality Control）**：剔除低质量位点、高缺失率样本和缺失率过低的位点等。
   - **基因型填充（Imputation）**：利用Beagle等软件对缺失的基因型进行填充。
   - **群体结构与亲缘关系分析（Population Structure & Kinship）**：计算主成分（PCA）和亲缘关系矩阵（Kinship），确保群体分类正确.
   - **关联分析（Association Mapping）**：利用GEMMA、EMMAX等软件，采用全基因组混合线性模型等进行全基因组关联分析。
   - **可视化与候选基因挖掘（Visualizaton & Candidate Genes）**：绘制曼哈顿图、QQ图、LD热图等寻找显著位点，并注释候选基因。

## 变异检测（Variant Calling）
### SNP caller
#### GATK (Genome Analysis Toolkit)
1. 原始数据质控与预处理  
**详见 WGS-analysis** https://github.com/ZhaoShuo0325/bioinfo-codes/blob/main/WGS_analysis.md  
   - 使用 `fastp` 软件去接头序列
   - 使用 `bwa` 软件进行序列比对，获得 bam 文件

2. 运行 GATK
**gatk CreateSequenceDictionary** 生成参考基因组 `.dict` 文件，必须且只需一次
```bash
gatk CreateSequenceDictionary -R ref.fa -O ref.dict
```

**gatk MarkDuplicates** 去除 PCR 重复序列（与 Picard 作用等价）
```bash
gatk MarkDuplicates -I ${sample}.bam -O ${sample}.markdup.bam -M ${sample}.markdup_metrics.txt
```

**gatk HaplotypeCaller** 单样本变异检测，建议分染色体进行
```bash
# --emit-ref-confidence GVCF 参数：记录全基因组所有位点的置信度
gatk HaplotypeCaller \
    -R $REF \
    -I ${sample}.markdup.bam \
    -O ${sample}_${chr}.g.vcf \
    -L $chr \
    --emit-ref-confidence GVCF \
    --native-pair-hmm-threads 32
```
可以将多条染色体合并  
```bash
# 获取所有染色体 gvcf 路径文件
list_file="${sample}_files.list"
for gvcf in $(ls "${sample}_chr"*.g.vcf.gz | sort -V); do
    echo "$gvcf" >> "$list_file"
done
# 合并 建立索引
bcftools concat -f "$list_file" --threads 32 -Oz -o "${sample}.g.vcf.gz"
bcftools index -t --threads 32 ${sample}.g.vcf.gz
```

**gatk CombineGVCFs** 多样本合并 GVCF 文件，推荐分染色体运行
```bash
# 生成传入参数文件 -V /path/to/gvcf
> $GVCF_LIST
for sample in $(cat $SAMPLE_LIST); do
    gvcf_path="$VCF_DIR/${sample}/${sample}_${chr}.g.vcf.gz"
        echo "-V $gvcf_path" >> "$GVCF_LIST"
done
# gatk CombineGVCFs 合并
gatk CombineGVCFs \
    -R $REF \
    -O combined_${chr}.g.vcf.gz \
    $(cat $GVCF_LIST)
```

**gatk GenotypeGVCFs** 联合基因分型(Joint Genotyping)  
将群体中所有样本的数据综合评估，转换成标准的 VCF 变异基因型矩阵文件  
```bash
gatk GenotypeGVCFs \
    -R $REF \
    -V combined_${chr}.g.vcf.gz \
    -O combined_${chr}_raw.vcf.gz
```

**gatk VariantFiltration** **VCFTools** 变异过滤工具  
过滤分为 `硬过滤` 和 `群体水平过滤`
```bash
# 硬过滤（Variant Hard Filtering）
# 过滤条件：SNP："QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" were removed
bcftools view -v snps -m 2 -M 2 combined_${chr}_raw.vcf.gz -Oz -o combined_${chr}_snp.vcf.gz
bcftools index -t combined_${chr}_snp.vcf.gz
gatk VariantFiltration \
    -R $REF \
    -V combined_${chr}_snp.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "SNP_HardFilter" \
    -O ${chr}_HardFilter.vcf.gz 
```
```bash
# 群体水平过滤
# 过滤条件：missing rate > 40% && MAF < 0.01 were removed
vcftools --gzvcf ${chr}_HardFilter.vcf.gz --maf 0.01 --max-missing 0.6 --recode --recode-INFO-all --stdout | bgzip > ${chr}_maf001_ms06.vcf.gz
```
最后将所有分染色体的文件 concat 在一起，并建立索引
```bash
bcftools concat *_maf001_ms06.vcf.gz -Oz -o all_maf001_ms06.vcf.gz
bcftools index -t all_maf001_ms06.vcf.gz
# 可选择只保留 GT 信息
bcftools annotate -x FORMAT all_maf001_ms06.vcf.gz -Oz -o all_maf001_ms06_GT_only.vcf.gz
bcftools index SGA_maf001_ms06_GT_only.vcf.gz
```


#### DeepVariant
1. 原始数据质控与预处理  
**详见 WGS-analysis** https://github.com/ZhaoShuo0325/bioinfo-codes/blob/main/WGS_analysis.md  
   - 使用 `fastp` 软件去接头序列
   - 使用 `bwa` 软件进行序列比对，获得 bam 文件
   - 使用 `Picard` 软件去除 PCR 重复序列
   
2. 运行 deepvariant  
**建议关闭小模型 `small_model_call_multiallelics=false`**
```bash
singularity exec deepvariant-1.9.0.sif \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref=$REF \
  --reads=$IN_BAM \
  --output_gvcf=$OUT/${sample}.g.vcf.gz \
  --output_vcf=$OUT/${sample}.vcf.gz \
  --num_shards=32 \
  --intermediate_results_dir="$TMP" \
  --make_examples_extra_args="small_model_call_multiallelics=false" \
  --vcf_stats_report=true \
  --logging_dir=$OUT/deepvariant_logs
```

3. 使用 GLnexus 合并多样本 gVCF 文件  
**DeepVariant 推荐使用 GLnexus 进行多样本合并**
```bash
# 生成 vcf 路径文件
> "$SAMPLE_FILE"
while read -r sample; do
    [ -z "$sample" ] && continue
    vcf_path="${sample}.g.vcf.gz"
    if [ -f "$vcf_path" ]; then
        echo "$vcf_path" >> "$SAMPLE_FILE"
    fi
done < "$SAMPLE_LIST"
```
```bash
#运行 GLnexus 建议分染色体合并
for i in {1..12}; do
    chr="chr$(printf "%02d" $i)"
    grep -w "$chr" ref.fa.fai | awk '{print $1 "\t0\t" $2}' > "${chr}.bed"
done
singularity exec glnexus-1.4.1.sif glnexus_cli \
    --config DeepVariant \
    --dir $DB_DIR \
    --threads 32 \
    --mem-gbytes 200 \
    --bed "${chr}.bed" \
    --list $SAMPLE_FILE > "${chr}_merged.bcf"
```

#### sentieon
**sentieon** 是高效率 GATK 替代软件，计算结果与 GATK 完全一致，但运行速度更快  
1. 重校正
```bash
# 对应 GATK BaseRecalibrator
# 分析原始 BAM 文件中测序 Reads 的碱基质量得分，生成质量统计表（重校正表）
sentieon driver -r $REF -t 16 -i ${sample}.bam --algo QualCal tmp/${sample}_recal_data.table
# 模拟重校正后的质量分布表
sentieon driver -r $REF -t 16 -i ${sample}.bam -q tmp/${sample}_recal_data.table --algo QualCal tmp/${sample}_recal_data.table.post
# 将校正前后的两份表进行对比，生成包含对比结果的 recal.csv 文件
sentieon driver -t 16 --algo QualCal --plot --before tmp/${sample}_recal_data.table --after tmp/${sample}_recal_data.table.post tmp/${sample}_recal.csv
# 绘制 BQSR 评估图表
sentieon plot bqsr -o tmp/${sample}_recal_plots.pdf tmp/${sample}_recal.csv
```

2. 应用BQSR 校正后 Call SNP  
```bash
# 通过 -q 动态引入第一步计算出的校正表，避免产生中间 bam 文件
sentieon driver -r $fasta -t $nt -i ${sample}.bam -q tmp/${sample}_recal_data.table --algo Haplotyper --emit_mode gvcf ${sample}_gvcf.gz
```

3. 群体 Joint Calling
```bash
> $GVCF_LIST
for sample in $(cat $SAMPLE_LIST); do
    gvcf_path="${sample}_${chr}.g.vcf.gz"
    if [ -f "$gvcf_path" ]; then
        echo "-v $gvcf_path" >> "$GVCF_LIST"
    else
        echo "Warning: $gvcf_path not found!"
    fi
done
sentieon driver \
    -t 32 \
    -r $REF \
    --algo GVCFtyper \
    $(cat $GVCF_LIST) \
    all_samples.vcf.gz
```

4. 过滤同 GATK

### SV caller
#### vg toolkit
1. 图泛基因组构建与 genotyping  
**详见vg pipeline**

2. SV 群体 VCF 文件合并  
**SURVIVOR**  
```bash
# 生成 VCF 文件路径文件
# 注意：必须是未压缩的 VCF 文件
> "$file_list"
while read -r sample; do
    [ -z "$sample" ] && continue
    echo "${sample}_q30.vcf" >> "$file_list"
done < "$SAMPLE_LIST"
```
```bash
# SURVIVOR 合并
SURVIVOR merge "$file_list" 1000 1 1 1 0 50 "survivor_merged.vcf"
bgzip "survivor_merged.vcf"
# 去除非 GT 字段 排序 填补0/0 过滤MAF 建立索引
bcftools annotate -x INFO,^FORMAT/GT survivor_merged.vcf.gz -Oz -o survivor_merged.GT_only.vcf.gz
bcftools sort survivor_merged.GT_only.vcf.gz -Oz -o survivor_merged_sorted.vcf.gz
# 将合并后缺失的 GT 添为 0/0
zcat survivor_merged_sorted.vcf.gz | awk '{if($0 ~ /^#/) print; else {gsub(/\.\/\./, "0/0"); print}}' | bgzip > modified_SV.vcf.gz
# 过滤 MAF < 0.05
vcftools --gzvcf modified_SV.vcf.gz --maf 0.05 --recode --recode-INFO-all --stdout | bgzip > modified_SV_maf005.vcf.gz
bcftools index -t modified_SV_maf005.vcf.gz
```

## 关联分析
### PCA 分析
**PLINK** 可以将 VCF 文件转换为 BED BIM FAM 等格式，用于处理大型数据集，进行 PCA 分析等
```bash
# 将 VCF 文件转换为 BED BIM FAM 格式，指定输出文件前缀
# --set-missing-var-ids @:# 给没有ID的位点自动命名
plink --vcf modified_SV.vcf.gz --make-bed --out modified_SV --set-missing-var-ids @:#
```

```bash
# 进行 PCA 分析，指定输出文件前缀和主成分数量
plink --bfile modified_SV --pca 10 --out modified_SV_pca
```
```bash
# 生成协变量文件
awk '{print 1, $3, $4, $5}' modified_SV_pca.eigenvec > gemma_cov.txt
```

### GEMMA
**GEMMA** 是全基因组混合线性模型软件，可以进行全基因组关联分析，包括单个位点关联分析和群体水平关联分析  
```bash
# 生成亲缘关系矩阵，校正群体结构
gemma -bfile $PLINK_DIR/modified_SV \
      -p "$PHENO_FILE" \
      -gk 1 \
      -outdir ./output \
      -o gwas_kinship
```
```bash
# 单个位点关联分析
gemma -bfile $PLINK_DIR/modified_SV \
      -k ./output/gwas_kinship.cXX.txt \
      -p $PHENO_FILE \
      -c $COV_FILE \
      -lmm 1 \
      -outdir ./output \
      -o gwas_all
```
注：确保表型文件与 fam 文件中的样本 ID 一致且对应

### 结果可视化
`Manhattan Plot`：展示全基因组范围内所有遗传位点与目标性状的关联显著性，快速定位显著关联的候选位点  
`Q-Q Plot`：评估 GWAS 统计模型的准确性，检查是否存在假阳性  
`LD Heatmap`：展示全基因组范围内遗传变异位点之间的连锁不平衡程度，帮助识别候选基因  
使用 qqplot 绘制曼哈顿图和 Q-Q 图，详见 `Manhattan_QQ_Plot.R` `Manhattan_local_region.R`。
使用 LDBlockShow 绘制 LD 热图:
```bash
# 将关联结果文件转换为 LDBlockShow 可读格式，只保留 chr ps p_wald 信息
awk 'BEGIN {OFS="\t"} NR==1 {next} {
    ch = $1;
    sub(/^chr/, "", ch);
    if (ch ~ /^[0-9]+$/ && ch < 10 && length(ch) == 1) {
        chr_fixed = sprintf("chr0%d", ch);
    } else {
        chr_fixed = "chr" ch;
    }
    print chr_fixed, $3, $NF;
}' gwas_all.assoc.txt > gwas_3_columns.txt
```
```bash
# 如果是 SV-VCF 文件，将 REF 和 ALT 保留为单字符
zcat modified_SV_maf005.vcf.gz | awk '
BEGIN {OFS="\t"} 
/^#/ {print; next} 
{
    if (length($4) > 1) $4 = substr($4, 1, 1);
    if (length($5) > 1) {
        $5 = substr($5, 1, 1);
    }
    print
}' | bgzip -c > modified_SV_maf005_snpstyle.vcf.gz
bcftools index -t modified_SV_maf005_snpstyle.vcf.gz
```
```bash
# 使用 LDBlockShow 绘制 LD 热图
LDBlockShow -InVCF modified_SV_maf005_snpstyle.vcf.gz \
            -OutPut $OUT_PREFIX \
            -Region 10:53300000-54300000 \
            -SeleVar 2 \
            -InGWAS gwas_3_columns.txt \
            -InGFF DM8.1_gene.gff3 \
            -OutPng
```

## 筛选候选基因
### 定位
1. 根据 LD 热图分析显著性最高的位点附近连锁紧密的 SV 位点，进行候选基因定位
```bash
# Step1
awk '$1=="chr10" && $2>=53730000 && $2<=53888000' gwas_3_columns.txt | sort -k3,3g | head -n 5 > top5_pos
# Step2
awk '{print $1 "\t" $2-1 "\t" $2}' top5_pos > regions.bed
bcftools view -R regions.bed modified_SV_maf005.vcf.gz -Oz -o candidate_SV.vcf.gz
bcftools index -t candidate_SV.vcf.gz
rm regions.bed
# Step3
awk '!/^#/ {
    chrom=$1; pos=$2; ref=$4; alt=$5;
    len_ref = length(ref);
    len_alt = length(alt);
    if (len_ref > len_alt) { end = pos + len_ref - 1; } else { end = pos; }
    print chrom "\t" pos-1 "\t" end;
}' <(bcftools view candidate_SV.vcf.gz) > candidate_top5.bed
# Step4
awk '
BEGIN { OFS="\t" }
NR==FNR {
    b_chrom[NR] = $1;
    b_start[NR] = $2;
    b_end[NR] = $3;
    total = NR;
    next
}
!/^#/ && $3=="gene" {
    g_chrom = $1; g_start = $4; g_end = $5; g_attr = $9;
    for (i = 1; i <= total; i++) {
        if (g_chrom == b_chrom[i]) {
            s = (b_start[i] > g_start) ? b_start[i] : g_start;
            e = (b_end[i] < g_end) ? b_end[i] : g_end;
            if (s <= e) {
                print g_attr;
            }
        }
    }
}
' candidate_top5.bed /public/home/zhaoshuo/work1/data/reference/DM8.1_gene.gff3 | sort -u > candidate_gene
# Step5
while read -r line; do
    [ -z "$line" ] && continue
    [[ "$line" =~ ID=([^;]+) ]] && gene_id="${BASH_REMATCH[1]}" || gene_id="$line"

    awk -v t="$gene_id" '!/^#/ && $3=="gene" && $9 ~ t {print $1, $4-1, $5, t; exit}' \
        OFS="\t" /public/home/zhaoshuo/work1/data/reference/DM8.1_gene.gff3 > "${gene_id}.bed"

    if [ -s "${gene_id}.bed" ]; then
        bedtools getfasta -fi /public/home/zhaoshuo/work1/data/reference/DM8.1_genome.fasta -bed "${gene_id}.bed" -fo "${gene_id}.fa" -name
        echo "已成功生成: ${gene_id}.fa"
    else
        echo "警告: 仍未匹配到坐标，请检查 ${gene_id} 是否存在于 GFF3 中"
    fi
    
    rm -f "${gene_id}.bed"
done < candidate_gene
```