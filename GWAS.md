<!--
 * @Author: Shuo Zhao && 18904530325@163.com
 * @Date: 2026-08-12 10:37:32
 * @LastEditors: Shuo Zhao && 18904530325@163.com
 * @LastEditTime: 2026-08-14 16:01:52
 * @FilePath: /Code_Notes/GWAS.md
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
#运行 GLnexus
singularity exec glnexus-1.4.1.sif glnexus_cli \
  --config DeepVariant \
  --dir $DB_DIR \
  --threads 32 \
  --mem-gbytes 120 \
  --list $SAMPLE_FILE > "total_merged.bcf"
```

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
# 去除非 GT 字段 排序 建立索引
bcftools annotate -x INFO,^FORMAT/GT survivor_merged.vcf.gz -Oz -o survivor_merged.GT_only.vcf.gz
bcftools sort survivor_merged.GT_only.vcf.gz -Oz -o survivor_merged_sorted.vcf.gz
bcftools index -t survivor_merged_sorted.vcf.gz
# 将合并后缺失的 GT 添为 0/0
zcat survivor_merged_sorted.vcf.gz | awk '{if($0 ~ /^#/) print; else {gsub(/\.\/\./, "0/0"); print}}' | bgzip > modified_SV.vcf.gz
```

## 群体结构分析
### PCA 分析
**PLINK** 可以将 VCF 文件转换为 BED BIM FAM 等格式，用于处理大型数据集，进行 PCA 分析等
```bash
# 将 VCF 文件转换为 BED BIM FAM 格式，指定输出文件前缀
# --set-missing-var-ids @:# 给没有ID的位点自动命名
plink --vcf modified_SV.vcf.gz --make-bed --out modified_SV --set-missing-var-ids @:#
```
```bash
# LD 修剪，去除高度相关的位点，减少连锁导致的过多权重
# --indep-pairwise 50kb 1 0.2 ：窗口大小 50kb，步长 1个变异位点，相关性阈值 0.2
plink --bfile modified_SV --indep-pairwise 50kb 1 0.2 --out modified_SV_pruned
plink --bfile modified_SV --extract modified_SV_pruned.prune.in --make-bed --out modified_SV_pruned
```
```bash
# 进行 PCA 分析，指定输出文件前缀和主成分数量
plink --bfile modified_SV_pruned --pca 10 --out modified_SV_pca
```