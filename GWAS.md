<!--
 * @Author: Shuo Zhao && 18904530325@163.com
 * @Date: 2026-08-12 10:37:32
 * @LastEditors: Shuo Zhao && 18904530325@163.com
 * @LastEditTime: 2026-08-13 09:56:38
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
#### DeepVariant
1. 准备 bam 文件
**确保双端 fq 文件经过过滤**
```bash
bwa mem -t 16 -R "@RG\tID:${sample}\tSM:${sample}\tPL:illumina" $REF $F1 $F2 | samtools view -@ 16 -b - | samtools sort -@ 16 -o ${sample}.bam
samtools index -@ 16 ${sample}.bam
```

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
