<!--
 * @Author: Shuo Zhao && 18904530325@163.com
 * @Date: 2026-09-19 13:22:42
 * @LastEditors: Shuo Zhao && 18904530325@163.com
 * @LastEditTime: 2026-09-20 21:36:52
 * @FilePath: /Code_Notes/QTL_analysis/QTL_analysis.md
 * @Description: 
 * 
-->

# QTL-Analysis  
## 基本概念    
**QTL (Quantitative Trait Locus)：** 即数量性状位点，代表染色体上影响数量性状的某个区段，区段内可能会有一个甚至多个影响数量性状的功能基因。QTL定位指检查分子标记与QTL间的连锁关系，并估算QTL的表型效应。其本质是利用功能基因与分子标记间的连锁和重组，实现对功能基因位置的定位。  

**QTL 分析的流程：**  
   - **构建分离群体：** 亲本的选择要有稳定的表型差异，一般为 F2 群体或 BC1 群体；  
   - **遗传标记检测和表型测定：**  利用分子标记（如 SNP ）确定分离群体的基因型来源；  
   - **连锁分析：** 观察LOD曲线的峰值，进行QTL的初步定位；  
   - **精细定位：** 通过多轮 PCR 验证缩小区间，使用统计分析鉴定候选基因；  

**Bin Map：** 将染色体上紧密连锁、在群体中具有相同重组断点的多个标记合并为一个 Bin 区间；Bin图 就是由这些 Bin 在染色体上排列组合而成的遗传图谱或可视化图。通过滑动窗口等算法，对分子标记进行打包，将原始 SNP 压缩成几千个甚至几百个高质量的 Bin 标记，同时不丢失真实的重组断点信息。  

## 遗传标记检测与预处理   
1. SNP 检测
**详见 GWAS** https://github.com/ZhaoShuo0325/bioinfo-codes/blob/main/GWAS/GWAS.md  
   - 原始数据预处理  
   - 使用 GATK call SNP、联合分型并过滤  
2. 构建亲本变异数据集 (VariantSet)  
   - 分样本合并 SNP VCF 文件  
```bash
# 按 SAMPLE_LIST合并样本
# 1. 只有亲本和 F1 的 SNP VCF 文件 (prefix=AEF1)
# 2. 包含所有样本的，包括亲本、F1 和 F2 群体 (prefix=total)
> $GVCF_LIST
for sample in $(cat $SAMPLE_LIST); do
    gvcf_path="${sample}_${chr}.g.vcf.gz"
    if [ -f "$gvcf_path" ]; then
        echo "-v $gvcf_path" >> "$GVCF_LIST"
done
# 合并
sentieon driver \
    -t 32 \
    -r $REF \
    --algo GVCFtyper \
    $(cat $GVCF_LIST) \
    total_${chr}_raw.vcf.gz

# Extract SNP
bcftools view -v snps -m 2 -M 2 total_${chr}_raw.vcf.gz -Oz -o total_${chr}_snp.vcf.gz
bcftools index -t total_${chr}_snp.vcf.gz
# Hard Filter
gatk VariantFiltration \
    -R $REF \
    -V total_${chr}_snp.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "SNP_HardFilter" \
    -O total_${chr}_SNPHardFilter.vcf.gz
```
   - 按染色体合并，执行标准化操作  
```bash
# 只需对 AEF1 文件进行处理
list_file="files.list"
find "$VCF_DIR" -maxdepth 1 -name "AEF1_chr*_SNPHardFilter.vcf.gz" | sort -V > "$list_file"
# 合并染色体
bcftools concat -f "$list_file" --threads 32 -Oz -o AEF1_SNPHardFilter.vcf.gz
bcftools index -t --threads 32 AEF1_SNPHardFilter.vcf.gz
# 过滤字段
bcftools annotate -x INFO,^FORMAT/GT $OUT -Oz -o AEF1_SNPHardFilter.GT_only.vcf.gz
bcftools index -t $VCF_DIR/AEF1_SNPHardFilter.GT_only.vcf.gz
```
   - 构建双亲标准变异数据集
```bash
# 过滤标准：变异需在双亲中为0/0, 1/1 或 1/1, 0/0，且在 F1 为 0/1
# Usage: python 00_filter_parents.py <in_vcf.gz> <out_vcf>
python 00_filter_parents.py AEF1_SNPHardFilter.GT_only.vcf.gz AEF1_SNP_VariantSet.vcf
bgzip -@ 8 -f AEF1_SNP_VariantSet.vcf
bcftools index -t AEF1_SNP_VariantSet.vcf.gz
```
   - 提取每个 F2 样本的 SNP
```bash
# 分样本循环，每个样本依次按染色体提取，最后 concat
for sample in $(cat $SAMPLE_LIST); do
> $list_file
while read -r chr; do
    bcftools view -s $sample $chrom_vcf -Ou | \
    bcftools annotate -x INFO,^FORMAT/GT -Oz -o $out_vcf
    bcftools index -t $out_vcf
    echo $out_vcf >> $list_file
done < $CHR
bcftools concat -f "$list_file" --threads 32 | \
bcftools view -e 'GT="missing"' --threads 32 -Oz -o "${sample}_SNP.vcf.gz"
bcftools index -t "${sample}_SNP.vcf.gz"
done
```

## 绘制 Bin 图
1. 滑动窗口计算杂合度
```bash
# 以 1Mb 为窗口 100kb 为步长分 Bin，计算杂合度
# Usage: python 01_windows_1M100K_v4.11.py <input_parents_file> <input_snp_file> <output_win_file>
python 01_windows_1M100K_v4.10.py AEF1_SNP_VariantSet.vcf.gz ${sample}_SNP.vcf.gz ${sample}_win.tsv
```
2. 判断 F2 群体中的亲本基因型
```bash
# Usage: python 02_snp_geno.py <parents_vcf> <f2_snp_file> <output_geno_file>
python 02_snp_geno.py AEF1_SNP_VariantSet.vcf.gz ${sample}_SNP.vcf.gz ${sample}_geno.tsv
```
3. 判断 Bin 的亲本基因型
```bash
# Usage: python 03_bin_geno.py <01-win_file> <02-snp_geno_file> <output_detail_file> <output_geno_file>
python 03_bin_geno.py ${sample}_win.tsv ${sample}_geno.tsv ${sample}_bin_geno_details.tsv ${sample}_bin_geno.tsv
```