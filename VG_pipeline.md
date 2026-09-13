<!--
 * @Author: Shuo Zhao && 18904530325@163.com
 * @Date: 2026-08-18 10:34:48
 * @LastEditors: Shuo Zhao && 18904530325@163.com
 * @LastEditTime: 2026-09-13 19:54:59
 * @FilePath: /Code_Notes/VG_pipeline.md
 * @Description: 
 * 
-->

# VG (Variation Graph Toolkit)  
**vg** 是用于构建图泛基因组，并在图中直接进行高精度序列比对、变异检测和基因分型的工具。https://github.com/vgteam/vg  

## vg 的文件格式
   - `.vg`：vg 的原生格式，存储图本身的拓扑结构及序列，在大规模下效率可能较低，尤其是当包含大量路径时。  
   - `.xg`：索引化的静态图格式，用于快速进行节点查询和路径遍历。  
   - `.gfa`：GFA 是一种标准的文本格式，可在 vg 与其他泛基因组工具之间交换图结构，可编辑。  
   - `.gbz`：GBZ 是一种高度压缩的格式，与上述格式相比，它存储路径所需的空间更小，但不允许对图进行常规编辑。  
   - `.dist`：用来存储图上节点间距离与路径长度分布的二进制索引文件。  
   - `.snarls`：用来存储变异图中所有 Snarl 结构及其层级关系的二进制文件。  
   - `.min`：用于储存基于 snarl 树结构计算图上两点间的最短路径的文件。  
   - `.zipcodes`：利用 snarl 树的层级结构为图中的节点分配坐标编码，能够高效地查询图上任意两个位点之间的距离，最新。  
   - `.gam`：GAM 是 vg 专有的二进制比对格式，类似于.sam，用于记录 reads 在变异图中的比对路径。  
   - `.pack`：是 vg 专有的二进制压缩文件，用于记录测序 reads 图中各个节点和边上的覆盖深度，以及比对质量权重。  

## VCF 文件的准备  
1. VCF 文件对齐和排序
```bash
for sample in $(cat $SAMPLE_LIST); do
    bcftools norm --threads 32 -m -any -f $REF ${sample}.vcf.gz -Oz | bcftools sort -Oz -o ${sample}.norm.vcf.gz
done
```
2. 拆分染色体
```bash
for chr in $(cat ${CHR}); do
    bcftools view ${sample}.norm.vcf.gz ${chr} -Oz -o ${sample}.${chr}.vcf.gz
    tabix -p vcf ${sample}.${chr}.vcf.gz
done
```
3. 按染色体合并样本
```bash
for chr in $(cat $CHR); do
    ls *.${chr}.vcf.gz > ${chr}_vcf.list
    MERGED_VCF=${chr}_merged.vcf.gz
    bcftools merge -l ${chr}_vcf.list | \
    bcftools norm -m -any -N | \
    bcftools norm -d none -f $REF | \
    bcftools sort | \
    bgzip > $MERGED_VCF
    tabix -p vcf $MERGED_VCF
```
4. vg 去冗余  
vg 有自己专属的去冗余方法，详见：https://github.com/vgteam/giraffe-sv-paper/tree/master/scripts/sv/remap-to-dedup-merged-svs  
**注**：去冗余标准为 1.相同类型的变异，2.重叠率大于 80%，3.插入位点间距在 50bp 范围内。符合上述要求的 SV 属于同一个冗余簇 neardups-clusters.tsv。然后通过 Re-mapping 的方式在冗余簇中选择有代表性的 SV。  

5. 过滤群体频率
```bash
# 将所有的 ./. 改写为 0/0
zcat DupFilter_${chr}.vcf.gz | sed 's/\\.\\/\\./0\\/0/g' | bgzip > Dupfilter_${chr}_raw.vcf.gz
tabix -p vcf Dupfilter_${chr}_raw.vcf.gz
# 过滤掉频率小于 0.01 的变异
vcftools --gzvcf Dupfilter_${chr}_raw.vcf.gz --maf 0.05 --recode --recode-INFO-all --stdout | bgzip > Dupfilter_${chr}.maf001.vcf.gz
tabix -p vcf Dupfilter_${chr}.maf001.vcf.gz
# 将所有染色体 concat 到一个文件
bcftools concat --threads 16 Dupfilter_${chr}.maf001.vcf.gz* -Oz -o all_maf001.vcf.gz
bcftools index -t all_maf001.vcf.gz
```

## 构建图泛基因组  
1. 构建图泛基因组  
```bash
# vg construct 构建变异图
vg construct -t 32 \
    -r $REF \
    -v all_maf001.vcf.gz \
    -a -p > ${PREFIX}.vg
```
2. 为图建立索引  
```bash
# 建立 gbwt 索引
vg gbwt --num-jobs 32 \
    -x ${PREFIX}.vg \
    -v all_maf001.vcf.gz \
    -o ${PREFIX}.gbwt \
    --vcf-variants \
    --progress

# 生成 xg 文件
vg index -t 32 -x ${PREFIX}.xg ${PREFIX}.vg

# 为短读长 (sr-giraffe) 建立所需的全部索引文件
vg autoindex \
    --workflow sr-giraffe \
    --prefix ${PREFIX} \
    --ref-fasta $REF \
    --vcf all_maf001.vcf.gz \
    --threads 64 \
    --target-mem 200G

# 生成 snarls 文件
vg snarls -t 32 ${PREFIX}.gbz > ${PREFIX}.snarls
```

## 基于图泛基因组的 genotyping  
1. vg giraffe
```bash
# 检查 fq 文件是否质控
for sample in $(cat $SAMPLE_LIST); do
F1=$FQ_DIR/${sample}/*_1.clean.fq.gz
F2=$FQ_DIR/${sample}/*_2.clean.fq.gz

vg giraffe -t 64 \
        -Z $GBZ \
        -d $DIST \
        -m $MIN \
        -z $ZIP \
        -f $F1 \
        -f $F2 \
        -o gam > \${sample}_giraffe.gam
done
```
2. vg pack  
```bash
for sample in $(cat $SAMPLE_LIST); do
vg pack -t 32 -Q 5 \
        -x $GBZ \
        -g $GAM \
        -o "${sample}.pack
done
```
3. vg call
```bash
for sample in $(cat $SAMPLE_LIST); do
REF_PATHS=$(printf -- '-p chr%02d ' {1..12})
vg call -t 32 $GBZ \
        -r $SNARLS \
        -k $PACK \
        -s ${sample} \
        $REF_PATHS > ${sample}_raw.vcf
bgzip -@ 16 ${sample}_raw.vcf
done
```

## 结果文件过滤与合并
1. 过滤
```bash
# 只保留 PASS 的变异
bcftools view -f PASS ${sample}_raw.vcf.gz -Oz -o ${sample}_pass.vcf.gz
tabix -f -p vcf ${sample}_pass.vcf.gz
# 只保留 QUAL 大于 30 的变异
bcftools view -i 'QUAL >= 30' ${sample}_raw.vcf.gz -Oz -o ${sample}_q30.vcf.gz
tabix -f -p vcf ${sample}_q30.vcf.gz
```
2. 合并去冗余
```bash
# SURVIVOR 只接受未压缩文件
gunzip -k ${sample}_q30.vcf.gz
# 合并
> "$file_list"
while read -r sample; do
    [ -z "$sample" ] && continue
    echo "${sample}_q30.vcf" >> "$file_list"
done < "$SAMPLE_LIST"
SURVIVOR merge "$file_list" 1000 1 1 1 0 50 "${PREFIX}.vcf"
```