<!--
 * @Author: Shuo Zhao && 18904530325@163.com
 * @Date: 2026-09-23 09:24:06
 * @LastEditors: Shuo Zhao && 18904530325@163.com
 * @LastEditTime: 2026-09-27 22:42:54
 * @FilePath: /Code_Notes/Genome_Assembly/genome_assembly.md
 * @Description: 
 * 
-->

# Genome assembly  
**组装方法：** 使用 `hifiasm` 进行 HiFi 与 Hi-C 联合分型组装 (Haplotype-resolved)，适合高杂合的二倍体马铃薯。

## 原始数据质控  
   - `PacBio HiFi` .bam file -> .fq file  
   - `HiC` .fq file  

**1. HiFi.bam 文件转换为 .fq 格式**
```bash
samtools fastq -@ 16 $bam_file | bgzip -@ 8 > ${sample}_hifi.fq.gz
```
**2. HiC.fq 文件质控**
```bash
# 过滤参数：
# -q 20 ：单碱基质量值小于 20 为低质量
# -u 50 : 低质量碱基的比例超过 50%，则这条 Read 会被整体丢弃
# -n 5 : Read 中包含超过 5 个无法识别的 N 碱基, 则这条 Read 会被整体丢弃
# -l 150 : 只保留长度为 150 bp 的 Read （有接头即会被丢弃，测序足够深时使用）
# --adapter_sequence : 接头序列，未知可选参数 --detect_adapter_for_pe 自动检测接头
fastp -i $F1 -I $F2 \
     -o ${sample}_r1.clean.fq.gz -O ${sample}_r2.clean.fq.gz \
     --thread 16 \
     --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
     --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
     -q 20 -u 50 -n 5 -l 150
```
**3. seqkit 原始数据统计**
```bash
seqkit stats -a -T ${sample}_hifi.fq.gz > ${sample}_QC_report.tsv
seqkit stats -a -T ${sample}_r1.clean.fq.gz | tail -n +2 >> ${sample}_QC_report.tsv
seqkit stats -a -T ${sample}_r2.clean.fq.gz | tail -n +2 >> ${sample}_QC_report.tsv
```

## 使用 hifiasm 进行组装  
**hifiasm:** 能够将 HiFi Reads 拼接成较长的 Contigs，并自带纠错和去冗余功能。在输入 Hi-C 数据时，可以直接在图构建阶段区分两条单倍型。https://github.com/chhylp123/hifiasm   
**1. 运行 hifiasm 分型组装**  
```bash
# hifiasm
hifiasm -o ${sample}.assembly -t 32 $hifi_file --h1 $hic_r1 --h2 $hic_r2
```
hifisam 根据 Hi-C 数据，初步拼接两条单倍型 contig 序列，生成 gfa 格式的图文件。  
```bash
# 将 gfa 文件转化为 fa 格式
awk '/^S/{print ">"\$2;print \$3}' ${sample}.assembly.hic.hap1.p_ctg.gfa > ${sample}.assembly.hifiasm.H1.fa
awk '/^S/{print ">"\$2;print \$3}' ${sample}.assembly.hic.hap2.p_ctg.gfa > ${sample}.assembly.hifiasm.H2.fa
# 将两条单倍型合并为一个文件
cat ${sample}.assembly.hifiasm.H1.fa ${sample}.assembly.hifiasm.H2.fa > ${sample}.assembly.Hapall.fa
```
**2. Juicer Hi-C数据比对与互作信息**
准备  
```bash
# Juicer 要求工作目录下要有 fastq 文件，和参考基因组
mkdir -p fastq
ln -sf $HIC_DIR/${sample}/${sample}_r1.clean.fq.gz fastq/${sample}_R1.fastq.gz
ln -sf $HIC_DIR/${sample}/${sample}_r2.clean.fq.gz fastq/${sample}_R2.fastq.gz
# 以 hifiasm 组装的 Hapall 作为参考
ln -sf $ASM_DIR/${sample}/${sample}.assembly.Hapall.fa ${sample}.fa
bwa index ${sample}.fa
```
识别酶切位点，提取每个 contig 长度信息
```bash
# 使用 Juicer 软件中脚本识别酶切位点（一般为 MboI）
python2 $JUICER/misc/generate_site_positions.py MboI ${sample} ${sample}.fa
awk 'BEGIN{OFS="\t"}{print \$1,\$NF}' ${sample}_MboI.txt > ${sample}.chrom.sizes
```
生成 Hi-C 互作比对文件
```bash
sh $JUICER/CPU/juicer.sh \
    -t 32 \
    -s MboI \
    -g ${sample} \
    -d ${sample} \
    -D $JUICER/CPU/common \
    -z ${sample}.fa \
    -p ${sample}.chrom.sizes \
    -y ${sample}_MboI.txt
```
互作结果文件存放在 aligned 文件夹下 merged_nodups.txt，记录 Hi-C 与基因组的比对信息。  
**3. Ragtag 参考基因组辅助的染色体挂载**  
```bash
# 利用 Ragtag 将组装的两套单倍型分别挂载到参考基因组（DM8）上
ragtag.py scaffold -t 32 $REF ${sample}.assembly.hifiasm.H1.fa -o ${sample}/ragtag_H1
ragtag.py scaffold -t 32 $REF ${sample}.assembly.hifiasm.H2.fa -o ${sample}/ragtag_H2
```
从 contig 提升至 scaffold  
```bash
# 使用 Juicer 软件自带脚本将 agp 文件转换为 assembly 格式
python2 $JUICER/juicebox_scripts/juicebox_scripts/agp2assembly.py ${sample}/ragtag_H1/ragtag.scaffold.agp ${sample}/ragtag_H1/${sample}.H1.assembly
python2 $JUICER/juicebox_scripts/juicebox_scripts/agp2assembly.py ${sample}/ragtag_H2/ragtag.scaffold.agp ${sample}/ragtag_H2/${sample}.H2.assembly
```
```bash
# 将 assembly1 assembly2 合并排序
python 04_assembly_merge_sort.py ${sample}/ragtag_H1/${sample}.H1.assembly ${sample}/ragtag_H2/${sample}.H2.assembly ${sample}/${sample}.assembly.Hapall.assembly
```

## 3d-dna 染色体挂载
```bash
# 使用 3d-dna 软件中脚本实现可视化 hic 接触图
# input1: 合并排序后 assembly 文件； input2: juicer 生成的互作文件； output: hic 接触热图，后续在 Juicebox 手工调整
$threeDDNA/visualize/run-assembly-visualizer.sh -q 0 ${sample}/${sample}.assembly.Hapall.assembly $juicer_dir/${sample}/aligned/merged_nodups.txt
```
   - juicebox 需要 .hic 文件和 ${sample}.assembly.Hapall.assembly 文件  
   - 具体操作方法见 https://www.youtube.com/watch?v=Nj7RhQZHM18&t=378s