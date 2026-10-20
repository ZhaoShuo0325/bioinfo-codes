<!--
 * @Author: Shuo Zhao && 18904530325@163.com
 * @Date: 2026-10-19 14:26:45
 * @LastEditors: Shuo Zhao && 18904530325@163.com
 * @LastEditTime: 2026-10-20 20:55:00
 * @FilePath: /Code_Notes/Genome_Assembly/genome_annotation.md
 * @Description: 
 * 
-->

# Genome annotation  
**注释工具：** 使用 `EDTA` 注释重复序列和转座子， `Helixer` 从头预测基因结构， `Eviann` 进行基因结构整合和精细注释  
**Learn From:** Zhigui Bao https://github.com/baozg/assembly-annotation-pipeline  

## 原始数据  
   - `H1/H2.fa` 组装并评估正确的基因组文件  
   - `clean.fq.gz` 过滤后的转录组测序文件  

## 重复序列/转座子注释 EDTA  
   - **`EDTA` 功能：** 预测三类转座子：LTR、TIR 和 Helitron  
   - **结果文件：** 全基因组转座子注释文件(.gff3)、屏蔽后基因组文件(.masked)、非冗余转座子库(.TElib.fa)和统计文件  
推荐使用 snakemake 工作流，可分三步并行提高效率（参考 https://github.com/baozg/assembly-annotation-pipeline）  
```bash
# 修改后的 01.EDTA.smk 支持指定样本，也可以并行多样本
snakemake -s 01.EDTA.smk --config sample="$sample" --cores 32 --rerun-incomplete --nolock -p
```
```bash
# 也可直接运行 EDTA.pl
EDTA.pl --genome ${sample}.H1.fa --species others  --sensitive 1 --anno 1 -t 32
```

## 从头预测基因结构 Helixer
   - **`Helixer` 功能：** 基于深度学习基因注释，无需转录组或同源蛋白证据，直接预测基因位置
   - **结果文件：** 基因注释文件(.gff/gtf)
推荐使用 GPU 运行，可极大提升效率
```bash
# 这里展示 GPU 运行
# --lineage land_plant 参数表明注释物种的类型，确保 ~/.local/share/Helixer/models/land_plant/ 有该文件
singularity run --nv $HELIXER Helixer.py --fasta-path ${sample}.H1.fa  \
--lineage land_plant --gff-output-path ${sample}.H1.gff > ${sample}.H1.gpu.log 2>&1
```

## 基于转录组的基因结构注释 EviAnn
   - **EviAnn 功能：** 基于转录组数据，对蛋白质编码基因、假基因和非编码 RNA 等进行注释  
   - **结果文件：** 基因注释文件(.gff/gtf)  
准备过滤后的双端转录组数据 clean.fq.gz，按如下格式准备文件  
```bash
# 生成转录组位置文件，实例展示的为二代转录组测序的文件格式
rnaseq_list="${sample}_rnaseq.txt"
> $rnaseq_list
for tissue in leaf root stem; do
    t_dir=$RNA_DIR/${sample}_${tissue}
    if [ -d "$t_dir" ]; then
        r1="$t_dir/${sample}_${tissue}_1.clean.fq.gz"
        r2="$t_dir/${sample}_${tissue}_2.clean.fq.gz"
        if [ -f "\$r1" ] && [ -f "\$r2" ]; then
            echo "$r1 $r2 fastq" >> $rnaseq_list
        fi
    fi
done
```
以 Helixer 预测的文件为骨架，EviAnn 能够结合转录组证据进一步完善  
```bash
# 如有两套单倍型，需单独建立文件夹，否则中间文件干扰
# -r 转录组位置文件； -c Helixer 预测文件； -e 参考基因组 CDS.fasta 文件； -p 近缘物种蛋白库 pep.fa 文件
mkdir -p ${sample}/H1
cd ${sample}/H1
eviann.sh -t 32 -g ${sample}.H1.fa -r $rnaseq_list -c $HELIXER_DIR/${sample}/${sample}.H1.gff -e $DM_CDS -p $related_pro
```

## 使用 PASApipeline 对注释文件进行 polish