<!--
 * @Author: Shuo Zhao && 18904530325@163.com
 * @Date: 2026-08-13 09:57:00
 * @LastEditors: Shuo Zhao && 18904530325@163.com
 * @LastEditTime: 2026-08-13 12:46:40
 * @FilePath: /Code_Notes/WGS_analysis.md
 * @Description: 
 * 
-->

# WGS-Data-Analysis-Workfolw

## Fastp
**Fastp** 快速预处理和质控短读长原始 fastq 数据
```bash
# 以双端测序数据为例
fastp -i r1.fq.gz -I r2.fq.gz -o r1.clean.fq.gz -O r2.clean.fq.gz -w 24 --detect_adapter_for_pe
```

## BWA
**BWA (Burrows-Wheeler Aligner)** 用于将二代测序高通量 Reads 快速比对到大型参考基因组，在变异检测、RNA-seq等分析中连接原始数据与下游分析
1. 参考基因组比对并排序
```bash
# 准备参考基因组 ref.fa 确保建立 bwa index
bwa mem -t 16 -R "@RG\tID:${sample}\tSM:${sample}\tPL:illumina" $REF r1.clean.fq.gz r2.clean.fq.gz | samtools view -@ 16 -b - | samtools sort -@ 16 -o ${sample}.bam
samtools index -@ 16 ${sample}.bam
```

2. 比对结果统计（seqkit、samtools、pandepth）
```bash
SAMPLE_TSV="\$SAMPLE_OUT/${sample}_stats.tsv"

    # 1. seqkit stats (计算 Total_Reads 和 Total_Length)
    stats_output=$(seqkit stats -j 24 "$fq1" "$fq2" 2>/dev/null | tail -n +2)
    if [ -n "$stats_output" ]; then
        total_reads=$(echo "$stats_output" | awk '{gsub(/,/, "", \$4); sum+=\$4} END {print sum}')
        total_length=$(echo "$stats_output" | awk '{gsub(/,/, "", \$5); sum+=\$5} END {print sum}')
    fi

    # 2. samtools flagstat (提取 Mapping_Rate)
    if [ -f "$bam_file" ]; then
        flagstat_out=$(samtools flagstat -@ 24 "$bam_file" 2>/dev/null)
        mapping_rate=$(echo "$flagstat_out" | grep "mapped (" | grep -oP '\(\K[0-9.]+(?=%)')
    fi

    # 3. pandepth stats (计算 Avg_Depth)
    out_prefix="$SAMPLE_OUT/pandepth_temp"
    stat_file="${out_prefix}.chr.stat.gz"
    
    pandepth -i "$bam_file" -o "$out_prefix" -t 24 >/dev/null 2>&1
    
    if [ -f "$stat_file" ]; then
        avg_depth=$(zcat "$stat_file" 2>/dev/null | tail -n 1 | grep -oP 'MeanDepth:\s*\K[0-9.]+')
    fi
    
    # 4. 根据 Total_Length 和 Mapping_Rate 计算 Mapping_Base
    if [ "$total_length" != "N/A" ] && [ "$mapping_rate" != "N/A" ]; then
        mapping_base=$(awk -v len="$total_length" -v rate="$mapping_rate" 'BEGIN {printf "%.0f", len * rate / 100}')
    fi

    echo -e "${sample}\t\$total_reads\t\$total_length\t\$mapping_base\t\${mapping_rate}%\t\$avg_depth" > "\$SAMPLE_TSV"
```

3. 标记并去除 PCR 重复
```bash
picard MarkDuplicates \
    I=${sample}.bam \
    O=${sample}_rmdup.bam \
    M=${sample}_dup_metrics.txt \
    REMOVE_DUPLICATES=true
samtools index -@ 16 ${sample}_rmdup.bam
```
