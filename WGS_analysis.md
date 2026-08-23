<!--
 * @Author: Shuo Zhao && 18904530325@163.com
 * @Date: 2026-08-13 09:57:00
 * @LastEditors: Shuo Zhao && 18904530325@163.com
 * @LastEditTime: 2026-08-19 14:36:14
 * @FilePath: /Code_Notes/WGS_analysis.md
 * @Description: 
 * 
-->

# WGS-Data-Analysis-Workfolw

## 原始数据质控
**Fastp** 快速预处理和质控短读长原始 fastq 数据
```bash
# 以双端测序数据为例
fastp -i r1.fq.gz -I r2.fq.gz -o r1.clean.fq.gz -O r2.clean.fq.gz -w 24 --detect_adapter_for_pe
```

## 比对与预处理
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
| Sample | Total_Reads | Total_Length(bp) | Mapping_Base | Mapping_Rate(%) | Avg_Depth |
| :--- | :---: | :---: | :---: | :---: | :---: |
| Sample1 | 92000952 | 13784281238 | 13687791269 | 99.30% | 18.14 |
| Sample2 | 93512932 | 14010639164 | 13906960434 | 99.26% | 18.34 |
| Sample3 | 96113584 | 14313934010 | 14213736472 | 99.30% | 18.92 |
| Sample4 | 110854190 | 16528258961 | 16448923318 | 99.52% | 21.88 |
| Sample5 | ... | ... | ... | ... | ... |

3. 标记并去除 PCR 重复  
**Picard** 用于标记或去除由于 PCR 扩增产生的重复序列，减少假阳性变异
```bash
# 不建议直接使用picard MarkDuplicates，脚本默认内存 -Xmx1g 可能不够
java -Xmx32g -jar /public/home/zhaoshuo/miniconda3/envs/bio_env/share/picard-2.20.4-0/picard.jar MarkDuplicates \
    I=${sample}.bam \
    O=${sample}_rmdup.bam \
    M=${sample}_dup_metrics.txt \
    REMOVE_DUPLICATES=true
samtools index -@ 16 ${sample}_rmdup.bam
```
```bash
# 生成质量评估报告
java -Xmx32g -jar /public/home/zhaoshuo/miniconda3/envs/bio_env/share/picard-2.20.4-0/picard.jar CollectWgsMetrics \
    I=${sample}_rmdup.bam \
    O=${sample}_metrics.txt \
    R=$REF
```

## Q & A
1. **Q：** 运行 Picard 时发生报错 `Value was put into PairInfoMap more than once`  
**A：** 是因为多个测序 reads 有相同的名字，可利用 samtools 强硬去除
```bash
# -f 0x2：丢弃单端比对（Single-end）、比对质量极差、未正确配对或被怀疑是二次比对的 reads
samtools view -@ 16 -b -f 0x2 input.bam > output.bam
```

2. **Q：** 运行完 Picard 去重后，bam 文件大幅缩水  
**A：** 检查 fastp 日志文件，查看是否是测序质量的原因  
```bash
Filtering result:
reads passed filter: 435652476
reads failed due to low quality: 172308
reads failed due to too many N: 0
reads failed due to too short: 280
reads with adapter trimmed: 129867062
bases trimmed due to adapters: 7067296008

Duplication rate: 56.3309%
```
或使用 `samtools flagstat` 查看生成的 bam 文件质量  
```bash
229273259 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
5030601 + 0 supplementary
0 + 0 duplicates
227812497 + 0 mapped (99.36% : N/A)
224242658 + 0 paired in sequencing
112121329 + 0 read1
112121329 + 0 read2
71848558 + 0 properly paired (32.04% : N/A)
222481395 + 0 with itself and mate mapped
300501 + 0 singletons (0.13% : N/A)
14843027 + 0 with mate mapped to a different chr
4881004 + 0 with mate mapped to a different chr (mapQ>=5)
```