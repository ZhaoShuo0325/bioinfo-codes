<!--
 * @Author: Shuo Zhao && 18904530325@163.com
 * @Date: 2026-08-09 09:59:18
 * @LastEditors: Shuo Zhao && 18904530325@163.com
 * @LastEditTime: 2026-08-11 09:42:31
 * @FilePath: /Code_Notes/BCFTools.md
 * @Description: 
 * 
-->

# BCFTools
**BCFTools** 是用于查看、处理、分析变异文件（VCF 和 BCF）的工具。https://github.com/samtools/bcftools

## VCF 文件结构
1. **基本结构**
   - `头部信息`：以 `##` 开头，记录参考基因组、软件命令及字段定义
   - `列名信息`：以 `#CHROM` 开头，定义每一列含义，包括 `CHROM`、`POS`、`ID`、`REF`、`ALT`、`QUAL`、`FILTER`、`INFO`、`FORMAT`、`SAMPLE`
   - `变异信息`：记录具体的变异位点与各样本的基因型数据

2. **FORMAT字段**
   - `GT` Genotype：基因型信息，0/0，1/0，0/1，1/1，./.等，/ 表示未分型，| 表示已分型
   - `GQ` Genotype Quality：基因型质量
   - `DP` Read Depth：深度信息，表示该位点的测序深度
   - `AD` Allele Depth：支持各等位基因的Reads数量 (e.g. 10,2 表示支持REF的Reads为10，支持ALT的为2)
   - `AC` Allele Count：等位基因个数

## bcftools view
`bcftools view [options] <in.vcf.gz> [region1 [...]]`
1. **输出类型参数 `-O`**
   - `-Oz`：输出compressed.vcf.gz文件
   - `-Ov`：输出uncompressed.vcf文件
   - `-Ou`：用于管道传输

2. **头部信息参数 `-H/-h`**
   - `-h` --header-only：输出文件的头部信息
   - `-H` --no-header：隐藏头部信息

3. **指定区域 `-R/-r`**
   - `-R` --regions-file FILE：指定区域文件
   - `-r` --regions：指定区域 (chr01:1000-2000)

4. **指定样本 `-S/-s`**
   - `-S` --samples-file [^]LIST：提取指定样本（[^]排除指定样本）
   - `-s` --samples [^]FILE：提取指定样本（[^]排除指定样本）

5. **按变异类型筛选 `-V/-v`**
   - `-V` --exclude-types LIST：排除指定变异类型 (e.g. snps,indels)
   - `-v` --types LIST：选择指定变异类型 (e.g. snps,indels)

6. **根据 QUAL DP INFO FORMAT 等字段筛选 `-i/-e`**
   - `-i` --include EXPR：选择指定类型 (e.g. 'QUAL >= 30 && DP > 10')
   - `-e` --exclude EXPR：排除指定类型 (e.g. 'QUAL < 30')

```bash
# 输出sample1 sample2的chr01:1000-2000区域的snps变异，且QUAL >= 30 && DP > 10的行
bcftools view -H -r chr01:1000-2000 -s sample1 sample2 -v snps -i 'QUAL >= 30 && DP > 10' input.vcf.gz
# 过滤QUAL >= 30的行并输出到output.vcf.gz文件
bcftools view -i 'QUAL >= 30' input.vcf.gz -Oz -o output.vcf.gz
```

## bcftools index
`bcftools index [options] <in.bcf>|<in.vcf.gz>`
1. **建立csi/tbi索引 `-c/-t`**
   - `-c` --csi：建立csi索引
   - `-t` --tbi：建立tbi索引

2. **统计信息 `-s/-n`**
   - `-s` --stats：统计变异信息 (#Contig #Length #VariantCount)
   - `-n` --nrecords：统计总变异数 (Better than `bcftools view -H input.vcf.gz | wc -l`)

## bcftools norm
`bcftools norm [options] <in.vcf.gz>`
1. **参考基因组 `-f`**
   - `-f` --fasta-ref FILE：指定参考基因组文件

2. **多等位基因位点处理 `-m`**
   - `-m -`：把一行中多等位基因拆成多行
   - `-m +`：把多行相同位置的多等位基因合并成一行
   - `-m -snps|indels|both|any`：指定处理多等位基因的类型

3. **参考碱基校验与修复 `-c`**
   - `-c` --check-ref e：报错并退出 exit (default)
   - `-c` --check-ref w：警告并继续 warn
   - `-c` --check-ref x：剔除并继续 exclude
   - `-c` --check-ref s：强行修正为参考基因组的碱基 set

## bcftools sort
`bcftools sort [OPTIONS] <FILE.vcf>`
   - `--write-index`：自动为输出文件生成索引

```bash
# 对input.vcf.gz文件拆分多等位基因，并进行排序，输出到output.vcf.gz，自动建立索引
bcftools norm --threads 32 -m -any -f ref.fa input.vcf.gz -Ou | \
bcftools sort --write-index -Oz -o output.vcf.gz
```

## bcftools merge
`bcftools merge [options] <A.vcf.gz> <B.vcf.gz> [...]`
1. **等位基因合并处理 `-0/-m`**
   - `-0` --missing-to-ref：把./.补成0/0
   - `-m` --merge STRING：相同位置等位基因合并策略 (e.g. snps, indels, both)

2. **输入文件管理 `-l/--force-sample`**
   - `-l` --file-list FILE：指定输入文件路径列表
   - `--force-sample`：强制合并所有样本，给重复的样本名重新编号

3. **INFO字段处理方法 `-i`**
   - `-i` --info-rules TAG:METHOD：标签名:处理方法,标签名:处理方法 (e.g. DP:sum,AF:avg) 默认:sum

```bash
# 合并sample1.vcf.gz和sample2.vcf.gz文件，DP求和，AF求平均值，输出到merged.vcf.gz文件
bcftools merge -i DP:sum,AF:avg sample1.vcf.gz sample2.vcf.gz -Oz -o merged.vcf.gz
# 合并sample_list.txt中所有文件，强制合并所有样本，输出到merged.vcf.gz文件，并自动建立索引
bcftools merge -l sample_list.txt --force-sample --write-index -Oz -o merged.vcf.gz
```