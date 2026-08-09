<!--
 * @Author: Shuo Zhao && 18904530325@163.com
 * @Date: 2026-08-09 09:59:18
 * @LastEditors: Shuo Zhao && 18904530325@163.com
 * @LastEditTime: 2026-08-09 11:43:43
 * @FilePath: /Code_Notes/BCFTools.md
 * @Description: 
 * 
-->

# BCFTools

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
   - -r --regions：指定区域 (chr01:1000-2000)

4. **指定样本 `-S/-s`**
   - `-S` --samples-file [^]LIST：提取指定样本（[^]排除指定样本）
   - -s --samples [^]FILE：提取指定样本（[^]排除指定样本）

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