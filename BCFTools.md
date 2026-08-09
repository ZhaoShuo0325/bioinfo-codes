<!--
 * @Author: Shuo Zhao && 18904530325@163.com
 * @Date: 2026-08-09 09:59:18
 * @LastEditors: Shuo Zhao && 18904530325@163.com
 * @LastEditTime: 2026-08-09 11:08:30
 * @FilePath: /Code_Notes/BCFTools.md
 * @Description: 
 * 
-->

# BCFTools

## BCFTools view
`bcftools view [options] <in.vcf.gz> [region1 [...]]`
1. **输出类型参数 `-O`**
   - `-Oz`：输出compressed.vcf.gz文件
   - `-Ov`：输出uncompressed.vcf文件
   - `-Ou`：用于管道传输

2. **头部信息参数 `-H/-h`**
   - `-h` --header-only：输出文件的头部信息
   - `-H` --no-header：隐藏头部信息

3. **指定区域 `-R/-r`**
   - -R --regions-file FILE：指定区域文件
   - -r --regions：指定区域 (chr01:1000-2000)

4. **指定样本 `-S/-s`**
   - -S --samples-file [^]LIST：提取指定样本（[^]排除指定样本）
   - -s --samples [^]FILE：提取指定样本（[^]排除指定样本）

5. **按变异类型筛选 `-V/-v`**
   - -V --exclude-types LIST：排除指定变异类型 (e.g. snps,indels)
   - -v --types LIST：选择指定变异类型 (e.g. snps,indels)

6. **根据QUAL DP INFO FORMAT等字段筛选 `-i/-e`**
   - -i --include EXPR：选择指定类型 (e.g. 'QUAL >= 30 && DP > 10')
   - -e --exclude EXPR：排除指定类型 (e.g. 'QUAL < 30')

```bash
# 输出sample1 sample2的chr01:1000-2000区域的snps变异，且QUAL >= 30 && DP > 10的行
bcftools view -H -r chr01:1000-2000 -s sample1 sample2 -v snps -i 'QUAL >= 30 && DP > 10' input.vcf.gz
# 过滤QUAL >= 30的行并输出到output.vcf.gz文件
bcftools view -i 'QUAL >= 30' input.vcf.gz -Oz -o output.vcf.gz
```
