<!--
 * @Author: Shuo Zhao && 18904530325@163.com
 * @Date: 2026-08-09 09:59:18
 * @LastEditors: Shuo Zhao && 18904530325@163.com
 * @LastEditTime: 2026-08-09 10:17:39
 * @FilePath: /Code_Notes/BCFTools.md
 * @Description: 
 * 
-->

# BCFTools

## BCFTools view
**bcftools view [options] <in.vcf.gz> [region1 [...]]**
1. 输出类型参数 -O
-Oz：输出compressed.vcf.gz文件
-Ov：输出uncompressed.vcf文件
-Ou：用于管道传输

2. 头部信息参数 -H/-h
-h --header-only：输出文件的头部信息
-H --no-header：隐藏头部信息

3. 指定区域 -R/-r
-R --regions-file：指定区域文件
-r --regions：指定区域 (chr01:1000-2000)