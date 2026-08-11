<!--
 * @Author: Shuo Zhao && 18904530325@163.com
 * @Date: 2026-08-10 14:37:43
 * @LastEditors: Shuo Zhao && 18904530325@163.com
 * @LastEditTime: 2026-08-11 09:42:37
 * @FilePath: /Code_Notes/SURVIVOR.md
 * @Description: 
 * 
-->

# SURVIVOR
**SURVIVOR** 适用于处理、比较和模拟结构变异的工具，能够将多个样本的SV进行合并去冗余。https://github.com/fritzsedlazeck/SURVIVOR

## SURVIVOR merge
1. `File with VCF names and path`：输入多个VCF文件的路径和名称 (vcf_list.txt) SURVIVOR不支持输入压缩格式！
2. `max distance between breakpoints`：断点最大偏差距离，如相差1000bp内视为同一个SV
3. `Minimum number of supporting caller`：支持该变异的最少软件数，如至少2个软件支持才保留
4. `Take the type into account`：是否要求同一类型才合并，1必须类型一致，0忽略类型
5. `Take the strands of SVs into account`：是否要求变异链方向相同才合并，1必须方向一致，0忽略方向
6. `Estimate distance based on the size of SV`：是否根据SV的实际大小动态估算允许的断点距离，1动态估算，0固定距离
7. `Minimum size of SVs to be taken into account`：最小SV长度，一般设置为50
8. `Output VCF filename`：输出文件路径

```bash
# 使用实例
ls cuteSV.vcf Sniffles.vcf pbsv.vcf > vcf_list.txt
SURVIVOR merge vcf_list.txt 1000 1 1 1 0 50 survivor_output.vcf
bgzip survivor_output.vcf
bcftools sort survivor_output.vcf.gz -Oz -o survivor_output_sorted.vcf.gz
bcftools index -t survivor_output_sorted.vcf.gz
```
