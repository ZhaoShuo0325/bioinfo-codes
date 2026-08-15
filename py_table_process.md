<!--
 * @Author: Shuo Zhao && 18904530325@163.com
 * @Date: 2026-08-15 15:56:00
 * @LastEditors: Shuo Zhao && 18904530325@163.com
 * @LastEditTime: 2026-08-15 16:31:10
 * @FilePath: /Code_Notes/py_table_process.md
 * @Description: 
 * 
-->

# Python 处理表格文件

## 安装库
   - `pandas`：用于数据分析、批量处理、行列过滤、统计计算，底层基于Numpy，速度快。  
   - `openpyxl`：用于读写 Excel 文件（.xlsx），支持复杂的 Excel 格式化操作。  

```bash
# 在当前环境安装 pandas 和 openpyxl库
pip install pandas openpyxl
```

## 读取 Excel 文件  
1. 使用 `read_excel()` 函数读取 Excel 文件  
Pandas 使用 openpyxl 等库来读取 Excel 文件  
```python
import pandas as pd

df = pd.read_excel(
    "data.xlsx", # 文件路径
    sheet_name=0 # 0 表示第 1 个 sheet，或 sheet 名称（sheet_name='Sheet1'）
    )
print(df.head()) # 打印前 5 行数据
```
同时读取多个 sheets 将返回一个有序字典  
```python
xls = pd.ExcelFile('data.xlsx')
dict_dfs = pd.read_excel(xls, sheet_name=None) # sheet_name=None 读取所有 sheet
df_sheet1 = dict_dfs['Sheet1'] # 读取指定 sheet
print(df_sheet1.head())
```