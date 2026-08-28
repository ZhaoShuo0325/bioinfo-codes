#!/usr/bin/env Rscript
# ==============================================================================
# Script Name: plot_local_region.R
# Description: 用于绘制 GWAS 局部区域关联图（局部曼哈顿图）的 R 脚本。
# 
# Usage:
#   Rscript plot_local_region.R [options]
#
# Options:
#   -i, --input <file>     输入 GWAS 关联分析结果文件路径 
#                          [默认: ../03_GWAS/output4/gwas_all.assoc.txt]
#   -o, --output <file>    输出局部图的路径/名称 
#                          [默认: ../03_GWAS/output4/Local_Manhattan_Chr10.png]
#   -c, --chr <num>        目标染色体编号 [默认: 10]
#   -s, --start <num>      起始物理位置 (bp) [默认: 53300000]
#   -e, --end <num>        终止物理位置 (bp) [默认: 54300000]
#   -a, --alpha <num>      显著性分子 (0.05 或 0.01) [默认: 0.05]
#   -h, --help             显示此帮助信息并退出
# ==============================================================================

library(ggplot2)
library(optparse)

# 1. 定义命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "../03_GWAS/output4/gwas_all.assoc.txt",
              help = "输入 GWAS 关联分析结果文件路径 [默认: %default]", metavar = "file"),
  make_option(c("-o", "--output"), type = "character", default = "../03_GWAS/output4/Local_Manhattan_Chr10.png",
              help = "输出局部图的路径 [默认: %default]", metavar = "file"),
  make_option(c("-c", "--chr"), type = "numeric", default = 10,
              help = "目标染色体 [默认: %default]", metavar = "number"),
  make_option(c("-s", "--start"), type = "numeric", default = 53300000,
              help = "起始位置 (bp) [默认: %default]", metavar = "number"),
  make_option(c("-e", "--end"), type = "numeric", default = 54300000,
              help = "终止位置 (bp) [默认: %default]", metavar = "number"),
  make_option(c("-a", "--alpha"), type = "numeric", default = 0.05,
              help = "显著性阈值 (0.05 或 0.01) [默认: %default]", metavar = "number")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 2. 读取并清洗数据
cat("正在读取数据:", opt$input, "\n")
gwas_res <- read.table(opt$input, header = TRUE, stringsAsFactors = FALSE)
gwas_res <- gwas_res[!is.na(gwas_res$p_wald), ]

n_snps <- nrow(gwas_res)
sig_threshold <- opt$alpha / n_snps

# 3. 提取目标区间子集
cat(sprintf("正在提取 Chr %s : %d - %d bp 的位点...\n", opt$chr, opt$start, opt$end))
local_res <- subset(gwas_res, chr == opt$chr & ps >= opt$start & ps <= opt$end)

# 转换物理位置为 Mb 单位
local_res$Mb <- local_res$ps / 1e6
local_res$neg_log_p <- -log10(local_res$p_wald)

# 4. 使用 ggplot2 绘制局部区域散点图
cat("正在生成局部区域图...\n")
p <- ggplot(local_res, aes(x = Mb, y = neg_log_p)) +
  geom_point(color = "#2A4494", size = 1.8, alpha = 0.8) +
  geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed", color = "black", linewidth = 0.8) +
  labs(
    title = paste0("Regional Association Plot (Chr ", opt$chr, ": ", opt$start/1e6, " - ", opt$end/1e6, " Mb)"),
    x = paste0("Chromosome ", opt$chr, " Position (Mb)"),
    y = expression(-log[10](italic(p)))
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# 5. 保存图片
ggsave(opt$output, plot = p, width = 8, height = 5, dpi = 300)
cat("局部放大图保存至:", opt$output, "\n")