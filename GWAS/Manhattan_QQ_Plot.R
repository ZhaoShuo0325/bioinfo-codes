#!/usr/bin/env Rscript
# ==============================================================================
# Script Name: qqman_plot.R
# Description: 用于绘制 GWAS 曼哈顿图（Manhattan Plot）与 Q-Q 图的 R 脚本。
# 
# Usage:
#   Rscript plot_gwas.R [options]
#
# Options:
#   -i, --input <file>    输入 GWAS 关联分析结果文件路径 
#                         [默认: ../03_GWAS/output4/gwas_all.assoc.txt]
#   -o, --output <file>   输出曼哈顿图的路径/名称 
#                         [默认: ../03_GWAS/output4/Manhattan_plot_beauty.png]
#   -a, --alpha <num>     显著性分子，可选 0.05 或 0.01 
#                         [默认: 0.05]
#   -t, --title <str>     自定义曼哈顿图的标题 
#                         [默认: 自动生成 "SV-based GWAS (N = ...)" ]
#   -h, --help            显示此帮助信息并退出
#
# Examples:
#   1. 使用默认参数直接运行：
#      Rscript qqman_plot.R
#
#   2. 自定义输入、输出、显著性线（0.01/N）及图表标题：
#      Rscript qqman_plot.R -i "../03_GWAS/output4/gwas_all.assoc.txt" \
#                          -o "../03_GWAS/output4/My_Manhattan.png" \
#                          -a 0.01 \
#                          -t "GWAS of Phenotype m1186"
# ==============================================================================
library(qqman)
library(optparse)

# 1. 定义命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "../03_GWAS/output4/gwas_all.assoc.txt",
              help = "输入 GWAS 关联分析结果文件路径 [默认: %default]", metavar = "file"),
  make_option(c("-o", "--output"), type = "character", default = "../03_GWAS/output4/Manhattan_plot_beauty.png",
              help = "输出曼哈顿图的路径/名称 [默认: %default]", metavar = "file"),
  make_option(c("-a", "--alpha"), type = "numeric", default = 0.05,
              help = "显著性分子，可选 0.05 或 0.01 [默认: %default]", metavar = "number"),
  make_option(c("-t", "--title"), type = "character", default = NULL,
              help = "曼哈顿图的标题", metavar = "string")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 2. 读取并清洗数据
cat("正在读取数据:", opt$input, "\n")
gwas_res <- read.table(opt$input, header = TRUE, stringsAsFactors = FALSE)
gwas_res <- gwas_res[!is.na(gwas_res$p_wald), ]

n_snps <- nrow(gwas_res)                      
sig_threshold <- opt$alpha / n_snps            # 动态计算阈值 (0.05/N 或 0.01/N)

cat("总位点数 (N):", n_snps, "\n")
cat(paste0("-log10(", opt$alpha, "/N) 阈值:"), -log10(sig_threshold), "\n")

# 3. 动态处理图表标题
if (is.null(opt$title)) {
  plot_title <- paste0("SV-based GWAS (N = ", format(n_snps, big.mark = ","), ")")
} else {
  plot_title <- opt$title
}

# 4. 绘制美化后的曼哈顿图
cat("正在生成曼哈顿图...\n")
png(opt$output, width = 2400, height = 1200, res = 300)

max_y <- ceiling(max(-log10(gwas_res$p_wald), na.rm = TRUE)) + 1
if (max_y < 8) max_y <- 8

manhattan(
  gwas_res,
  chr = "chr",
  bp = "ps",
  p = "p_wald",
  snp = "rs",
  col = c("#2A4494", "#4EB1D3"), 

  genomewideline = FALSE,
  suggestiveline = FALSE,

  ylim = c(0, max_y),
  cex = 0.6,          
  cex.axis = 0.9,     
  chrlabs = c(1:12, ""), 
  main = plot_title
)
abline(h = -log10(sig_threshold), col = "black", lty = 2, lwd = 1.2)
dev.off()

# 5. 绘制对应的 Q-Q 图（自动根据曼哈顿图输出路径生成同目录下的 QQ 图）
qq_output <- sub("Manhattan", "QQ", opt$output)
cat("正在生成 Q-Q 图...\n")
png(qq_output, width = 1200, height = 1200, res = 300)
qq(
  gwas_res$p_wald, 
  main = "Q-Q Plot of GWAS p-values",
  col = "#2A4494",       
  pch = 19,              
  cex = 0.6,
  cex.axis = 0.9
)
dev.off()

cat(" 曼哈顿图已保存至:", opt$output, "\n")
cat(" Q-Q 图已保存至:", qq_output, "\n")