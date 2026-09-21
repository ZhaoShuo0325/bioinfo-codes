#!/usr/bin/env Rscript
# Usage:
#    Rscript 06_4_run_qtl.R <control_yaml_file> <output_pdf>
#    - control_yaml_file: Path to the YAML control file.
#    - output_pdf: Path to the output PDF file.
#

args <- commandArgs(trailingOnly = TRUE)

control_file <- args[1]
output_pdf <- args[2]

library(qtl2)

cross <- read_cross2(control_file)
pr <- calc_genoprob(cross, map = cross$pmap)

out <- scan1(pr, cross$pheno)

print(find_peaks(out, map = cross$pmap, threshold = 3))

pdf(output_pdf, width = 11, height = 6)
for (i in colnames(out)) {
  plot(out, map = cross$pmap, lodcolumn = i, 
       main = paste("QTL Scan for Phenotype:", i),
       col = "#2b5c8f", lwd = 1.5,
       ylab = "LOD Score", xlab = "Chromosome")
  abline(h = 3, col = "red", lty = 2, lwd = 1.2)
}
dev.off()