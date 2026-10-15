#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
base_file <- ifelse(length(args) > 0, args[1], "results/Height.QC.gz")
out_file  <- ifelse(length(args) > 1, args[2], "results/Height.QC.Transformed")

dat <- fread(base_file)
if ("OR" %in% colnames(dat)) {
  dat[, BETA := log(OR)]
}
fwrite(dat, out_file, sep = "\t")
cat("Transformed effect sizes written to", out_file, "\n")
