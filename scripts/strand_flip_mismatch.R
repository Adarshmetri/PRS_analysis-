#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
})

args <- commandArgs(trailingOnly = TRUE)
bim_file   <- ifelse(length(args) > 0, args[1], "data/EUR.bim")
base_qc    <- ifelse(length(args) > 1, args[2], "results/Height.QC.gz")
snplist    <- ifelse(length(args) > 2, args[3], "results/EUR.QC.snplist")
out_prefix <- ifelse(length(args) > 3, args[4], "results/EUR")

bim <- fread(bim_file) %>%
  setnames(., colnames(.), c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")) %>%
  .[, c("B.A1", "B.A2") := list(toupper(B.A1), toupper(B.A2))]

height <- fread(base_qc) %>%
  .[, c("A1", "A2") := list(toupper(A1), toupper(A2))]

qc <- fread(snplist, header = FALSE)

info <- merge(bim, height, by = c("SNP", "CHR", "BP")) %>%
  .[SNP %in% qc[[1]]]

complement <- function(x) {
  c("A" = "T", "T" = "A", "C" = "G", "G" = "C")[x]
}

info.match <- info[A1 == B.A1 & A2 == B.A2, SNP]
com.snps   <- info[complement(B.A1) == A1 & complement(B.A2) == A2, SNP]

# Strand-flip
bim[SNP %in% com.snps, c("B.A1", "B.A2") := list(complement(B.A1), complement(B.A2))]

# Recode
recode.snps <- info[B.A1 == A2 & B.A2 == A1, SNP]
bim[SNP %in% recode.snps, c("B.A1", "B.A2") := list(B.A2, B.A1)]

# Recode + Strand-flip
com.recode <- info[complement(B.A1) == A2 & complement(B.A2) == A1, SNP]
bim[SNP %in% com.recode, c("B.A1", "B.A2") := list(complement(B.A2), complement(B.A1))]

fwrite(bim[, c("SNP", "B.A1")], paste0(out_prefix, ".a1"), col.names = FALSE, sep = "\t")

mismatch <- bim[!(SNP %in% info.match |
                  SNP %in% com.snps |
                  SNP %in% recode.snps |
                  SNP %in% com.recode), SNP]

fwrite(data.table(SNP = mismatch), paste0(out_prefix, ".mismatch"), col.names = FALSE)
cat("Identified", length(mismatch), "mismatching SNPs.\n")
