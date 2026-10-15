#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
})

args <- commandArgs(trailingOnly = TRUE)
target_raw <- ifelse(length(args) > 0, args[1], "data/EUR")
res_dir    <- ifelse(length(args) > 1, args[2], "results")

# Calculate 6 PCs for population stratification if not already present
eigenvec <- file.path(res_dir, "EUR.eigenvec")
if (!file.exists(eigenvec)) {
  system(paste0("plink --bfile ", res_dir, "/EUR.QC --indep-pairwise 200 50 0.25 --out ", res_dir, "/EUR"))
  system(paste0("plink --bfile ", res_dir, "/EUR.QC --extract ", res_dir, "/EUR.prune.in --pca 6 --out ", res_dir, "/EUR"))
}

p.threshold <- c(0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)

phenotype <- fread(paste0(target_raw, ".height"))
covariate <- fread(paste0(target_raw, ".cov"))
pcs <- fread(eigenvec, header = FALSE) %>%
  setnames(., colnames(.), c("FID", "IID", paste0("PC", 1:6)))

pheno <- merge(phenotype, covariate, by = c("FID", "IID")) %>%
  merge(., pcs, by = c("FID", "IID"))

# Fit Null model (without PRS)
null.model <- lm(Height ~ ., data = pheno[, !c("FID", "IID"), with = FALSE])
null.r2 <- summary(null.model)$r.squared

prs.result <- data.table()

for (i in p.threshold) {
  prof_file <- file.path(res_dir, paste0("EUR.", i, ".profile"))
  if (!file.exists(prof_file)) next
  
  prs <- fread(prof_file)[, c("FID", "IID", "SCORE")]
  pheno.prs <- merge(pheno, prs, by = c("FID", "IID"))
  
  model <- lm(Height ~ ., data = pheno.prs[, !c("FID", "IID"), with = FALSE])
  model.sum <- summary(model)
  
  model.r2 <- model.sum$r.squared
  prs.r2   <- model.r2 - null.r2
  prs.coef <- model.sum$coefficients["SCORE", ]
  
  prs.result <- rbind(prs.result, data.table(
    Threshold = i,
    R2        = prs.r2,
    P         = as.numeric(prs.coef[4]),
    BETA      = as.numeric(prs.coef[1]),
    SE        = as.numeric(prs.coef[2])
  ))
}

best_row <- prs.result[which.max(R2)]
cat("\n========================================\n")
cat("Best-fit P-value Threshold:", best_row$Threshold, "\n")
cat("Variance Explained (PRS R2):", round(best_row$R2, 5), "\n")
cat("P-value:", best_row$P, "\n")
cat("========================================\n\n")

fwrite(prs.result, file.path(res_dir, "prs_evaluation_results.csv"))
saveRDS(prs.result, file.path(res_dir, "prs_result.rds"))
