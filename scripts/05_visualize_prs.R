#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
res_dir    <- ifelse(length(args) > 0, args[1], "results")
target_raw <- ifelse(length(args) > 1, args[2], "data/EUR")

prs.result <- readRDS(file.path(res_dir, "prs_result.rds"))

prs.result[, print.p := round(P, digits = 3)]
prs.result[print.p == 0, print.p := format(P, digits = 2)]
prs.result[, print.p := sub("e", "*x*10^", print.p)]

# 1. Bar Plot of Model Fit (R2)
p_bar <- ggplot(data = prs.result, aes(x = factor(Threshold), y = R2)) +
  geom_bar(aes(fill = -log10(P)), stat = "identity", width = 0.6) +
  geom_text(aes(label = print.p), vjust = -1.2, hjust = 0.5, size = 3.8, parse = TRUE) +
  scale_y_continuous(limits = c(0, max(prs.result$R2) * 1.25)) +
  scale_fill_gradient2(
    low = "dodgerblue", mid = "darkorange", high = "firebrick",
    midpoint = median(-log10(prs.result$P)),
    name = expression(-log[10](italic(P)))
  ) +
  labs(
    x = expression(italic(P) - value ~ threshold ~ (italic(P)[T])),
    y = expression(paste("PRS model fit:  ", Delta, R ^ 2)),
    title = "Polygenic Risk Score Model Fit Across P-Value Thresholds"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

ggsave(file.path(res_dir, "EUR.height.bar.png"), plot = p_bar, width = 8, height = 6, dpi = 300)

# 2. Scatter Plot for Best Threshold
best_th <- prs.result[which.max(R2), Threshold]
best_profile <- fread(file.path(res_dir, paste0("EUR.", best_th, ".profile")))
pheno <- fread(paste0(target_raw, ".height"))
merged <- merge(pheno, best_profile, by = c("FID", "IID"))

p_scatter <- ggplot(merged, aes(x = SCORE, y = Height)) +
  geom_point(alpha = 0.6, color = "dodgerblue4") +
  geom_smooth(method = "lm", color = "firebrick", fill = "pink", se = TRUE) +
  labs(
    x = paste0("Best PRS (Threshold = ", best_th, ")"),
    y = "Observed Height",
    title = "Observed Phenotype vs Best Polygenic Risk Score"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(res_dir, "EUR.height.scatter.png"), plot = p_scatter, width = 7, height = 6, dpi = 300)
cat("Visualizations saved to", res_dir, "\n")
