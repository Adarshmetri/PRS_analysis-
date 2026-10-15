#!/usr/bin/env bash
set -euo pipefail

BASE_RAW="${1:-data/Height.gwas.txt.gz}"
OUT_DIR="${2:-results}"

mkdir -p "${OUT_DIR}"
echo "==> [Step 1] Base Data QC on ${BASE_RAW}"

# 1. Filter MAF > 0.01 and INFO > 0.8
echo "Filtering MAF > 0.01 and INFO > 0.8..."
gunzip -c "${BASE_RAW}" | \
awk 'NR==1 || ($11 > 0.01 && $10 > 0.8) {print}' | \
gzip > "${OUT_DIR}/Height.filtered.gz"

# 2. Remove duplicate SNPs (Column 3 = SNP ID)
echo "Removing duplicate SNPs..."
gunzip -c "${OUT_DIR}/Height.filtered.gz" | \
awk '{seen[$3]++; if(seen[$3]==1){ print}}' | \
gzip > "${OUT_DIR}/Height.nodup.gz"

# 3. Remove ambiguous SNPs (A/T, T/A, G/C, C/G)
echo "Removing strand-ambiguous SNPs..."
gunzip -c "${OUT_DIR}/Height.nodup.gz" | \
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' | \
gzip > "${OUT_DIR}/Height.QC.gz"

rm -f "${OUT_DIR}/Height.filtered.gz" "${OUT_DIR}/Height.nodup.gz"
echo "==> Base Data QC complete: ${OUT_DIR}/Height.QC.gz"
