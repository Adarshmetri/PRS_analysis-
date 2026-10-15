#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${1:-results}"
echo "==> [Step 3] Clumping and PRS calculation with PLINK"

# 1. Log transform OR to BETA
Rscript scripts/transform_effects.R "${OUT_DIR}/Height.QC.gz" "${OUT_DIR}/Height.QC.Transformed"

# 2. Clumping (r2 < 0.1, 250kb window)
plink \
    --bfile "${OUT_DIR}/EUR.QC" \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump "${OUT_DIR}/Height.QC.Transformed" \
    --clump-snp-field SNP \
    --clump-field P \
    --out "${OUT_DIR}/EUR"

awk 'NR!=1{print $3}' "${OUT_DIR}/EUR.clumped" > "${OUT_DIR}/EUR.valid.snp"

# 3. Extract P-values & threshold ranges
awk '{print $3,$8}' "${OUT_DIR}/Height.QC.Transformed" > "${OUT_DIR}/SNP.pvalue"

cat << 'EOF' > "${OUT_DIR}/range_list"
0.001 0 0.001
0.05 0 0.05
0.1 0 0.1
0.2 0 0.2
0.3 0 0.3
0.4 0 0.4
0.5 0 0.5
EOF

# 4. Generate PRS across thresholds
plink \
    --bfile "${OUT_DIR}/EUR.QC" \
    --score "${OUT_DIR}/Height.QC.Transformed" 3 4 12 header \
    --q-score-range "${OUT_DIR}/range_list" "${OUT_DIR}/SNP.pvalue" \
    --extract "${OUT_DIR}/EUR.valid.snp" \
    --out "${OUT_DIR}/EUR"

echo "==> PRS profile scores generated in ${OUT_DIR}/EUR.*.profile"
