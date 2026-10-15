#!/usr/bin/env bash
set -euo pipefail

echo "=========================================================="
echo " Starting Polygenic Risk Score (PRS) Pipeline (Choi et al.)"
echo "=========================================================="

DATA_DIR="data"
RES_DIR="results"

mkdir -p "${DATA_DIR}" "${RES_DIR}"

# Step 1: Base Data QC
bash scripts/01_qc_base.sh "${DATA_DIR}/Height.gwas.txt.gz" "${RES_DIR}"

# Step 2: Target Data QC
bash scripts/02_qc_target.sh "${DATA_DIR}/EUR" "${RES_DIR}"

# Step 3: PRS Clumping & Profile Generation
bash scripts/03_run_prs_plink.sh "${RES_DIR}"

# Step 4: Regression & Best Threshold Selection
Rscript scripts/04_evaluate_prs.R "${DATA_DIR}/EUR" "${RES_DIR}"

# Step 5: Visualization
Rscript scripts/05_visualize_prs.R "${RES_DIR}" "${DATA_DIR}/EUR"

echo "=========================================================="
echo " Pipeline Finished Successfully! Results saved in ${RES_DIR}/"
echo "=========================================================="
