#!/usr/bin/env bash
set -euo pipefail

TARGET_PREFIX="${1:-data/EUR}"
OUT_DIR="${2:-results}"

mkdir -p "${OUT_DIR}"
echo "==> [Step 2] Target Data QC on ${TARGET_PREFIX}"

# 1. Variant and Sample Missingness, MAF, HWE
plink \
    --bfile "${TARGET_PREFIX}" \
    --maf 0.01 \
    --hwe 1e-6 \
    --geno 0.01 \
    --mind 0.01 \
    --write-snplist \
    --make-just-fam \
    --out "${OUT_DIR}/EUR.QC"

# 2. LD Pruning for sample-level QC
plink \
    --bfile "${TARGET_PREFIX}" \
    --keep "${OUT_DIR}/EUR.QC.fam" \
    --extract "${OUT_DIR}/EUR.QC.snplist" \
    --indep-pairwise 200 50 0.25 \
    --out "${OUT_DIR}/EUR.QC"

# 3. Heterozygosity rate check
plink \
    --bfile "${TARGET_PREFIX}" \
    --extract "${OUT_DIR}/EUR.QC.prune.in" \
    --keep "${OUT_DIR}/EUR.QC.fam" \
    --het \
    --out "${OUT_DIR}/EUR.QC"

# Filter heterozygosity outliers (+/- 3 SD)
Rscript -e "
library(data.table)
dat <- fread('${OUT_DIR}/EUR.QC.het')
valid <- dat[F <= mean(F) + 3*sd(F) & F >= mean(F) - 3*sd(F)]
fwrite(valid[, c('FID', 'IID')], '${OUT_DIR}/EUR.valid.sample', sep='\t')
"

# 4. Resolve strand flipping and mismatches with Base QC
Rscript scripts/strand_flip_mismatch.R "${TARGET_PREFIX}.bim" "${OUT_DIR}/Height.QC.gz" "${OUT_DIR}/EUR.QC.snplist" "${OUT_DIR}/EUR"

# 5. Sex check
plink \
    --bfile "${TARGET_PREFIX}" \
    --extract "${OUT_DIR}/EUR.QC.prune.in" \
    --keep "${OUT_DIR}/EUR.valid.sample" \
    --check-sex \
    --out "${OUT_DIR}/EUR.QC"

Rscript -e "
library(data.table)
valid <- fread('${OUT_DIR}/EUR.valid.sample')
dat <- fread('${OUT_DIR}/EUR.QC.sexcheck')[FID %in% valid\$FID]
fwrite(dat[STATUS == 'OK', c('FID', 'IID')], '${OUT_DIR}/EUR.QC.valid', sep='\t')
"

# 6. Relatedness filtering (cutoff 0.125)
plink \
    --bfile "${TARGET_PREFIX}" \
    --extract "${OUT_DIR}/EUR.QC.prune.in" \
    --keep "${OUT_DIR}/EUR.QC.valid" \
    --rel-cutoff 0.125 \
    --out "${OUT_DIR}/EUR.QC"

# 7. Generate final QC target binary fileset
plink \
    --bfile "${TARGET_PREFIX}" \
    --make-bed \
    --keep "${OUT_DIR}/EUR.QC.rel.id" \
    --extract "${OUT_DIR}/EUR.QC.snplist" \
    --exclude "${OUT_DIR}/EUR.mismatch" \
    --a1-allele "${OUT_DIR}/EUR.a1" \
    --out "${OUT_DIR}/EUR.QC"

echo "==> Target Data QC complete: ${OUT_DIR}/EUR.QC.{bed,bim,fam}"
