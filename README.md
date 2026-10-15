# Polygenic Risk Score (PRS) Pipeline

An end-to-end, reproducible pipeline for calculating Polygenic Risk Scores (PRS) using PLINK (Clumping + Thresholding), based on [Choi et al., *Nature Protocols* (2020)](https://www.nature.com/articles/s41596-020-0353-1) and the [Choi et al. PRS Tutorial](https://choishingwan.github.io/PRS-Tutorial/).

## Features
- **Base Data QC**: Filters MAF (> 0.01), INFO score (> 0.8), duplicate SNPs, and strand-ambiguous SNPs.
- **Target Data QC**: Missingness filtering, HWE, extreme heterozygosity check (±3 SD), strand-flipping, sex discordance check, and relatedness filtering ($\hat{\pi} > 0.125$).
- **PRS Scoring (C+T)**: Accurate effect size transformation ($\beta = \ln(\text{OR})$), LD clumping ($r^2 < 0.1$, 250kb), and scoring across multiple $P$-value thresholds.
- **Phenotypic Association**: Population stratification correction using 6 genotype PCs, incremental variance explained ($\Delta R^2$) estimation, and identification of the optimal threshold.
- **Publication Figures**: Generates high-resolution bar charts and scatter plots.

---

## Directory Structure

```text
prs-pipeline/
├── .github/
│   └── workflows/
│       └── ci.yml               # GitHub Actions CI workflow
├── data/                        # Raw base and target datasets (ignored by git)
├── scripts/
│   ├── 01_qc_base.sh            # Base GWAS Summary Stats QC
│   ├── 02_qc_target.sh          # Target genotype/phenotype QC (PLINK)
│   ├── strand_flip_mismatch.R   # Identifies strand flips and mismatches
│   ├── 03_run_prs_plink.sh      # LD clumping & PRS calculations
│   ├── transform_effects.R      # Log-transforms OR to Beta
│   ├── 04_evaluate_prs.R        # Population PCs, regression, best threshold
│   └── 05_visualize_prs.R       # Generates ggplot2 figures
├── results/                     # Intermediate QC files, profiles, and plots
├── .gitignore                   # Keeps large genomics data off GitHub
├── environment.yml              # Conda environment recipe
├── run_pipeline.sh              # Master pipeline runner (Linux / macOS / Git Bash)
├── run_pipeline.bat             # Windows CMD launcher
└── README.md
```

---

## Quickstart

### 1. Set Up Dependencies

Using Conda / Mamba:
```bash
conda env create -f environment.yml
conda activate prs-pipeline
```

Or install tools individually:
- **PLINK 1.9**: [https://www.cog-genomics.org/plink/1.9/](https://www.cog-genomics.org/plink/1.9/)
- **R packages**: `data.table`, `magrittr`, `ggplot2`


### 2. Run Pipeline

- **On Linux / macOS / Git Bash**:
  ```bash
  chmod +x run_pipeline.sh scripts/*.sh scripts/*.R
  ./run_pipeline.sh
  ```

- **On Windows Command Prompt**:
  ```cmd
  run_pipeline.bat
  ```

---

## Outputs

All output files are saved in `results/`:
- `Height.QC.gz`: Quality-controlled base summary statistics.
- `EUR.QC.{bed,bim,fam}`: Cleaned target genotype files.
- `EUR.<threshold>.profile`: Polygenic scores for each target sample across thresholds.
- `prs_evaluation_results.csv`: Table of incremental $R^2$, regression $\beta$, standard error, and $P$-values across all thresholds.
- `EUR.height.bar.png`: Publication-ready bar chart of PRS model fit ($\Delta R^2$).
- `EUR.height.scatter.png`: Scatter plot with regression line of best PRS score vs observed height.

---

## Citation

If you use this pipeline or tutorial data in your research, please cite:

> Choi, S.W., Mak, T.S. & O’Reilly, P.F. Tutorial: a guide to performing polygenic risk score analyses. *Nat Protoc* **15**, 2759–2772 (2020). [https://doi.org/10.1038/s41596-020-0353-1](https://doi.org/10.1038/s41596-020-0353-1)
