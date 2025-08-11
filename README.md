# Pathway Overlap – Nextflow pipeline

A Nextflow pipeline to evaluate pathway/gene-set enrichment across GWAS traits using MAGMA and PRSet, generate empirical p-values via two null models (BiReWire and KeepPathSize), and assess biological relevance with OpenTargets evidence and tissue-specific expression.

Key capabilities:
- Harmonize GWAS inputs (with optional SNP de-duplication).
- Run pathway enrichment via MAGMA or PRSet.
- Build null distributions using two strategies and compute empirical p-values and standardized effect sizes.
- Compare top pathways against OpenTargets gene–disease evidence using size-matched controls; create summary plots.
- Evaluate tissue specificity using precomputed tissue expression specificity tables.
- Perform correlation analyses across ranking methods and tissues (optional).

---

## Requirements

- Nextflow (>= 22.10) and Java 11+
- R (>= 4.1) with packages:
  - tidyverse, data.table, plyr, jsonlite, GSA, MatchIt, ggplot2, patchwork, gridExtra
- MAGMA (if `--run_magma true`) and a reference PLINK bfile (e.g., 1000G EUR)
- PRSet (if `--run_prset true`) and compatible GWAS formats
- LSF scheduler (default config) or edit `nextflow.config` to use your scheduler/local
- Data dependencies:
  - GMT gene set file(s) for real and randomized gene sets
  - OpenTargets association JSONs (associationByDatatypeDirect) if running OpenTargets analyses
  - Tissue specificity data tables (GeneExpressionLandscape specificity CSVs)

---

## Repository layout

- `main.nf` – pipeline workflow
- `nextflow.config` – default parameters and process settings (LSF by default)
- `json_files/GWAS_input.json` – example GWAS input configuration
- `modules/` – Nextflow modules (magma, prset, empirical_pval, tissuespecificity, opentargets)
- `scripts/` – R scripts used by modules
- `magma_tissuespec/` – example outputs for tissue specificity by trait

---

## Inputs

### 1) GWAS list JSON
Provide one or more traits with paths and column mappings. See `json_files/GWAS_input.json` for a full example. Minimal schema per trait:

```json
{
  "trait": "t2d",
  "gwas_file": "/path/to/gwas.txt",
  "rsid_col": "ID",
  "chr_col": "CHR",
  "pos_col": "POS",
  "pval_col": "PVAL",
  "n_col": "N",
  "effect_allele": "A1",
  "other_allele": "A2",
  "binary_target": "F",
  "summary_statistic_name": "BETA",
  "summary_statistic_type": "beta"
}
```

Ensure column names match your files. Compressed inputs supported if the corresponding tools accept them.

### 2) Gene sets
- Real gene sets (GMT): `params.geneset_real`
- Randomized sets directories: `params.gmt_dirs.birewire` and `params.gmt_dirs.keeppathsize`

### 3) MAGMA reference
- `params.bfile`: prefix to PLINK files used by MAGMA

### 4) OpenTargets data
- The OpenTargets statistics script expects association JSON files (e.g., associationByDatatypeDirect). Update paths in scripts or pass through the pipeline module if customized.

### 5) Tissue specificity data
- `params.tissue_expression_data` should point to the directory/file(s) with gene-level specificity metrics used by the tissue specificity module.

---

## Configuration

Key parameters in `nextflow.config` (edit to your environment):

- General
  - `params.traits_config` – path to GWAS JSON list
  - `params.outdir` – output root directory
  - `params.num_random_sets` – number of permutations per trait/method
  - `params.enrichment_method` – "magma" or "prset"
  - `params.randomization_methods` – `[ 'birewire', 'keeppathsize' ]`
  - `params.gene_files` – background and real gene region files (used by tools)
  - `params.gmt_dirs` – directories for randomized GMTs
  - `params.geneset_real` – real GMT file
  - `params.bfile` – MAGMA reference PLINK prefix
  - `params.ukb_dir` – directory for UKB-derived helper data if used

- Workflow toggles
  - `params.run_magma` – run MAGMA
  - `params.run_prset` – run PRSet
  - `params.run_empirical` – compute empirical p-values and standardized effect sizes
  - `params.run_tissue_specificity` – run tissue specificity analysis
  - `params.run_opentargets` – run OpenTargets size-matched comparisons
  - `params.run_ot_viz` – generate OpenTargets plots (default false)
  - `params.run_ot_correlation` – run OpenTargets correlation analyses across ranking methods
  - `params.run_tissue_correlation` – run tissue specificity correlation analyses

- OpenTargets
  - `params.opentargets_supported_traits` – whitelist of traits
  - `params.opentargets_n_values` – top-N pathway sizes to test (comma-separated string, e.g., "10,20,50,100")
  - `params.opentargets_sig_threshold` – significance threshold used in plots (default 0.00625)

- Schedulers
  - Default is LSF with queues and project codes per process. If you are not on LSF, remove the per-process `executor/queue/clusterOptions` blocks and configure for your scheduler (e.g., `executor = 'local'`, `slurm`, etc.).

---

## Usage

Basic run (uses values from `nextflow.config`):

```
nextflow run main.nf -resume
```

Override common parameters at runtime:

```
nextflow run main.nf \
  --traits_config ./json_files/GWAS_input.json \
  --outdir /path/to/results \
  --enrichment_method magma \
  --num_random_sets 1000 \
  --run_magma true --run_prset true --run_empirical true \
  --run_opentargets true --run_ot_viz true \
  --run_tissue_specificity true --run_ot_correlation true --run_tissue_correlation true \
  -resume
```

Run a subset of stages, e.g., MAGMA + empirical only:

```
nextflow run main.nf --run_magma true --run_prset false --run_empirical true --run_opentargets false --run_tissue_specificity false -resume
```

Note: If not using LSF, edit `nextflow.config` to your scheduler before running.

---

## What the pipeline does

1) GWAS preprocessing
- Optional SNP de-duplication per trait

2) Enrichment analyses
- MAGMA gene-set analysis (requires `params.bfile` and real/random GMTs)
- PRSet analysis (if enabled)

3) Null models and empirical p-values
- Randomized gene sets via BiReWire and KeepPathSize
- Empirical p-values calculated by `scripts/calc_empirical.r`
  - Outputs per trait/method: `<trait>_<tool>_empirical_pvalues.txt`
  - Includes standardized effect sizes and pathway sizes when available

4) OpenTargets size-matched comparison (optional)
- Uses empirical results from both null models
- Computes metrics and produces summaries and plots
  - Detailed metrics: `<trait>_<tool>_*_detailed_advantage.csv`
  - Combined plots: `<trait>_combined_advantage_plots.pdf`, `<trait>_combined_boxplots.pdf`

5) Tissue specificity analysis (optional)
- Computes pathway-level tissue expression summaries
- Summary TSV: `<trait>_<tool>_size_matched_analysis_summary.tsv`
- Per-trait outputs organized under `<outdir>`; examples in `magma_tissuespec/`

6) Correlation analyses (optional)
- OpenTargets correlation across ranking methods (enable with `--run_ot_correlation true`)
  - Correlates pathway rankings/metrics across methods to assess concordance
  - Produces summary tables and plots (see `modules/opentargets` outputs)
- Tissue specificity correlation (enable with `--run_tissue_correlation true`)
  - Runs `scripts/tissue_correlation_stats.R`
  - Computes Spearman correlations between pathway ranks and tissue specificity metrics
  - Outputs per-trait CSV summaries and multi-page PDFs per method/subset

---

## Outputs (typical)

- Empirical results: `results/<trait>/<tool>/<trait>_<tool>_empirical_pvalues.txt`
- OpenTargets:
  - `<trait>_<tool>_gene_disease_associations.csv`
  - `<trait>_<tool>_*_detailed_advantage.csv`
  - `<trait>_combined_advantage_plots.pdf`, `<trait>_combined_boxplots.pdf`
  - OpenTargets correlation summaries/plots (if enabled)
- Tissue specificity:
  - `<trait>_<tool>_all_detailed_metrics.csv`
  - `<trait>_<tool>_size_matched_analysis_summary.tsv`
  - Tissue correlation summaries and plots (if enabled):
    - `<trait>_<tool>_tissue_correlation_summary.csv`
    - `<trait>_<tool>_best_method_by_tissue.csv`
    - `<trait>_<tool>_<method>_<subset>_tissue_correlation_*.pdf`

Note: Exact paths depend on how modules in `main.nf` stage outputs into `params.outdir`.

---

## Customization tips

- Update absolute paths in `nextflow.config` to your environment (bfile, gene sets, OpenTargets data, tissue data).
- If your OpenTargets JSONs live elsewhere, modify the path used in `scripts/size_matched_OT_stats_optimized.R` or expose it as a Nextflow parameter.
- For tissue correlation, ensure `params.tissue_expression_data` points to the correct CSV with a `Name` column for gene symbols and tissue columns.
- For non-LSF systems, simplify the `process {}` blocks in `nextflow.config` and set a global `executor`.
- Adjust `params.opentargets_supported_traits` to restrict which traits are processed in the OpenTargets stage.

---

## Troubleshooting

- No random result files found: ensure randomized GMT directories are correct and that randomization steps have run (BiReWire/KeepPathSize).
- Missing columns errors in R scripts: confirm your GWAS JSON column mappings and that tool-specific outputs contain the expected fields.
- MAGMA failures: verify `params.bfile` and that MAGMA binaries are available on PATH or referenced in the module.
- Correlation step warnings about missing metrics: check tissue file columns and that empirical results include required fields (e.g., `empirical_pval`, `std_effect_size`).
- Scheduler issues: adapt `nextflow.config` executors and resource directives to your environment.

---

## Citation

If you use this pipeline, please cite MAGMA, PRSet, and OpenTargets appropriately, and reference this repository.

## License

See `LICENSE`.