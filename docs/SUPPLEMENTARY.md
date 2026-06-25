# Supplementary Analyses & Figure Generation

This document maps each **supplementary analysis / figure / table** to the script
that produces it, the command to run it, and its inputs/outputs.

The core pipeline (`nextflow run main.nf`) automates real + null enrichment,
empirical statistics, and the validation **correlation stats** (OpenTargets, GTEx
tissue, MalaCards, DoRothEA). The items below are **not fully push-button**: some
run as Nextflow modules, others are standalone scripts a user runs by hand after
the core pipeline finishes.

> **How to use this doc:** fill in the `Manuscript item` column with the actual
> figure/table number from the paper (marked `TODO`), then verify each command
> once on your own results so a labmate can copy-paste it.

## Legend

- 🟢 **Automated** — runs as part of `nextflow run main.nf` (listed here for traceability).
- 🟡 **Standalone** — run manually with the command shown, after the core pipeline.
- ✏️ **Has hardcoded paths** — edit the marked lines before running on a new system.

---

## Quick reference table

| Manuscript item | Status | Script | Produces |
|---|---|---|---|
| TODO (fig/table #) | 🟡 ✏️ | [scripts/pathsize_rank_correlation.R](../scripts/pathsize_rank_correlation.R) | Pathway-size vs rank: heatmap, scatter, dotplot + summary table |
| TODO | 🟡 | [scripts/validation/opentargets/opentargets.R](../scripts/validation/opentargets/opentargets.R) | OpenTargets comparison figure (PDF) |
| TODO | 🟢 | [scripts/validation/opentargets/OT_correlation_stats.R](../scripts/validation/opentargets/OT_correlation_stats.R) | OpenTargets rank-correlation summary CSVs |
| TODO | 🟢 | [scripts/validation/tissue/tissue_correlation_stats.R](../scripts/validation/tissue/tissue_correlation_stats.R) | GTEx tissue rank-correlation summary CSVs |
| TODO | 🟢 | [scripts/validation/malacards/malacards_correlation.R](../scripts/validation/malacards/malacards_correlation.R) | MalaCards rank-correlation summary CSV |
| TODO | 🟢 | [scripts/validation/dorothea/dorothea_correlation.R](../scripts/validation/dorothea/dorothea_correlation.R) | DoRothEA rank-correlation summary CSV |
| TODO | 🟢 | [scripts/core/calc_fpr.r](../scripts/core/calc_fpr.r) | False-positive-rate table per trait/tool/method |
| TODO | 🟡 ✏️ | [scripts/validation/malacards/malacards_convert_gene_ID.R](../scripts/validation/malacards/malacards_convert_gene_ID.R) | Prep: MalaCards symbol→Ensembl mapping |
| TODO | 🟡 ✏️ | [scripts/validation/dorothea/dorothea.r](../scripts/validation/dorothea/dorothea.r) | Prep: DoRothEA pathway scores CSV |
| TODO | 🟡 | [scripts/multipath_sldsc.sh](../scripts/multipath_sldsc.sh), [scripts/multipath_sldsc2.sh](../scripts/multipath_sldsc2.sh), [scripts/setup_sldsc.sh](../scripts/setup_sldsc.sh), [scripts/run_partitioned_h2.sh](../scripts/run_partitioned_h2.sh) | Partitioned heritability / stratified LDSC |

---

## Details

### Pathway-size vs rank correlation (figures + table)
- **Script:** [scripts/pathsize_rank_correlation.R](../scripts/pathsize_rank_correlation.R) — 🟡 standalone, ✏️ hardcoded paths
- **Run:**
  ```bash
  Rscript scripts/pathsize_rank_correlation.R
  ```
- **✏️ Edit before running:** `base_path` (line ~6) and the per-tool excluded traits (line ~10).
- **Inputs:** `${base_path}/combined_empirical_pvalues_{magma,prset,gsamixer,pascal}.xlsx`
- **Outputs (written to working dir):**
  - `pathsize_rank_correlation_results_all_tools.txt`
  - `pathsize_correlation_heatmap_all_tools.png`
  - `pathsize_correlation_scatter_gsr_vs_original.png`
  - `pathsize_correlation_dotplot_by_tool.png`

### OpenTargets comparison figure
- **Script:** [scripts/validation/opentargets/opentargets.R](../scripts/validation/opentargets/opentargets.R) — 🟡 standalone
- **Run:**
  ```bash
  Rscript scripts/validation/opentargets/opentargets.R <trait> <tool_base> <birewire_results> <keeppathsize_results>
  ```
- **Inputs:** per-trait birewire + keeppathsize empirical p-value files (from core pipeline).
- **Outputs:** results table (`.tsv`) + comparison figure (`.pdf`) in working dir.

### OpenTargets correlation stats
- **Script:** [scripts/validation/opentargets/OT_correlation_stats.R](../scripts/validation/opentargets/OT_correlation_stats.R) — 🟢 via [modules/opentargets/opentargets_stats_correlation.nf](../modules/opentargets/opentargets_stats_correlation.nf)
- **Run manually (if needed):**
  ```bash
  Rscript scripts/validation/opentargets/OT_correlation_stats.R <trait> <tool_base> <birewire_results> <keeppathsize_results> <gmt_file> <opentargets_json_dir>
  ```
- **Outputs:** `{trait}_{tool}_gene_disease_associations.csv`, `{trait}_{tool}_rank_correlation_summary.csv`

### GTEx tissue correlation stats
- **Script:** [scripts/validation/tissue/tissue_correlation_stats.R](../scripts/validation/tissue/tissue_correlation_stats.R) — 🟢 via [modules/tissuespecificity/tissue_correlation.nf](../modules/tissuespecificity/tissue_correlation.nf)
- **Run manually (if needed):**
  ```bash
  Rscript scripts/validation/tissue/tissue_correlation_stats.R <trait> <tool_base> <birewire_results> <keeppathsize_results> <gmt_file> <tissue_file>
  ```
- **Outputs:** `{trait}_{tool}_tissue_correlation_summary.csv`, `..._best_method_by_tissue.csv`, `..._best_tissue_overall_by_method_subset.csv`

### MalaCards correlation stats
- **Script:** [scripts/validation/malacards/malacards_correlation.R](../scripts/validation/malacards/malacards_correlation.R) — 🟢 via [modules/malacards/malacards_correlation.nf](../modules/malacards/malacards_correlation.nf)
- **Run manually (if needed):**
  ```bash
  Rscript scripts/validation/malacards/malacards_correlation.R <trait> <tool_base> <malacards_path> <birewire_results> <keeppathsize_results> <gmt_file>
  ```
- **Prerequisite:** run the gene-ID conversion prep script first (below).
- **Outputs:** `{trait}_{tool}_malacards_rank_correlation_summary.csv`

### DoRothEA correlation stats
- **Script:** [scripts/validation/dorothea/dorothea_correlation.R](../scripts/validation/dorothea/dorothea_correlation.R) — 🟢 via [modules/dorothea/dorothea_correlation.nf](../modules/dorothea/dorothea_correlation.nf)
- **Run manually (if needed):**
  ```bash
  Rscript scripts/validation/dorothea/dorothea_correlation.R <trait> <tool_base> <dorothea_path> <birewire_results> <keeppathsize_results>
  ```
- **Prerequisite:** run the DoRothEA pathway-score prep script first (below).
- **Outputs:** `{trait}_{tool}_dorothea_rank_correlation_summary.csv`

### False positive rate (FPR)
- **Script:** [scripts/core/calc_fpr.r](../scripts/core/calc_fpr.r) — 🟢 via [modules/fpr/fpr_calculation.nf](../modules/fpr/fpr_calculation.nf)
- **Run manually (if needed):**
  ```bash
  Rscript scripts/core/calc_fpr.r <trait> <tool_base> <rand_method> <random_dir>
  ```
- **Outputs:** FPR results CSV.

### Prep — MalaCards gene-ID conversion
- **Script:** [scripts/validation/malacards/malacards_convert_gene_ID.R](../scripts/validation/malacards/malacards_convert_gene_ID.R) — 🟡 standalone, ✏️ hardcoded paths
- **✏️ Edit before running:** `input_path` and `out_dir` (lines ~6–7).
- **Run:**
  ```bash
  Rscript scripts/validation/malacards/malacards_convert_gene_ID.R
  ```
- **Outputs:** `malacards_mapping_summary.csv`, `malacards_symbol_to_ensembl_master.csv` (in `out_dir`).

### Prep — DoRothEA pathway scores
- **Script:** [scripts/validation/dorothea/dorothea.r](../scripts/validation/dorothea/dorothea.r) — 🟡 standalone, ✏️ hardcoded paths
- **✏️ Edit before running:** the hardcoded GMT path (line ~46).
- **Run:**
  ```bash
  Rscript scripts/validation/dorothea/dorothea.r
  ```
- **Outputs:** DoRothEA pathway scores CSV (consumed by `dorothea_correlation.R`).

### Partitioned heritability / stratified LDSC
- **Scripts:** [scripts/setup_sldsc.sh](../scripts/setup_sldsc.sh), [scripts/multipath_sldsc.sh](../scripts/multipath_sldsc.sh), [scripts/multipath_sldsc2.sh](../scripts/multipath_sldsc2.sh), [scripts/run_partitioned_h2.sh](../scripts/run_partitioned_h2.sh) — 🟡 standalone
- **Run:** TODO — confirm execution order and arguments, then document here.
- **Outputs:** TODO.

---

## Known manual steps / caveats

- Several standalone scripts contain **hardcoded absolute paths** (marked ✏️ above)
  pointing at the original compute environment. A new user must edit these.
- Standalone scripts write outputs to the **current working directory** unless an
  output path is hardcoded — `cd` into your desired results folder before running.
- The pascalx tool scripts (`scripts/tool_specific/pascalx/`) are wired into
  [modules/pascalx/pascalx.nf](../modules/pascalx/pascalx.nf); see [docs/TOOL_DETAILS.md](TOOL_DETAILS.md).
