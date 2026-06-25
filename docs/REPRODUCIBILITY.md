# Reproducibility Checklist

Step-by-step for a new scientist reproducing the two core analyses:
**(1) the GSR adjustment** and **(2) external-benchmark validation**. Work top to
bottom; nothing downstream runs until the inputs above it are in place.

---

## A. Software

- [ ] **Nextflow** ≥ 22.10 (DSL2) and **Java** 11+
- [ ] **R** ≥ 4.1 with packages:
  `tidyverse, data.table, plyr, jsonlite, GSA, MatchIt, ggplot2, patchwork, gridExtra, RColorBrewer, pheatmap`
- [ ] **Singularity/Apptainer** ≥ 3.0 (for GSA-MiXeR and PascalX containers)
- [ ] Tool binaries / containers: **MAGMA v1.10**, **PRSice/PRSet v2.3+**,
  **GSA-MiXeR** (`gsa-mixer.sif`), **PascalX** (`pascalx_singularity.sif`)

> Record exact versions of every tool and validation database — they affect results.

## B. Configure paths (single source of truth = `nextflow.config`)

Every absolute path in the repo is a `/path/to/...` placeholder. Edit
[`../nextflow.config`](../nextflow.config) and
[`../json_files/GWAS_input.json`](../json_files/GWAS_input.json) only — module and
script bodies no longer contain hardcoded data paths.

- [ ] `traits_config`, `outdir`, `num_random_sets` (study = 1000)
- [ ] References: `bfile`, `gtf_reference`, `gene_files`, `geneset_real`
- [ ] PRSet: `prsice_bin`, `prset_background`, `ukb_dir`
- [ ] GSA-MiXeR: `mixer_sif`, `mixer_ref_*`, `mixer_go_*`
- [ ] PascalX: `pascalx_sif`, `pascalx_ref_panel`, `pascalx_genome_annot`,
  and the `singularity.runOptions` bind mounts at the bottom of the config
- [ ] Validation: `opentargets_json_dir`, `tissue_expression_data`,
  `malacards_path`, `dorothea_scores_path`
- [ ] In `json_files/GWAS_input.json`: set each `gwas_file` and verify every `*_col`
  matches your summary-stat headers (case-sensitive)

See [`../data/README.md`](../data/README.md) for what each file is and where to get it.

## C. Scheduler

The config is preconfigured for **LSF** with the allocation `-P acc_paul_oreilly`
(the original study's account), set per-process in the `process { withName: ... }`
blocks of `nextflow.config`.

- [ ] Replace `acc_paul_oreilly` with your own allocation (find/replace in the config)
- [ ] For **SLURM**: change `executor = 'lsf'` → `'slurm'`, `queue`, and
  `clusterOptions` (e.g. `--account=...`) in each block
- [ ] For **local testing**: set `executor = 'local'` and `maxForks` (and reduce
  `num_random_sets`)

## D. Generate the null gene sets (GSR setup)

- [ ] `generate_random_gmts = true` (default) — the pipeline builds
  `${outdir}/randomized_gene_sets/random_birewire/` and `.../random_keeppathsize/`
  (`GeneSet.random{1..N}.gmt`) before any enrichment runs. **~12 h** for N=1000.
- [ ] (Optional) pre-generate with `scripts/core/generate_*_gmts.R`, set
  `generate_random_gmts = false`, and point `gmt_dirs` at them.

## E. Run

```bash
nextflow run main.nf -resume                 # full pipeline
nextflow run main.nf --run_magma true \
  --run_prset false --run_gsamixer false \
  --run_pascalx false -resume                # MAGMA-only smoke test
```

- [ ] Core GSR output: `{trait}_{tool}_{method}_empirical_pvalues.txt`
  (see [METHODS_GSR.md](METHODS_GSR.md))
- [ ] Validation output under `${outdir}/{opentargets,tissue,malacards,dorothea}_correlation/{trait}/`
  (see [VALIDATION.md](VALIDATION.md))

## F. Known study-specific behaviors to be aware of

- **PRSet trait exclusion:** SCZ, IBD, and AD are skipped for PRSet
  ([`../main.nf`](../main.nf) — the `prset_dedup_data` filter). Remove traits from
  that filter if you have suitable PRSet inputs for them.
- **PRSet `--set-perm`:** real model uses `2`, random models use `10`
  ([`../modules/prset/prset.nf`](../modules/prset/prset.nf)). Confirm this matches
  your intent before reusing.
- **OpenTargets Top-N:** `OT_correlation_stats.R` evaluates the top `500` pathways
  (hardcoded `top_ns <- c(500)`); edit that line to change the cutoff.
- **Validation prep scripts** (`dorothea.r`, `malacards_convert_gene_ID.R`) are
  one-off helpers with their own hardcoded paths — edit before running; they are not
  part of the Nextflow workflow.
- **MalaCards** data is license-restricted and cannot be redistributed.
