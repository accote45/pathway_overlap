# Data Directory — Required Inputs Manifest

The pipeline expects a number of reference, genotype, and validation files that are
**not committed to this repository** (they are large and/or license-restricted).
This file lists every input a new user must obtain, the expected format, and where
to get it. After downloading, point the matching `params.*` entry in
[`../nextflow.config`](../nextflow.config) at the file's location.

> Every absolute path in `nextflow.config`, `json_files/GWAS_input.json`, and this
> manifest is a `/path/to/...` placeholder. Replace them with your own paths.

---

## 1. Core gene sets

| File | `params` key | Format | Source |
|------|--------------|--------|--------|
| `c2.all.v2023.2.Hs.symbols.gmt_filtered.txt` | `geneset_real` | GMT (tab-delimited: name, description, gene symbols…) | [MSigDB C2](https://www.gsea-msigdb.org/gsea/msigdb/) v2023.2.Hs, symbols. "filtered" = restricted to genes present in the gene-coordinate reference (see note below). |

**Filtering note:** the `_filtered` GMT keeps only genes that exist in your gene
coordinate file (`gene_files.msigdbgenes`). Filter the raw MSigDB GMT to that gene
universe before running, so the real and randomized gene sets share one gene space.

Place at: `data/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt` (or update `params.geneset_real`).

## 2. Randomized (null) gene sets

Generated automatically by the pipeline when `params.generate_random_gmts = true`
(default). They are written to `${outdir}/randomized_gene_sets/random_birewire/`
and `.../random_keeppathsize/` as `GeneSet.random{1..N}.gmt`.

- **BiReWire** — degree-preserving rewiring of the gene×pathway bipartite graph
  (preserves both pathway sizes and per-gene membership frequency).
- **KeepPathSize** — random gene sets that preserve only pathway sizes.

To pre-generate manually instead, see [`../scripts/core/generate_birewire_gmts.R`](../scripts/core/generate_birewire_gmts.R)
and `generate_keeppathsize_gmts.R`, then set `generate_random_gmts = false` and point
`params.gmt_dirs` at the output directories.

## 3. Gene coordinate references

| File | `params` key | Format | Source |
|------|--------------|--------|--------|
| `msigdbgenes.regions` | `gene_files.msigdbgenes`, `pascalx_genome_annot` | MAGMA gene-loc (gene, chr, start, stop, strand, symbol) | Derived from gene annotation (e.g. NCBI/Ensembl) restricted to MSigDB genes |
| `proteincoding.regions` | `gene_files.realgenes` | MAGMA gene-loc | All protein-coding genes |
| `Homo_sapiens.GRCh37.75.gtf.gz` | `gtf_reference` | Ensembl GTF (GRCh37 build 75) | [Ensembl GRCh37 archive](https://grch37.ensembl.org/) — used by both MAGMA gene mapping and PRSet `--gtf` |

## 4. LD reference panels

| File | `params` key | Tool | Source |
|------|--------------|------|--------|
| `g1000_eur.{bed,bim,fam}` | `bfile` | MAGMA | 1000 Genomes EUR, e.g. [MAGMA reference data](https://ctg.cncr.nl/software/magma) |
| `EUR.1KG.GRCh37.*` | `pascalx_ref_panel` | PascalX | 1000 Genomes EUR (GRCh37) PLINK panel |
| `1000G.EUR.QC.@.{bim,bin,annot.gz}` | `mixer_ref_*` | GSA-MiXeR | Distributed with [GSA-MiXeR](https://github.com/precimed/gsa-mixer) (`@` = chromosome) |

> **PascalX panel placement (important):** unlike the other rows, `pascalx_ref_panel`
> is a **container-internal** path, not a host path. Place the `EUR.1KG.GRCh37.*`
> files in `data/pascalx_reference/` on the host — `nextflow.config` bind-mounts that
> directory to `/pascalx_ref` inside the container (see the `singularity.runOptions`
> block), and `pascalx_ref_panel` defaults to `/pascalx_ref/EUR.1KG.GRCh37`. Set the
> host location by changing the bind mount, not `pascalx_ref_panel`. The same applies
> to `pascalx_genome_annot` (`/data/...`, mounted from `data/`).

## 5. GWAS summary statistics

Listed per-trait in [`../json_files/GWAS_input.json`](../json_files/GWAS_input.json).
Each entry maps a trait to a summary-stat file and its **column names** (RSID, CHR,
POS, P, N, effect/other allele, effect-size column + type, SE, and optional Neff).

- Replace each `gwas_file` placeholder with your file path.
- Update the `*_col` fields to match your file's exact (case-sensitive) headers.
- `summary_statistic_type` is `beta` or `or`; `binary_target` is `T`/`F`.
- Traits used in the published study: `t2d, cad, ad, mdd, scz, ibd, breast, bmi`
  plus UK Biobank lab-value traits.

## 6. PRSet validation cohort (UK Biobank)

PRSet scores polygenic risk in a held-out cohort. Files live under `params.ukb_dir`:

| File | Purpose |
|------|---------|
| `ukb18177_chr1.22.{bed,bim,fam}` | Target genotypes (`--target`) |
| `ukb18177-qc.snplist` | QC'd SNP list (`--extract`) |
| `ukb_test_samples.txt` | Held-out test samples (`--keep`) |
| `ukb_phenofile_forprset.txt` | Phenotypes; column `${trait}_resid` per trait (`--pheno`) |

Also set `params.prsice_bin` (PRSice executable) and `params.prset_background`
(gene-coordinate background file, with the `:gene` suffix). This is a
**PRSice-format** gene background (consumed as `--background <file>:gene`) and is a
**separate file** from the MAGMA gene-loc `msigdbgenes.regions` in section 3 — don't
reuse one for the other.

> **Access:** UK Biobank individual-level data is restricted and cannot be
> redistributed. A new user needs their own approved UKB application + the prep
> scripts in [`../scripts/ukb_data_preparation/`](../scripts/ukb_data_preparation/).

## 7. Tool containers

| File | `params` key | Source |
|------|--------------|--------|
| `gsa-mixer.sif` | `mixer_sif` | Build from [GSA-MiXeR](https://github.com/precimed/gsa-mixer) |
| `pascalx_singularity.sif` | `pascalx_sif` | See [`../scripts/tool_specific/pascalx/to_create_sif.sh`](../scripts/tool_specific/pascalx/to_create_sif.sh) |
| GSA-MiXeR annotation CSVs | `mixer_go_base`, `mixer_go_full`, `mixer_go_full_test` | Distributed with GSA-MiXeR reference data |

---

## Validation benchmark databases

These power the external-benchmark comparison (the second core analysis).

### OpenTargets — `params.opentargets_json_dir`
Directory of `associationByDatatypeDirect` JSON files (one per disease EFO/MONDO ID).
Download from the [OpenTargets data downloads](https://platform.opentargets.org/downloads).
Trait → disease-ID mapping is in [`../scripts/validation/opentargets/OT_correlation_stats.R`](../scripts/validation/opentargets/OT_correlation_stats.R)
(`trait_mapping`). **Record the OpenTargets release version** for reproducibility.

### GTEx tissue specificity — `params.tissue_expression_data`
Directory of per-tissue specificity tables (gene `Name` column + expression/specificity
values). Derived from [GTEx](https://gtexportal.org/home/datasets). Note the GTEx
version and how specificity was computed.

### MalaCards — `params.malacards_path`
Per-trait disease-gene CSVs (`{trait}_genes.csv`, columns `ENSEMBL,score`).
⚠️ **License-restricted: MalaCards data cannot be redistributed.** A new user must
obtain their own access from [MalaCards](https://www.malacards.org/). The ID-conversion
helper is [`../scripts/validation/malacards/malacards_convert_gene_ID.R`](../scripts/validation/malacards/malacards_convert_gene_ID.R)
(edit its hardcoded input/output paths before running — it is a one-off prep script,
not part of the Nextflow workflow).

### DoRothEA pairwise scores — `params.dorothea_scores_path`
A **precomputed** `dorothea_pairwise_scores.csv` (`pathway1,pathway2,score`) — this is
**not** the standard DoRothEA TF-target database. It is generated by
[`../scripts/validation/dorothea/dorothea.r`](../scripts/validation/dorothea/dorothea.r)
(a one-off prep script with a hardcoded GMT path — edit it before running). Build this
file once, then point `params.dorothea_scores_path` at it.
