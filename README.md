# Pathway Overlap Pipeline – Comprehensive Pathway Enrichment Analysis

A Nextflow DSL2 pipeline for pathway enrichment analysis across GWAS traits with empirical p-value computation via null models. The pipeline coordinates **three enrichment tools** (MAGMA, PRSet, GSA-MiXeR) and validates results against multiple biological databases.

**Core workflow:**
1. **Real pathway enrichment** → Run enrichment tools on actual gene sets  
2. **Null model generation** → Run same tools on randomized gene sets (BiReWire, KeepPathSize)
3. **Empirical validation** → Compare real vs null distributions to compute standardized effect sizes
4. **Biological validation** → Cross-reference top pathways with OpenTargets, MalaCards, tissue specificity, and DoRothEA

---

## Key Features

- **Three enrichment methods**: MAGMA gene-set analysis, PRSet polygenic risk scores, GSA-MiXeR mixed-model analysis
- **Dual randomization strategies**: BiReWire (degree-preserving) and KeepPathSize (size-preserving) null models
- **Empirical statistics**: Standardized effect sizes and empirical p-values from 1000+ permutations per method
- **Multi-database validation**: OpenTargets gene-disease evidence, MalaCards disease associations, tissue expression specificity, DoRothEA regulatory networks
- **False positive rate analysis**: Assess method-specific FPR across randomized gene sets
- **Flexible GWAS input**: JSON-configured column mapping for heterogeneous GWAS formats
- **Trait filtering**: Automatic filtering for tool compatibility (e.g., SCZ excluded from PRSet)

---

## Requirements

- **Nextflow** (>= 22.10) and Java 11+
- **R** (>= 4.1) with packages: `tidyverse`, `data.table`, `plyr`, `jsonlite`, `GSA`, `MatchIt`, `ggplot2`, `patchwork`, `gridExtra`
- **MAGMA** v1.10 (if `--run_magma true`) and reference PLINK bfile (e.g., 1000G EUR)  
- **PRSet/PRSice** (if `--run_prset true`) with compatible GWAS formats
- **GSA-MiXeR** Singularity container (if `--run_gsamixer true`)
- **Scheduler**: LSF (default) or adapt `nextflow.config` for your environment
- **Data dependencies**:
  - GMT gene set files (real + randomized via BiReWire/KeepPathSize)
  - OpenTargets association JSONs (`associationByDatatypeDirect`)
  - Tissue expression specificity data (GTEx-derived)
  - MalaCards disease-gene association files
  - DoRothEA regulatory interaction scores

---

## Repository Structure

```
├── main.nf                     # Main pipeline workflow
├── nextflow.config            # Configuration and resource settings (LSF default)
├── json_files/
│   └── GWAS_input.json       # GWAS trait configuration example
├── modules/                   # Nextflow process modules
│   ├── magma/                # MAGMA gene-set analysis
│   ├── prset/                # PRSet polygenic risk scoring  
│   ├── gsamixer/             # GSA-MiXeR enrichment
│   ├── empirical_pval/       # Empirical p-value calculation
│   ├── opentargets/          # OpenTargets validation
│   ├── malacards/            # MalaCards validation
│   ├── tissuespecificity/    # Tissue expression analysis
│   ├── dorothea/             # DoRothEA regulatory networks
│   └── fpr/                  # False positive rate analysis
├── scripts/                   # R analysis scripts
│   ├── calc_empirical.r      # Core empirical p-value computation
│   ├── *_correlation_*.R     # Validation correlation analyses
│   └── gsamixer/             # GSA-MiXeR helper scripts
└── README.md

---

## Configuration

### Essential Parameters (`nextflow.config`)

**Data Paths:**
```groovy
params.traits_config = "/path/to/GWAS_input.json"     // GWAS trait configuration
params.outdir = "/path/to/results"                   // Output directory
params.geneset_real = "/path/to/real_pathways.gmt"   // Real gene sets
params.gmt_dirs = [                                  // Randomized gene sets
  birewire: "/path/to/random_birewire_10k",
  keeppathsize: "/path/to/random_keeppathsize_10k"  
]
```

**Tool References:**
```groovy
params.bfile = "/path/to/plink_reference"            // MAGMA PLINK reference
params.ukb_dir = "/path/to/ukb_data"                // PRSet UKB reference data
params.mixer_sif = "/path/to/gsa-mixer.sif"         // GSA-MiXeR container
params.opentargets_json_dir = "/path/to/ot_jsons"   // OpenTargets data
```

**Workflow Control:**
```groovy
params.run_magma = true                    // Enable MAGMA analysis
params.run_prset = true                   // Enable PRSet analysis  
params.run_gsamixer = true                // Enable GSA-MiXeR analysis
params.run_empirical = true               // Calculate empirical p-values
params.run_ot_correlation = true          // OpenTargets validation
params.run_tissue_correlation = true      // Tissue specificity analysis
params.run_malacards_correlation = true   // MalaCards validation
params.run_dorothea_correlation = true    // DoRothEA regulatory analysis
params.run_fpr_analysis = false          // False positive rate analysis
```

**Method Configuration:**
```groovy
params.randomization_methods = ['birewire', 'keeppathsize']
params.num_random_sets = 1000                       // Permutations per method
params.opentargets_supported_traits = ["t2d", "cad", "mdd", "scz", "ibd", "breast", "ad"]
params.malacards_traits = 'bmi,cad,t2d,mdd,ad,scz,ibd,breast'  // Comma-separated
```

**⚠️ Important Trait Filtering:**
- **PRSet excludes SCZ**: Hardcoded filter in `main.nf:184` removes schizophrenia from PRSet analysis
- **OpenTargets whitelist**: Only traits in `opentargets_supported_traits` will be processed
- **MalaCards subset**: Limited to traits specified in `malacards_traits` parameter

---

## GWAS Input Format

Configure traits in `json_files/GWAS_input.json` with column mappings:

```json
[
  {
    "trait": "t2d",
    "gwas_file": "/path/to/t2d_gwas.txt",
    "rsid_col": "SNP",
    "chr_col": "CHR", 
    "pos_col": "POS",
    "pval_col": "P",
    "n_col": "N",
    "neff_col": "Neff",
    "effect_allele": "A1",
    "other_allele": "A2", 
    "binary_target": "T",
    "summary_statistic_name": "BETA",
    "summary_statistic_type": "beta",
    "se_col": "SE"
  }
]
```

**Column Mapping Notes:**
- Ensure column names exactly match your GWAS file headers
- `neff_col` optional (falls back to `n_col` if missing)
- `binary_target`: "T" for case-control, "F" for quantitative traits
- **PRSet automatically excludes SCZ trait** (hardcoded filter)

---

## Usage

**Basic execution:**
```bash
nextflow run main.nf -resume
```

**Custom configuration:**
```bash
nextflow run main.nf \
  --traits_config ./my_gwas_config.json \
  --outdir /scratch/results \
  --num_random_sets 1000 \
  --run_magma true --run_prset true --run_gsamixer true \
  --run_empirical true --run_ot_correlation true \
  -resume
```

**Subset analysis (MAGMA + empirical only):**
```bash
nextflow run main.nf \
  --run_magma true --run_prset false --run_gsamixer false \
  --run_empirical true --run_ot_correlation false \
  --run_tissue_correlation false --run_malacards_correlation false \
  -resume
```

---

## Pipeline Workflow

### 1. GWAS Preprocessing
- **Column harmonization** via JSON configuration
- **SNP deduplication** for PRSet (removes ALL instances of duplicate RSIDs)
- **Format preparation** for each enrichment tool
- **Trait filtering**: SCZ automatically excluded from PRSet analysis

### 2. Enrichment Analysis (Real Gene Sets)
- **MAGMA**: SNP → Gene → Pathway enrichment with PLINK reference
- **PRSet**: GWAS → PRS → Pathway scoring with UKB validation cohort  
- **GSA-MiXeR**: GWAS → Chr-split → Baseline + Full model with Singularity container

### 3. Null Model Generation
- **BiReWire**: Degree-preserving network randomization (genes keep same pathway count, different assignments)
- **KeepPathSize**: Size-preserving randomization (pathways keep size, random gene sampling)
- **1000 permutations** per trait-method combination
- **Parallel execution** across randomized GMT files

### 4. Empirical Validation  
- **Empirical p-values**: Real results vs null distributions
- **Standardized effect sizes**: `(real_beta - null_mean) / null_sd`
- **Tool-agnostic calculation** via `scripts/calc_empirical.r`
- **Output**: `{trait}_{tool}_{method}_empirical_pvalues.txt`

### 5. Biological Validation

**OpenTargets Analysis:**
- Size-matched pathway comparison against gene-disease associations
- **⚠️ Restricted to `params.opentargets_supported_traits` whitelist**
- Uses `associationByDatatypeDirect` JSON files 
- Outputs: pathway rankings, advantage metrics, correlation summaries

**Tissue Specificity Analysis:**
- Correlates pathway rankings with tissue expression specificity
- Independent validation across all traits
- Uses GTEx-derived tissue expression data
- Outputs: correlation summaries, tissue-specific plots

**MalaCards Integration:**  
- Disease-gene association validation
- **⚠️ Restricted to `params.malacards_traits` list**
- Correlation analysis
- Outputs: pathway scores, correlation summaries

**DoRothEA Regulatory Networks:**
- Transcription factor regulatory pathway analysis  
- Correlates pathway rankings with regulatory interaction scores
- Outputs: regulatory relevance summaries

### 6. False Positive Rate Analysis (Optional)
- **Method-specific FPR** calculation across randomized gene sets
- **Pathway size effects** on false discovery rates
- **Quality control** for enrichment methods

---

## Key Outputs

```
results/
├── {tool}_real/{trait}/                    # Real enrichment results
├── {tool}_random/{method}/{trait}/         # Null model results (1000 per trait)  
├── empirical_pvalues/{tool}/{trait}/       # Empirical statistics
├── opentargets_correlation/{trait}/        # OpenTargets validation
├── tissue_correlation/{trait}/             # Tissue specificity analysis  
├── malacards_correlation/{trait}/          # MalaCards validation
├── dorothea_correlation/{trait}/           # DoRothEA regulatory analysis
└── fpr_analysis/{tool}/{method}/           # False positive rates
```

**Critical Output Files:**
- **Empirical results**: `{trait}_{tool}_{method}_empirical_pvalues.txt` 
  - Columns: `pathway_name`, `p_value`, `beta_value`, `empirical_pval`, `std_effect_size`
- **OpenTargets**: `{trait}_{tool}_rank_correlation_summary.csv`
- **Tissue**: `{trait}_{tool}_tissue_correlation_summary.csv` 
- **Delta-rank**: `{trait}_{tool}_delta_rank_{validation}_correlation_summary.csv`

---

## Validation Database Details

### OpenTargets
- **Purpose**: Gene-disease association evidence validation
- **Data**: `associationByDatatypeDirect` JSON files (animal_model, known_drug, literature)
- **Method**: Size-matched pathway comparison with evidence density scoring
- **⚠️ Restriction**: Limited to `params.opentargets_supported_traits` whitelist

### MalaCards  
- **Purpose**: Human disease-gene association validation
- **Data**: Disease-specific gene score files with ENSEMBL IDs
- **Method**: Pathway-level score aggregation and ranking correlation
- **⚠️ Coverage**: Restricted to traits in `params.malacards_traits`

### Tissue Specificity
- **Purpose**: Tissue-relevant pathway identification  
- **Data**: GTEx-derived tissue expression specificity matrices
- **Method**: Spearman correlation between pathway ranks and tissue specificity
- **Output**: Multi-tissue correlation profiles per trait

### DoRothEA
- **Purpose**: Transcriptional regulatory pathway analysis
- **Data**: Curated transcription factor-target interaction scores
- **Method**: Pathway regulatory activity scoring and ranking correlation
- **Focus**: Regulatory mechanism validation

---

## Method-Specific Notes

### MAGMA
- Requires PLINK reference files (`params.bfile`)
- Uses 35kb gene windows for SNP-gene mapping
- Outputs: `*.genes.raw` (gene-level) → `*.gsa.out` (pathway-level)

### PRSet  
- **⚠️ Excludes SCZ trait** (hardcoded in `main.nf:184`)
- Requires UKB validation cohort data 
- Uses deduplicated GWAS (removes ALL duplicate SNP instances)
- Outputs: `*.summary` (pathway results), `*.prsice` (PRS details)

### GSA-MiXeR
- Uses Singularity container for execution
- Requires baseline annotations and reference files
- Two-stage analysis: baseline model → full enrichment model  
- High memory requirements (40GB for full models)
- Outputs: `*.go_test_enrich.csv` (enrichment results), `*.json` (model parameters)

---

## Troubleshooting

### Common Issues

**Missing Random Files:**
```bash
# Check randomization file counts (should be 1000 per trait-method)
find results/{tool}_random/{method} -name "*out" | wc -l
```

**Column Mapping Errors:**
```bash
# Verify GWAS column names match JSON configuration  
head -1 /path/to/gwas.txt | tr '\t' '\n' | nl
```

**Memory Failures:**
- GSA-MiXeR full models require 40GB+ memory
- Empirical calculation needs 50GB for large trait sets
- Adjust memory limits in `nextflow.config` process blocks

**Channel Grouping Issues:**
- Random results must complete BEFORE empirical calculation  
- Each trait-method needs exactly 1000 random files for valid statistics
- Use `.view()` to debug channel contents: `channel.view { "DEBUG: $it" }`

### Scheduler Configuration
- **Default**: LSF with project codes and resource specifications
- **Adaptation**: Remove LSF-specific blocks, set global executor in `nextflow.config`
- **Resource tuning**: Adjust memory/time limits per process requirements

### Tool-Specific Debugging

**MAGMA Issues:**
- Verify `module load magma_gwas/1.10` availability
- Check PLINK reference file integrity
- Ensure gene annotation files are compatible

**PRSet Issues:** 
- Confirm UKB reference data paths and permissions
- Validate SNP deduplication doesn't remove too many variants
- Check phenotype file column names match trait specifications
- **Note**: SCZ trait automatically filtered out

**GSA-MiXeR Issues:**
- Verify Singularity container accessibility 
- Check baseline annotation file formats
- Monitor memory usage during full model execution

### Validation-Specific Issues

**OpenTargets Failures:**
- Verify `associationByDatatypeDirect` JSON files exist
- Check trait is in `params.opentargets_supported_traits` whitelist
- Ensure gene ID mapping between pathways and OpenTargets data

**Tissue Specificity Errors:**
- Verify `params.tissue_expression_data` path and format
- Check CSV has 'Name' column with gene symbols
- Ensure tissue column names match expected format

**MalaCards Issues:**
- Confirm trait is in `params.malacards_traits` list
- Check MalaCards data files have correct ENSEMBL ID format
- Verify disease-trait name mapping consistency

---

## Citation

If you use this pipeline, please cite:
- **MAGMA**: de Leeuw, C.A. et al. (2015) MAGMA: Generalized gene-set analysis of GWAS data. PLOS Computational Biology 11(4): e1004219
- **PRSice**: Choi, S.W. et al. (2019) PRSice-2: Polygenic Risk Score software for biobank-scale data. GigaScience 8(7): giz082  
- **GSA-MiXeR**: Frei, O. et al. (2019) Bivariate causal mixture model quantifies polygenic overlap between complex traits beyond genetic correlation. Nature Communications 10: 2417
- **OpenTargets**: Ochoa, D. et al. (2021) Open Targets Platform: supporting systematic drug-target identification and prioritisation. Nucleic Acids Research 49(D1): D1302-D1310

## License

MIT License - see `LICENSE` file for details.