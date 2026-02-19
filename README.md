# Pathway Overlap Pipeline – Comprehensive Pathway Enrichment Analysis with Empirical Validation

A **Nextflow DSL2 pipeline** for pathway enrichment analysis across GWAS traits with rigorous empirical p-value computation via null models. This pipeline coordinates **three complementary enrichment tools** (MAGMA, PRSet, GSA-MiXeR) and validates results against multiple databases (OpenTargets, MalaCards, GTEx tissue specificity, DoRothEA regulatory networks).

## Core Methodology

**Problem**: Standard GWAS pathway enrichment methods do not account for the sharing of genes across pathways.

**Solution**: This pipeline adjusts real enrichment statistics for empirical null distributions generated through two gene randomization strategies:

1. **Real pathway enrichment** → Run enrichment tools on actual curated gene sets  
2. **Null model generation** → Run identical analyses on K randomized gene sets per method
3. **Empirical validation** → Calculate standardized effect sizes: `(real_beta - null_mean) / null_sd`
4. **Biological validation** → Cross-reference top pathways with disease-relevant databases

**Key Innovation**: The pipeline uses **two complementary null models**:
- **GSR**: Preserves gene frequency (genes maintain degree), randomizes pathway membership
- **KeepPathSize**: Preserves pathway sizes, randomly samples genes with equal probability. This approach is in line with standard gene randomization and is provided for comparison to GSR.

This dual approach distinguishes between **pathway size effects** and **pathway overlap biases**.

---

## Key Features

✅ **Three pathway methods**: MAGMA, PRSet (polygenic risk scores), GSA-MiXeR
✅ **Dual randomization strategies**: GSR (pathway overlap preserving) + KeepPathSize (pathway size-preserving)  
✅ **Robust empirical statistics**: K+ permutations per trait-method combination  
✅ **Multi-database validation**: OpenTargets, MalaCards, GTEx tissues, DoRothEA networks  
✅ **False positive rate analysis**: Method-specific FPR quantification  
✅ **Flexible GWAS input**: JSON-configured column mapping for heterogeneous formats  
✅ **Production-ready**: LSF scheduler integration, containerized tools, automatic resume  

---
## Requirements

### Software Dependencies

- **Nextflow** ≥ 22.10 (DSL2 support required)
- **Java** 11+ (for Nextflow execution)
- **R** ≥ 4.1 with packages:
  ```r
  install.packages(c("tidyverse", "data.table", "plyr", "jsonlite", 
                     "GSA", "MatchIt", "ggplot2", "patchwork", 
                     "gridExtra", "RColorBrewer", "pheatmap"))
  ```
  
- **MAGMA** v1.10 (if `params.run_magma = true`)
  - Download: [https://ctg.cncr.nl/software/magma](https://ctg.cncr.nl/software/magma)
  - Requires `module load` capability or PATH configuration
  
- **PRSice/PRSet** v2.3+ (if `params.run_prset = true`)
  - Download: [https://www.prsice.info/](https://www.prsice.info/)
  - Requires UKB validation cohort data (see Data Dependencies)
  
- **GSA-MiXeR** Singularity container (if `params.run_gsamixer = true`)
  - Container: `gsa-mixer.sif` (contact authors or build from source)
  - Requires Singularity/Apptainer 3.0+

### Scheduler Requirements

**Default configuration**: LSF batch system with project codes
- If using a different scheduler (SLURM, SGE, PBS), modify executor blocks in `nextflow.config`
- Resource requirements per process documented in config file

### Data Dependencies

#### Required Data Files

1. **Reference gene sets** (real pathways):
   ```
   # MSigDB C2 canonical pathways (recommended)
   wget http://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.2.Hs/c2.all.v2023.2.Hs.symbols.gmt
   # Filter to protein-coding genes (script provided in repository)
   ```

2. **Randomized gene sets** (BiReWire + KeepPathSize):
   - **⚠️ CRITICAL**: BiReWire GMTs must be **pre-generated** (12h runtime for 1000 permutations)
   - KeepPathSize GMTs can be generated on-the-fly (2-4h runtime)
   - Expected structure:
     ```
     data/randomized_gene_sets/
     ├── random_birewire/
     │   ├── GeneSet.random1.gmt
     │   ├── GeneSet.random2.gmt
     │   └── ... (up to GeneSet.random(k).gmt)
     └── random_keeppathsize/
         ├── GeneSet.random1.gmt
         └── ... (same structure)
     ```

3. **MAGMA reference panel** (1000 Genomes EUR):
   ```bash
   # Download from MAGMA website: https://cncr.nl/research/magma/
   ```

4. **PRSet UKB reference data**:
   - Training genotypes: `ukb18177_eur_autosomes.{bed,bim,fam}`
   - QC SNP list: `ukb18177-qc.snplist`
   - Phenotype file: `ukb_phenofile_forprset.txt` (residualized phenotypes)
   - Train/test splits: `ukb_train_samples.txt`, `ukb_test_samples.txt`
   - GTF gene annotation: `Homo_sapiens.GRCh37.75.gtf.gz`

5. **GSA-MiXeR reference files**:
   - Reference panel: `1000G.EUR.QC.@.bim` (chromosome-split format)
   - Baseline annotations: `baseline_v1.1_genes_chr@.annot.gz`
   - LD scores: Provided with container

#### Validation Database Files

6. **OpenTargets gene-disease associations** (Required for `run_ot_correlation`):
   ```bash
   # Download associationByDatatypeDirect JSONs
   # These are large files (~50-100GB per disease)
   # Required traits: t2d, cad, ad, mdd, scz, ibd, breast
   # Example structure:
   # opentargets_data/
   # ├── EFO_0000311.json  # Schizophrenia (SCZ)
   # ├── EFO_0001645.json  # Coronary artery disease (CAD)
   # └── ...
   ```
   **Disease ID mapping**:
   ```
   t2d     → EFO_0001360
   cad     → EFO_0001645
   ad      → EFO_0000249
   mdd     → EFO_0003761
   scz     → EFO_0000311
   ibd     → EFO_0000555
   breast  → MONDO_0007254
   ```

7. **GTEx tissue specificity data** (Required for `run_tissue_correlation`):
   - Tissue-specific gene expression matrices
   - Format: CSV with 'Name' column (gene symbols) + tissue columns
   - Download: [GTEx Portal](https://gtexportal.org/home/datasets)
   - Expected path: `${params.tissue_expression_data}`

8. **MalaCards disease-gene associations** (Required for `run_malacards_correlation`):
   - Disease-specific CSV files with gene scores
   - Required columns: `ENSEMBL`, `score`
   - Trait restriction: `bmi, cad, t2d, mdd, ad, scz, ibd, breast`
   - Path pattern: `${params.malacards_path}/{trait}_genes.csv`

9. **DoRothEA regulatory interaction scores** (Required for `run_dorothea_correlation`):
   - **⚠️ NOT the standard DoRothEA TF-target database**
   - Requires **pairwise pathway scores**: `dorothea_pairwise_scores.csv`
   - Format: `pathway1, pathway2, score` (pre-computed regulatory similarity)
   - Contact authors for access or generate using DoRothEA confidence levels A-C

---
## Quick Start

### Step 1: Clone Repository
```bash
git clone https://github.com/accote45/pathway_overlap.git
cd pathway_overlap
```

### Step 2: Configure GWAS Inputs
Edit `json_files/GWAS_input.json` with your trait-specific column mappings:

```json
[
  {
    "trait": "t2d",
    "gwas_file": "/full/path/to/t2d_gwas.txt.gz",
    "rsid_col": "SNP",
    "chr_col": "CHR",
    "pos_col": "BP",
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

**Critical JSON fields**:
- `trait`: Short identifier (lowercase, used in filenames)
- `gwas_file`: **Absolute path** to GWAS summary statistics
- `*_col`: Column header names in your GWAS file (case-sensitive)
- `neff_col`: Effective sample size (falls back to `n_col` if missing)
- `binary_target`: "T" for case-control, "F" for quantitative
- `summary_statistic_type`: "beta", "or" (odds ratio), or "z"

**Column verification tip**:
```bash
# Show column names with numbers
head -1 /path/to/gwas.txt | tr '\t' '\n' | nl
```

### Step 3: Update Configuration Paths
Edit `nextflow.config` to point to your data:

```groovy
params {
  // Input/output
  traits_config = "${projectDir}/json_files/GWAS_input.json"
  outdir = "/scratch/user/pathway_results"
  
  // Reference data
  geneset_real = "/data/msigdb/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt"
  bfile = "/data/reference/g1000_eur"  // MAGMA reference
  ukb_dir = "/data/ukb/"               // PRSet reference
  mixer_sif = "/containers/gsa-mixer.sif"
  
  // Randomized gene sets
  gmt_dirs = [
    birewire: "/data/randomized_gene_sets/random_birewire",
    keeppathsize: "/data/randomized_gene_sets/random_keeppathsize"
  ]
  
  // Validation databases
  opentargets_json_dir = "/data/opentargets/associationByDatatypeDirect"
  tissue_expression_data = "/data/gtex/tissue_specificity"
  malacards_path = "/data/malacards"
  dorothea_scores_path = "/data/dorothea_pairwise_scores.csv"
}
```

### Step 4: Run Pipeline
```bash
# Full analysis (all tools + validation)
nextflow run main.nf -resume

# MAGMA only (faster for testing)
nextflow run main.nf \
  --run_magma true \
  --run_prset false \
  --run_gsamixer false \
  --run_empirical true \
  --run_ot_correlation false \
  -resume

# Custom output directory
nextflow run main.nf \
  --outdir /scratch/my_analysis \
  --num_random_sets 500 \
  -resume
```

**Execution tips**:
- Always use `-resume` to restart from last successful step
- Check `work/` directory for intermediate files if debugging
- Use `-with-report report.html` for execution statistics

---

## Pipeline Architecture

### Workflow Stages

```
┌─────────────────────────────────────────────────────────────┐
│ 1. RANDOMIZATION SETUP (Must Complete First)               │
│    ├─ BiReWire GMT generation (12h, pre-generated)         │
│    └─ KeepPathSize GMT generation (2-4h, can run on-the-fly)│
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ 2. REAL ENRICHMENT ANALYSIS (Parallel by Trait)            │
│    ├─ MAGMA: SNP→Gene→Pathway                              │
│    ├─ PRSet: GWAS→PRS→Pathway (excludes SCZ)               │
│    └─ GSA-MiXeR: GWAS→Chr-split→Baseline→Full              │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ 3. NULL MODEL ENRICHMENT (1000 permutations per trait)     │
│    ├─ MAGMA: 1000 × random_birewire + 1000 × keeppathsize  │
│    ├─ PRSet: Same structure                                 │
│    └─ GSA-MiXeR: Same structure                             │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ 4. EMPIRICAL STATISTICS CALCULATION                         │
│    ├─ Group 1000 random results per trait-method           │
│    ├─ Calculate: empirical_pval = rank(real) / 1001        │
│    └─ Calculate: std_effect = (real - μ_null) / σ_null     │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ 5. BIOLOGICAL VALIDATION (Parallel)                         │
│    ├─ OpenTargets: Gene-disease association                 │
│    ├─ MalaCards: Disease-gene association                   │
│    ├─ Tissue: GTEx tissue-specific expression               │
│    └─ DoRothEA: Regulatory pathway interaction              │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ 6. FALSE POSITIVE RATE ANALYSIS (Optional)                  │
│    └─ Calculate FPR across 1000 null permutations          │
└─────────────────────────────────────────────────────────────┘
```

### Critical Channel Architecture

The pipeline uses **grouped channels** to coordinate real/random results:

```nextflow
// 1. Generate 1000 random results per trait-method
random_sets_inputs = gene_results
    .combine(Channel.from(1..1000))  // Permutation numbers
    .map { trait, gene_result, perm -> 
        [trait, gene_result, "birewire", perm]
    }

// 2. Run enrichment on each random set
random_results = run_random_sets(random_sets_inputs)

// 3. Group by trait-method (collects all 1000 files)
random_grouped = random_results.groupTuple(by: [0, 2])
// Output: [trait, [file1, file2, ...1000], method]

// 4. Join real results with grouped random results
empirical_inputs = real_results
    .combine(random_grouped, by: [0, 2])
    .map { trait, real_file, method, random_files ->
        def random_dir = "${params.outdir}/.../random/${method}/${trait}"
        tuple(trait, "tool", real_file, random_dir)
    }

// 5. Calculate empirical statistics
empirical_results = calc_empirical_pvalues(empirical_inputs)
```

**Key pattern**: Random results MUST complete before empirical calculation starts.

### File Format Specifications

#### Tool Output Formats

**MAGMA** (`*.gsa.out`):
```
VARIABLE  TYPE  NGENES  BETA  BETA_STD  SE  P  FULL_NAME
pathway1  set   150     0.45  0.023     0.12  0.0001  Sample Pathway 1
pathway2  set   89      -0.12 0.008     0.09  0.18    Sample Pathway 2
```

**PRSet** (`*.summary`):
```
Set              Threshold  PRS.R2   Full.R2  Coefficient  P   Competitive.P
KEGG_PATHWAY_X   1         0.0123   0.045    0.0056       0.03  0.001
REACTOME_Y       1         0.0089   0.042    0.0043       0.08  0.023
```

**GSA-MiXeR** (`*_go_test_enrich.csv`):
```
GO,gene_set,trait1_enrich_estimate,trait1_enrich_se,trait1_pval
GO:0006915,apoptosis,0.234,0.045,0.0023
GO:0007049,cell_cycle,-0.123,0.067,0.067
```
**Note**: GSA-MiXeR p-values are NOT used—only enrichment estimates.

#### Empirical Results Format

**Output** (`{trait}_{tool}_{method}_empirical_pvalues.txt`):
```
pathway_name         p_value    beta_value  empirical_pval  std_effect_size  num_genes
KEGG_GLYCOLYSIS      0.0001     0.456       0.001          2.34             52
REACTOME_APOPTOSIS   0.0234     -0.123      0.156         -0.89            147
```

**Column definitions**:
- `p_value`: Original tool p-value (raw, not empirical)
- `beta_value`: Effect size from enrichment tool
- `empirical_pval`: Rank-based p-value from null distribution
- `std_effect_size`: `(real_beta - mean_null_beta) / sd_null_beta`
- `num_genes`: Pathway size

**Special handling**:
- **"Base" pathway**: Filtered out (GSA-MiXeR baseline model artifact)
- **"coding_genes" pathway**: Filtered out (GSA-MiXeR universe control)
- Empty pathway names: Skipped with warning

---
## Tool-Specific Implementation Details

### MAGMA Gene-Set Analysis

**Workflow**:
1. **SNP-to-gene mapping**: 35kb upstream/downstream windows
2. **Gene-level analysis**: P-value aggregation via SNP-wise mean model
3. **Pathway enrichment**: Competitive test controlling for gene length/density

**Key parameters**:
```bash
magma \
  --bfile ${params.bfile} \
  --pval ${gwas_file} use=SNP,P ncol=N \
  --gene-annot ${trait}.genes.annot \
  --out ${trait}_genebased

magma \
  --gene-results ${trait}_genebased.genes.raw \
  --set-annot ${pathway_gmt_file} \
  --out ${trait}_pathway_enrichment
```

**Resource requirements**:
- Gene annotation: 4GB RAM, 20min
- Gene analysis: 50GB RAM, 12h (for large GWAS)
- Pathway enrichment: 4GB RAM, 1h

**Common issues**:
- Missing SNPs: Check RSID consistency between GWAS and reference
- High memory usage: Split GWAS by chromosome if needed
- Module load failure: Verify `magma_gwas/1.10` is available

### PRSet Polygenic Risk Score Analysis

**Workflow**:
1. **SNP deduplication**: Removes ALL instances of duplicate RSIDs (conservative)
2. **PRS calculation**: Clumping with r²=0.1, 1000kb window
3. **Pathway scoring**: MSigDB gene sets tested for trait association
4. **Validation**: Uses held-out UKB test set for unbiased estimates

**⚠️ CRITICAL EXCLUSION**: SCZ trait is **hardcoded to be excluded** from PRSet analysis (line 184 in `main.nf`).

**Key parameters**:
```bash
PRSice_linux \
  --a1 ${effect_allele} --a2 ${other_allele} \
  --background msigdb.genes.txt:gene \
  --base ${gwas_file} \
  --binary-target ${binary_target} \
  --clump-kb 1000 --clump-p 1.0 --clump-r2 0.1 \
  --fastscore \
  --gtf Homo_sapiens.GRCh37.75.gtf.gz \
  --keep ukb_test_samples.txt \
  --msigdb ${pathway_gmt_file} \
  --pheno ukb_phenofile_forprset.txt \
  --target ${ukb_genotypes}
```

**Resource requirements**:
- Deduplication: 4GB RAM, 10min
- PRS calculation: 15GB RAM, 2-4h
- Memory scales with: # SNPs × # pathways × # samples

**Common issues**:
- Duplicate SNPs: Pipeline removes ALL duplicates (not just highest P)
- Phenotype misalignment: Ensure trait names match between GWAS JSON and UKB pheno file
- Missing genotypes: Check `ukb18177-qc.snplist` for QC'd SNP list

**Phenotype preparation**:
```r
# Continuous traits: Residualize for covariates
pheno$trait_resid <- residuals(lm(trait ~ Age + Sex + PC1:PC10, data=pheno))

# Binary traits: Use raw 0/1 coding (covariates handled in model)
```

### GSA-MiXeR Mixed-Model Enrichment

**Workflow**:
1. **Chromosome splitting**: Split GWAS by chromosome for parallel processing
2. **Baseline model**: Estimate SNP heritability with gene-level annotations
3. **Full model**: Test pathway enrichment controlling for baseline architecture

**Key parameters**:
```bash
# Baseline model
singularity exec gsa-mixer.sif python /tools/mixer/precimed/mixer.py plsa \
  --trait1-file ${trait}.chr@.sumstats.gz \
  --out ${trait}_base \
  --use-complete-tag-indices \
  --go-file baseline.txt

# Full pathway model  
singularity exec gsa-mixer.sif python /tools/mixer/precimed/mixer.py plsa --gsa-full \
  --trait1-file ${trait}.chr@.sumstats.gz \
  --out ${trait}_full \
  --go-file full_gene.txt \
  --go-file-test full_gene_set.txt
```

**Resource requirements**:
- Baseline model: 15GB RAM, 4h
- Full model: **40GB RAM**, 24h (critical bottleneck)
- Chromosome splitting: 1GB per chromosome

**Input format** (prepared by `scripts/gsamixer/prepare_sumstats.R`):
```
RSID  CHR  POS  A1  A2  N  NEFF  Z
rs123 1    1000 A   G   50000  48000  2.34
```

**Common issues**:
- Memory killed: Increase to 50GB for traits with >10M SNPs
- Missing NEFF: Script automatically falls back to N if NEFF column absent
- Singularity permissions: Ensure `--home $(pwd):/home` bind mount works

**Special pathway handling**:
- `coding_genes`: Universe control (all protein-coding genes)
- `base`: Baseline model (filtered in empirical calculation)
- Individual genes: Tested as single-gene "pathways" for fine-mapping

---
## Randomization Strategy

### Why GSR?

### GSR Method

**Algorithm**:
```
1. Construct gene-pathway bipartite network (binary matrix)
2. Iteratively swap edges to:
   - Maintain row sums (genes belong to same number of pathways)
   - Maintain column sums (pathways keep same size)
   - Ensure no duplicate edges
3. Converge to uniform distribution over all valid graphs
```

**Properties**:
- Gene X appears in K pathways → Still appears in K random pathways
- Pathway Y has M genes → Still has M random genes
- Gene-pathway associations completely randomized

**Use case**: Tests whether pathway enrichment mayyers beyond pathway overlap structure alone.

### KeepPathSize Method

**Algorithm**:
```
1. For each pathway P with N genes:
   - Sample N genes uniformly from all available genes
   - No replacement within pathway
   - Allow genes to appear in multiple pathways
2. Pathway sizes perfectly preserved
3. Gene frequencies determined by sampling
```

**Use case**: Tests whether pathway enrichment matters beyond pathway size alone.

### GMT File Structure

**Real pathways** (`c2.all.v2023.2.Hs.symbols.gmt_filtered.txt`):
```
KEGG_GLYCOLYSIS    NA    HK1    HK2    PFKM    ALDOA    ...
REACTOME_APOPTOSIS NA    BAX    BCL2   CASP3   CASP9   ...
```

**Randomized pathways** (`GeneSet.random{1-1000}.gmt`):
```
KEGG_GLYCOLYSIS    NA    RANDOM_GENE1    RANDOM_GENE2    ...  (52 genes total)
REACTOME_APOPTOSIS NA    RANDOM_GENE3    RANDOM_GENE4    ...  (147 genes total)
```

**Critical details**:
- **Column 1**: Pathway name (preserved in randomization)
- **Column 2**: Description field (set to "NA")
- **Columns 3+**: Gene symbols (randomized content)
- **File naming**: Must follow exact pattern `GeneSet.random{N}.gmt` (N=1-1000)

### Generating Randomized GMTs

**BiReWire** (pre-generate before pipeline):
```bash
Rscript scripts/core/generate_birewire_gmts.R \
  --input data/msigdb/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt \
  --output data/randomized_gene_sets/random_birewire_10k \
  --num_random 1000

# Expect 12h runtime, 15GB RAM
# Creates GeneSet.random1.gmt through GeneSet.random1000.gmt
```

**KeepPathSize** (can generate on-the-fly):
```bash
# Enable in nextflow.config
params.generate_random_gmts = true

# OR pre-generate manually
Rscript scripts/core/generate_keeppathsize_gmts.R \
  --input data/msigdb/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt \
  --output data/randomized_gene_sets/random_keeppathsize_10k \
  --num_random 1000 \
  --background data/msigdb.genes.txt  # Universe of genes to sample

# Expect 2-4h runtime, 8GB RAM
```

---

## Biological Validation

### OpenTargets Gene-Disease Evidence

**Purpose**: Test whether top-ranked pathways are enriched for genes with known disease associations.

**Method**:
1. Load disease-specific gene evidence from OpenTargets JSON files
2. For each pathway, calculate average gene evidence score

**Data requirements**:
```bash
# Download associationByDatatypeDirect JSON files
# Required traits: t2d, cad, ad, mdd, scz, ibd, breast

opentargets_data/
├── EFO_0000311.json  # Schizophrenia (SCZ)
├── EFO_0001645.json  # Coronary artery disease (CAD)
├── EFO_0001360.json  # Type 2 diabetes (T2D)
├── EFO_0003761.json  # Major depressive disorder (MDD)
├── EFO_0000249.json  # Alzheimer's disease (AD)
├── EFO_0000555.json  # Inflammatory bowel disease (IBD)
└── MONDO_0007254.json # Breast cancer
```

**JSON structure** (simplified):
```json
{
  "data": {
    "disease": {
      "evidences": {
        "rows": [
          {
            "target": {"id": "ENSG00000123456"},
            "score": 0.85,
            "dataTypeId": "literature",
            "evidenceCount": 42
          }
        ]
      }
    }
  }
}
```

**Trait restriction**: Only analyzes traits in `params.opentargets_supported_traits`.

**Output files**:
```
{trait}_{tool}_rank_correlation_summary.csv
{trait}_{tool}_detailed_opentargets_advantage.csv
{trait}_{tool}_opentargets_comparison.pdf
```

### MalaCards Disease-Gene Associations

**Purpose**: Correlate pathway rankings with MalaCards disease-gene relevance scores.

**Method**:
1. Load trait-specific gene scores from MalaCards CSV files
2. For each pathway, calculate mean gene score
3. Rank pathways by enrichment vs. MalaCards score
4. Report Spearman correlation between rankings

**Data format** (`{trait}_genes.csv`):
```
ENSEMBL,score
ENSG00000123456,85
ENSG00000789012,42
```

**Trait restriction**: Only analyzes traits in `params.malacards_traits` (comma-separated string).

**Output files**:
```
{trait}_{tool}_malacards_rank_correlation_summary.csv
{trait}_{tool}_malacards_correlation.pdf
```

**Key metrics**:
- `rho`: Spearman correlation between pathway ranks
- `p_value`: Significance of correlation
- `num_pathways`: Number of pathways with MalaCards evidence

### GTEx Tissue Specificity

**Purpose**: Identify whether top-ranked pathways show tissue-specific expression patterns.

**Method**:
1. Load GTEx tissue-specific expression data (τ scores)
2. For each pathway, calculate mean tissue specificity across genes
3. Correlate pathway rankings with tissue specificity
4. Report correlations across all tissues

**Data format** (CSV):
```
Name,Adipose,Brain,Heart,Liver,...
HK1,0.12,0.89,0.34,0.45,...
BCL2,0.05,0.67,0.12,0.23,...
```
- **Column 1**: Gene symbol (header = "Name")
- **Columns 2+**: Tissue-specific τ scores (0-1 scale)

**No trait restriction**: Runs on all traits.

**Output files**:
```
{trait}_{tool}_tissue_correlation_summary.csv
{trait}_{tool}_best_method_by_tissue.csv
{trait}_{tool}_tissue_correlation.pdf
```

**Key metrics**:
- `rho_tissue`: Spearman correlation for each tissue
- `best_method`: BiReWire vs. KeepPathSize winner per tissue
- `significant_tissues`: Count of tissues with p < 0.05

**Common issues**:
- Gene ID mismatch: Ensure gene symbols match between GTEx and GMT files
- Missing tissues: Check that all expected tissues have data columns

### DoRothEA Transcription Factor Regulation

**Purpose**: Test whether top-ranked pathways show coordinated transcriptional regulation.

**⚠️ CRITICAL**: This module does NOT use the standard DoRothEA TF-target database. It requires **pre-computed pairwise pathway regulatory similarity scores**.

**Method**:
1. Load pairwise pathway regulatory scores (pre-computed from DoRothEA TF data)
2. For each top pathway, calculate mean regulatory similarity to other pathways
3. Correlate pathway rankings with DoRothEA scores

**Data format** (`dorothea_pairwise_scores.csv`):
```
pathway,score
KEGG_GLYCOLYSIS,0.734
REACTOME_APOPTOSIS,0.623
GO_CELL_CYCLE,0.891
```
**Note**: This is NOT `pathway1,pathway2,score` format—it's aggregated scores per pathway.

**Output files**:
```
{trait}_{tool}_dorothea_rank_correlation_summary.csv
{trait}_{tool}_dorothea_correlation.pdf
```

**Key metrics**:
- `rho_dorothea`: Spearman correlation between rankings and DoRothEA scores
- `p_value`: Correlation significance

**Common issues**:
- File not found: Verify `params.dorothea_scores_path` in config
- Format mismatch: Ensure CSV has `pathway,score` format (not TF-target format)

---

## Output Directory Structure

```
results/
├── randomized_gene_sets/                # Generated GMT files
│   ├── random_birewire/
│   │   └── GeneSet.random{1-1000}.gmt  # BiReWire permutations
│   └── random_keeppathsize/
│       └── GeneSet.random{1-1000}.gmt  # KeepPathSize permutations
│
├── magma_real/{trait}/                  # Real MAGMA results
│   ├── {trait}_case.control.genes.annot # Gene annotations
│   ├── {trait}_case.control_genebased.genes.raw  # Gene-level results
│   └── {trait}_real_set.{method}.gsa.out         # Pathway results
│
├── magma_random/{method}/{background}/{trait}/   # MAGMA null models
│   └── {trait}_set_random{1-1000}.{method}.gsa.out
│
├── prset_real/{trait}/                  # Real PRSet results
│   └── {trait}_set.{method}.summary     # Pathway results
│
├── prset_random/{method}/{background}/{trait}/   # PRSet null models
│   └── {trait}_set_random{1-1000}.{method}.summary
│
├── gsamixer/{trait}/                    # GSA-MiXeR results
│   ├── {trait}.chr{1-22}.sumstats.gz    # Chromosome-split GWAS
│   ├── {trait}_base.json                # Baseline model
│   └── {trait}_full.go_test_enrich.csv  # Pathway enrichment
│
├── gsamixer_random/{method}/{trait}/    # GSA-MiXeR null models
│   └── {trait}_random{1-1000}_full.go_test_enrich.csv
│
├── empirical_pvalues/{tool}/{trait}/    # ⭐ MAIN RESULTS
│   ├── {trait}_{tool}_birewire_empirical_pvalues.txt
│   └── {trait}_{tool}_keeppathsize_empirical_pvalues.txt
│
├── opentargets_correlation/{trait}/     # Validation results
│   ├── {trait}_{tool}_rank_correlation_summary.csv
│   └── {trait}_{tool}_detailed_opentargets_advantage.csv
│
├── tissue_correlation/{trait}/
│   ├── {trait}_{tool}_tissue_correlation_summary.csv
│   └── {trait}_{tool}_best_method_by_tissue.csv
│
├── malacards_correlation/{trait}/
│   └── {trait}_{tool}_malacards_rank_correlation_summary.csv
│
├── dorothea_correlation/{trait}/
│   └── {trait}_{tool}_dorothea_rank_correlation_summary.csv
│
└── fpr_analysis/{tool}/{method}/        # False positive rates
    └── {trait}_{tool}_{method}_fpr_results.csv
```

## Interpreting Results

### Primary Results: Empirical P-values

**File**: `results/empirical_pvalues/{tool}/{trait}/{trait}_{tool}_{method}_empirical_pvalues.txt`

**Columns**:
```
pathway_name         p_value    beta_value  empirical_pval  std_effect_size  num_genes
KEGG_GLYCOLYSIS      0.0001     0.456       0.001          2.34             52
REACTOME_APOPTOSIS   0.0234     -0.123      0.156         -0.89            147
```

### Validation Results: OpenTargets

**File**: `results/opentargets_correlation/{trait}/{trait}_{tool}_rank_correlation_summary.csv`

**Columns**:
```
method,N,rho,p_value,advantage,advantage_p
BireWire,100,0.45,0.0023,0.67,0.0012
KeepPathSize,100,0.32,0.028,0.43,0.041
```

**Interpretations**:
- `rho`: Correlation between pathway rank and OpenTargets pathway score
  - Positive → Top pathways have more known disease genes

### Validation Results: Tissue Specificity

**File**: `results/tissue_correlation/{trait}/{trait}_{tool}_tissue_correlation_summary.csv`

**Interpretation**:
- Each row = one tissue × one method

### False Positive Rate Analysis

**File**: `results/fpr_analysis/{tool}/{method}/{trait}_{tool}_{method}_fpr_results.csv`

**Columns**:
```
pathway_name         fpr_0.05  fpr_0.01  fpr_0.001  mean_rank
KEGG_GLYCOLYSIS      0.023     0.008     0.001      234
REACTOME_APOPTOSIS   0.067     0.021     0.003      512
```

**Interpretations**:
- `fpr_0.05`: Proportion of 1000 null runs where p < 0.05
  - Ideal FPR = 0.05 (calibrated)
  - FPR >> 0.05 → Method inflates false positives for this pathway
  - FPR << 0.05 → Method is conservative
- `mean_rank`: Average rank across null runs (1 = most significant)
  - Use to identify pathways consistently ranked high by chance

---

### Scheduler Configuration

**Adapting from LSF to SLURM**:
```groovy
// In nextflow.config, replace LSF blocks with:
process {
  executor = 'slurm'
  queue = 'normal'
  
  withName: run_gene_analysis {
    cpus = 1
    memory = '50 GB'
    time = '12h'
    clusterOptions = '--account=my_account'
  }
  
  withName: gsamixer_plsa_full {
    cpus = 4
    memory = '40 GB'
    time = '24h'
    clusterOptions = '--account=my_account'
  }
}
```

**Testing locally (small dataset)**:
```groovy
process {
  executor = 'local'
  maxForks = 4  // Limit parallel processes
  
  withName: run_random_sets {
    maxForks = 2  // Only run 2 random sets at a time
  }
}
```

### Debugging Nextflow Channels

**View channel contents**:
```nextflow
// Add .view() operator to inspect channel
trait_data.view { "TRAIT DATA: $it" }

random_grouped.view { trait, files, method -> 
  "GROUPED: trait=${trait}, method=${method}, files=${files.size()}"
}
```

**Check work directory**:
```bash
# Find failed processes
find work/ -name ".exitcode" -exec grep -l "1" {} \;

# Inspect stdout/stderr
cat work/ab/cd1234.../.command.out
cat work/ab/cd1234.../.command.err

# Check command that was run
cat work/ab/cd1234.../.command.sh
```

**Resume strategies**:
```bash
# Resume from last successful step (always recommended)
nextflow run main.nf -resume

# Force re-run specific processes
rm -rf work/ab/cd1234.../  # Delete work directory for that task
nextflow run main.nf -resume

# Clean all work directories (nuclear option)
rm -rf work/
nextflow run main.nf  # Start fresh
```

### Data Integrity Checks

**Verify randomization quality**:
```r
library(GSA)

# Load real and random GMT files
real <- GSA.read.gmt("data/msigdb/c2.all.gmt_filtered.txt")
rand <- GSA.read.gmt("data/randomized_gene_sets/random_birewire_10k/GeneSet.random1.gmt")

# Check pathway size preservation
real_sizes <- sapply(real$genesets, length)
rand_sizes <- sapply(rand$genesets, length)
identical(sort(real_sizes), sort(rand_sizes))  # Should be TRUE

# Check gene degree distribution (BiReWire only)
real_degrees <- table(unlist(real$genesets))
rand_degrees <- table(unlist(rand$genesets))
cor(as.numeric(real_degrees), as.numeric(rand_degrees))  # Should be ~1.0
```

**Check empirical null distributions**:
```r
# Load empirical results
results <- read.table("results/empirical_pvalues/magma/t2d/t2d_magma_birewire_empirical_pvalues.txt", 
                      header = TRUE)

# Check distribution of empirical p-values (should be ~uniform if no signal)
hist(results$empirical_pval, breaks = 50, main = "Empirical P-value Distribution")

# Check for enrichment of low p-values
sum(results$empirical_pval < 0.05) / nrow(results)  # Should be > 0.05 if signal
```

---
## Advanced Usage

### Running Subset of Traits

**Option 1: Modify JSON**
```json
// json_files/GWAS_input_subset.json
[
  {"trait": "t2d", "gwas_file": "/path/to/t2d.txt", ...},
  {"trait": "cad", "gwas_file": "/path/to/cad.txt", ...}
]
```
```bash
nextflow run main.nf --traits_config json_files/GWAS_input_subset.json -resume
```

**Option 2: Filter in main.nf**
```nextflow
// Add at line ~75 in main.nf
trait_data = Channel.fromList(phenoConfig.collect { ... })
    .filter { trait -> trait[0] in ['t2d', 'cad', 'mdd'] }  // Only these traits
```

### Adjusting Randomization Count

**For testing (fewer permutations)**:
```bash
nextflow run main.nf --num_random_sets 100 -resume
```

**⚠️ Statistical power**: 
- Minimum: 100 permutations (empirical p-value precision = 0.01)
- Recommended: 1000 permutations (precision = 0.001)
- High-stakes: 10000 permutations (precision = 0.0001, but 10× runtime)

### Custom Gene Sets

**Using your own pathways**:
```bash
# Format GMT file (tab-delimited)
echo -e "MY_PATHWAY\tDescription\tGENE1\tGENE2\tGENE3" > custom_pathways.gmt

# Update config
params.geneset_real = "/path/to/custom_pathways.gmt"

# Generate random GMTs
Rscript scripts/core/generate_birewire_gmts.R \
  --input custom_pathways.gmt \
  --output data/randomized_gene_sets/custom_birewire \
  --num_random 1000
```

### Parallel Execution Tuning

**Control process-level parallelism**:
```groovy
// In nextflow.config
executor {
  queueSize = 500        // Max jobs in queue
  submitRateLimit = '10/1min'  // Submit rate (LSF)
}

process {
  withName: run_random_sets {
    maxForks = 50  // Run 50 random sets simultaneously
  }
}
```

**Control tool-specific parallelism**:
```groovy
process {
  withName: gsamixer_plsa_full {
    cpus = 8  // Use 8 CPUs per GSA-MiXeR job
  }
}
```

### Output Cleaning

**Remove intermediate files to save space**:
```bash
# Keep only empirical results and validation outputs
find results/ -type d -name "*_random" -exec rm -rf {} +
find results/ -type d -name "*_real" -exec rm -rf {} +
find results/ -type d -name "gsamixer" -exec rm -rf {} +

# Keep work directory for resume capability
# (delete only if you're 100% done)
```

### Containerization (Singularity)

**Convert to full Singularity workflow**:
```groovy
// In modules/{tool}/{tool}.nf
process run_gene_analysis {
  container 'library://myaccount/magma:1.10'
  
  script:
  """
  magma --bfile ${params.bfile} ...
  """
}
```

**Create MAGMA container**:
```bash
singularity pull library://sylabsed/examples/magma-gwas:1.10
# OR build from scratch (see scripts/tool_specific/pascalx/to_create_sif.sh for example)
```

## FAQ

### Q: How long does the full pipeline take?

**A**: Depends on number of traits, tools, and randomization sets:
- **BiReWire GMT generation**: 12h (one-time cost)
- **MAGMA** (per trait): 12-16h (gene analysis) + 1h × 2000 random sets
- **PRSet** (per trait): 4h + 2h × 2000 random sets
- **GSA-MiXeR** (per trait): 4h + 24h + 24h × 2000 random sets (bottleneck)
- **Empirical calculation**: 1-2h per trait-tool-method
- **Validation**: 30min per trait

**Estimated total** (3 traits, all tools, 1000 permutations × 2 methods):
- With pre-generated BiReWire: ~5-7 days (parallel execution on HPC)
- Without BiReWire pre-generation: +12h upfront

### Q: Why does PRSet exclude SCZ?

**A**: Hardcoded exclusion at `main.nf:184`:
```nextflow
prset_dedup_data = trait_data
    .filter { trait_tuple -> trait_tuple[0].toUpperCase() != "SCZ" }
```
**Reason**: SCZ GWAS may have technical issues with PRSet (e.g., specific allele coding, duplicate SNPs). Modify this filter if your SCZ GWAS is clean.

### Q: Can I use summary statistics from different genome builds?

**A**: **Not recommended** without liftover:
- MAGMA reference: hg19/GRCh37
- PRSet/GSA-MiXeR: Typically hg19/GRCh37
- If GWAS is hg38: Use [UCSC liftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) first

### Q: How do I cite this pipeline?

**A**: Cite the methods papers for each tool:
- **MAGMA**: de Leeuw et al. (2015) PLOS Comput Biol
- **PRSice**: Choi et al. (2019) GigaScience
- **GSA-MiXeR**: Frei et al. (2019) Nat Commun
- **BiReWire**: Gobbi et al. (2014) Bioinformatics
- **OpenTargets**: Ochoa et al. (2021) Nucleic Acids Res

**Pipeline citation** (add your own paper once published):
```
This analysis used the Pathway Overlap Pipeline (https://github.com/accote45/pathway_overlap),
which implements empirical null models via BiReWire and KeepPathSize randomization.
```

### Q: What's the difference between `p_value` and `empirical_pval`?

**A**:
- `p_value`: Raw p-value from enrichment tool (MAGMA/PRSet/GSA-MiXeR)
  - **Problem**: Inflated by gene connectivity and pathway size biases
  - **Use**: Ranking pathways within tool
  
- `empirical_pval`: Rank-based p-value from null distribution
  - **Calculation**: `rank(real_beta among 1001 betas) / 1001`
  - **Advantage**: Controls for connectivity and size biases
  - **Use**: Significance testing (apply Bonferroni/FDR correction)

**Recommendation**: Always report `empirical_pval` in publications.

### Q: Why do BiReWire and KeepPathSize give different results?

**A**: They test different null hypotheses:
- **BiReWire null**: Pathway composition doesn't matter (given gene connectivity)
- **KeepPathSize null**: Pathway composition doesn't matter (given pathway size)

**When both agree**: Strong evidence (robust to both biases)
**When both disagree**: Investigate pathway properties (size, hub genes, clustering)

### Q: Can I add a new enrichment tool?

**A**: Yes! Follow this pattern:

1. **Create module**: `modules/mytool/mytool.nf`
```nextflow
process run_real_mytool {
  input: tuple val(trait), path(gwas_file), ...
  output: tuple val(trait), path("${trait}_mytool_results.txt")
  script: """
  mytool --input ${gwas_file} --output ${trait}_mytool_results.txt
  """
}

process run_random_mytool {
  input: tuple val(trait), val(perm), ...
  output: tuple val(trait), path("${trait}_random${perm}.txt")
  script: """
  mytool --input ${gwas_file} --gmt ${gmt_dir}/GeneSet.random${perm}.gmt
  """
}
```

2. **Add tool config**: `scripts/core/calc_empirical.r`
```r
tool_config <- list(
  mytool = list(
    pathway_col = "pathway_name",
    pval_col = "p_value",
    beta_col = "effect_size",
    ngenes_col = "num_genes",
    required_cols = c("pathway_name", "p_value", "effect_size"),
    calc_pvalue = TRUE
  )
)
```

3. **Include in main.nf**:
```nextflow
include { run_real_mytool; run_random_mytool } from './modules/mytool/mytool.nf'

// Add to workflow
mytool_results = run_real_mytool(trait_data)
mytool_random = run_random_mytool(random_inputs)
```

### Q: How do I interpret negative beta values?

**A**:
- **Positive beta**: Pathway genes show stronger trait association than genome-wide average
- **Negative beta**: Pathway genes show weaker association (depletion)
- **Magnitude**: `|beta|` indicates strength (but use `std_effect_size` for comparisons)

**Negative betas are rare but valid**—may indicate protective pathways or confounders.

---
## Citation

If you use this pipeline in your research, please cite:

### Pipeline
```
[Add your paper citation once published]

Code repository: https://github.com/accote45/pathway_overlap
```

### Methods
**MAGMA**:
> de Leeuw, C. A., Mooij, J. M., Heskes, T., & Posthuma, D. (2015).  
> MAGMA: Generalized gene-set analysis of GWAS data.  
> *PLOS Computational Biology*, 11(4), e1004219.  
> [https://doi.org/10.1371/journal.pcbi.1004219](https://doi.org/10.1371/journal.pcbi.1004219)

**PRSice/PRSet**:
> Choi, S. W., & O'Reilly, P. F. (2019).  
> PRSice-2: Polygenic Risk Score software for biobank-scale data.  
> *GigaScience*, 8(7), giz082.  
> [https://doi.org/10.1093/gigascience/giz082](https://doi.org/10.1093/gigascience/giz082)

**GSA-MiXeR**:
> Frei, O., Holland, D., Smeland, O. B., et al. (2019).  
> Bivariate causal mixture model quantifies polygenic overlap between complex traits beyond genetic correlation.  
> *Nature Communications*, 10, 2417.  
> [https://doi.org/10.1038/s41467-019-10310-0](https://doi.org/10.1038/s41467-019-10310-0)

**BiReWire**:
> Gobbi, A., Iorio, F., Dawson, K. J., Wedge, D. C., Tamborero, D., et al. (2014).  
> Fast randomization of large genomic datasets while preserving alteration counts.  
> *Bioinformatics*, 30(17), i617-i623.  
> [https://doi.org/10.1093/bioinformatics/btu474](https://doi.org/10.1093/bioinformatics/btu474)

**OpenTargets**:
> Ochoa, D., Hercules, A., Carmona, M., et al. (2021).  
> Open Targets Platform: supporting systematic drug–target identification and prioritisation.  
> *Nucleic Acids Research*, 49(D1), D1302-D1310.  
> [https://doi.org/10.1093/nar/gkaa1027](https://doi.org/10.1093/nar/gkaa1027)

**MSigDB**:
> Liberzon, A., Birger, C., Thorvaldsdóttir, H., et al. (2015).  
> The Molecular Signatures Database (MSigDB) hallmark gene set collection.  
> *Cell Systems*, 1(6), 417-425.  
> [https://doi.org/10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004)

---

## License

This pipeline is released under the **MIT License**.

```
MIT License

Copyright (c) 2024 [Your Name/Institution]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## Contact & Support

**Issues**: [GitHub Issues](https://github.com/accote45/pathway_overlap/issues)  
**Discussions**: [GitHub Discussions](https://github.com/accote45/pathway_overlap/discussions)  
**Email**: [Your email or lab email]

**Acknowledgments**: This pipeline was developed at [Your Institution] with support from [Funding sources].

---
