# Pathway Overlap Pipeline - AI Copilot Instructions

## Big Picture Architecture

This is a **Nextflow DSL2 pipeline** for pathway enrichment analysis across GWAS traits with empirical p-value computation via null models. The pipeline coordinates four enrichment tools (MAGMA, PRSet, PascalX, GSA-MiXeR) and validates results against biological databases.

**Core workflow pattern:**
1. **Real pathway enrichment** → Run enrichment tools on actual gene sets
2. **Null model generation** → Run same tools on randomized gene sets (BiReWire, KeepPathSize)
3. **Empirical validation** → Compare real vs null distributions to compute standardized effect sizes
4. **Biological validation** → Cross-reference top pathways with OpenTargets, MalaCards, tissue specificity, DoRothEA

## Key Integration Points

### GWAS Input Harmonization
- Traits configured via JSON: `json_files/GWAS_input_nobmi.json` (BMI excluded from default config)
- Each trait specifies column mappings for heterogeneous GWAS formats
- Channel pattern: `trait_data` broadcasts to all enrichment workflows
- **Critical:** PRSet excludes SCZ, IBD, and AD traits (hardcoded filter in the PRSet workflow block)

### Tool-Specific Data Flows
```
MAGMA:   SNP → Gene → Pathway (requires PLINK bfile reference)
PRSet:   GWAS → SNP deduplication → PRS → Pathway  
PascalX: GWAS → gene scoring (Singularity) → Pathway enrichment
GSA-MiXeR: GWAS → chr-split → baseline → full model (uses Singularity)
```

### Randomization Strategy
- **BiReWire:** Preserves degree distribution, changes gene-pathway associations
- **KeepPathSize:** Maintains pathway sizes, randomly samples genes
- Both methods generate `params.num_random_sets` permutations per trait-method combination
- GMT files stored in `params.gmt_dirs[method]` directories (output under `${params.outdir}/randomized_gene_sets/`)

**GMT Generation and File Organization:**
- Controlled by `params.generate_random_gmts` (default: `true`) — the pipeline generates GMTs at runtime
- When `false`, uses existing files from `params.gmt_dirs`
- A `gmt_ready_signal` channel gates all random enrichment processes, ensuring GMTs complete first
```
${params.outdir}/randomized_gene_sets/
├── random_birewire/
│   ├── GeneSet.random1.gmt
│   ├── GeneSet.random2.gmt
│   └── ... (up to GeneSet.random{num_random_sets}.gmt)
└── random_keeppathsize/
    ├── GeneSet.random1.gmt
    └── ... (same pattern)
```

**Randomization Algorithm Details:**
- BiReWire: Uses degree-preserving network randomization - genes maintain same number of pathway memberships but in different pathways
- KeepPathSize: Pathways maintain exact size, genes randomly reassigned - pathway functional coherence is disrupted
- Both generate identical file naming: `GeneSet.random{1-N}.gmt`

## Critical Development Patterns

### Channel Architecture
The pipeline uses **grouped channels** to coordinate real/random results. This is the most complex aspect:

```nextflow
// 1. Gate random enrichment on GMT generation completing
random_sets_inputs = gene_analysis_results.gene_results
    .combine(gmt_ready_signal)          // Wait for GMTs to be written
    .map { trait, gene_result, ready -> tuple(trait, gene_result) }
    .combine(Channel.fromList(params.randomization_methods))
    .combine(Channel.from(1..params.num_random_sets))
    .map { trait, gene_result, rand_method, perm -> 
        tuple(trait, gene_result, rand_method, perm)
    }

// 2. Group random results by trait+method with deterministic size
random_results_grouped = random_sets_results
    .groupTuple(by: [0, 2], size: params.num_random_sets)

// 3. Broadcast real results to all randomization methods, then join
magma_real_for_empirical = real_geneset_results
    .map { trait, real_result_file, dummy -> tuple(trait, real_result_file) }
    .combine(Channel.fromList(params.randomization_methods))

magma_for_empirical = magma_real_for_empirical
    .combine(random_results_grouped, by: [0, 2])
    .map { trait, rand_method, real_result_file, random_files_list ->
        def random_dir = "${params.outdir}/magma_random/${rand_method}/${params.background}/${trait}"
        tuple(trait, "magma_${rand_method}", real_result_file, random_dir)
    }
```

**Critical Channel Patterns:**
- `by: [0, 2]` groups by trait (index 0) and randomization method (index 2)
- `size: params.num_random_sets` makes `groupTuple` deterministic — emits only when all permutations arrive
- Real results get a dummy `rand_method = "deduplicate"` for channel tuple consistency, then broadcast
- Empirical tool identifiers are `{tool}_{rand_method}` (e.g., `magma_birewire`, `prset_keeppathsize`)
- `gmt_ready_signal` is created by combining birewire and keeppathsize output dirs — gates all random enrichment

**Validation Channel Helper (`group_by_trait_tool`):**
Empirical results are post-processed before validation. The inline closure:
1. Strips `rand_method` suffix to get `base_tool`
2. Groups by `[trait, base_tool]`
3. Produces `[trait, base_tool, birewire_file, keeppathsize_file]` tuples for each validation module

### Module Structure
- `modules/*/` contain process definitions
- `scripts/` contain R analysis code called by processes
- Each tool has separate real/random modules due to different parameters

### R Script Conventions
- Tool-agnostic empirical calculation: `scripts/core/calc_empirical.r`
- Tool config mapping handles column differences:
  ```r
  tool_config <- list(
    magma = list(pathway_col = "FULL_NAME", pval_col = "P", ...),
    prset = list(pathway_col = "Set", pval_col = "P", ...),
    pascalx = list(pathway_col = "pathway", pval_col = "p_value", ...),
    gsamixer = list(pathway_col = "GO", pval_col = NULL, ...)
  )
  ```

**Script Parameter Patterns:**
```r
# Standard calling convention for empirical scripts
args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]        # e.g., "mdd", "cad"
full_tool <- args[2]    # e.g., "magma_birewire", "prset_keeppathsize", "gsamixer"
real_file <- args[3]    # Path to real enrichment results
random_dir <- args[4]   # Directory containing random files
gmt_path <- args[5]     # Optional: GMT file for validation
```

**File Format Expectations:**
- MAGMA outputs: `*.gsa.out` (tab-delimited, FULL_NAME column)
- PRSet outputs: `*.summary` (includes Competitive.P, PRS.R2 columns)
- PascalX outputs: pathway enrichment TSV (pathway, p_value columns)
- GSA-MiXeR: `*.go_test_enrich.csv` (JSON also created but CSV used for analysis)

## Environment & Dependencies

### Scheduler Configuration
- **Default:** LSF with project codes (`clusterOptions = '-P acc_paul_oreilly'`)
- **Memory patterns:** Base processes (4-15GB), GSA-MiXeR full models (40GB, 24h)
- **To adapt:** Remove LSF blocks in `nextflow.config`, set global executor

### External Tool Requirements
- **MAGMA:** Requires `module load magma_gwas/1.10` and reference PLINK files
- **PascalX:** Uses Singularity container (`params.pascalx_sif`); requires `module load apptainer/1.3.6`; reference panel and genome annotations mounted read-only via Singularity bind
- **GSA-MiXeR:** Uses Singularity container, requires baseline annotations
- **R packages:** tidyverse, data.table, GSA, MatchIt (see README for full list)

## Debugging Workflows

### Common Failure Points
1. **Missing random files:** Check `params.gmt_dirs` paths and file counts
2. **Column mapping errors:** Verify GWAS JSON column names match file headers  
3. **Memory failures:** GSA-MiXeR full models need 40GB+ for large trait sets
4. **Channel grouping issues:** Random results must complete before empirical calculation

**Specific Debugging Scenarios:**
```bash
# Verify all random files exist before empirical calculation
for trait in mdd cad t2d; do
  count=$(find results/magma_random/birewire/msigdbgenes/$trait -name "*.gsa.out" | wc -l)
  echo "$trait has $count random files (should be 1000)"
done

# Check GWAS column mapping issues
head -1 /path/to/gwas.txt | tr '\t' '\n' | nl  # Show column numbers

# Monitor GSA-MiXeR memory usage (common bottleneck)
watch 'ps aux | grep "mixer.py" | grep -v grep'
```

**Channel Debugging Patterns:**
- Use `.view()` operator to inspect channel contents: `channel.view { "DEBUG: $it" }`
- Check grouped channels: Random results should show `[trait, [file1, file2, ...], method]`
- Empirical inputs should show `[trait, tool, real_file, random_dir]` pattern

### Output Structure
```
results/
├── magma_real/{trait}/           # Real MAGMA enrichment results
├── magma_random/{method}/{background}/{trait}/  # MAGMA null distributions
├── pascalx_genes/{trait}/        # PascalX gene scores
├── pascalx_random/{method}/{background}/{trait}/ # PascalX null distributions
├── prset_random/{method}/{background}/{trait}/   # PRSet null distributions
├── gsamixer_random/{method}/{trait}/             # GSA-MiXeR null distributions
├── {trait}_{tool}_empirical_pvalues.txt         # Final empirical results
└── opentargets/, tissue/, malacards/, dorothea/ # Validation outputs
```

**Validation Workflow Details:**
1. **OpenTargets Analysis** (`modules/opentargets/`):
   - Size-matched pathway comparison against gene-disease associations
   - Requires `params.opentargets_supported_traits` whitelist
   - Uses JSON files from associationByDatatypeDirect
   - Outputs: `*_detailed_advantage.csv`, combined plots

2. **Tissue Specificity** (`modules/tissuespecificity/`):
   - Correlates pathway rankings with tissue expression specificity
   - Independent of OpenTargets, runs on all traits
   - Uses precomputed tissue expression data
   - Outputs: correlation summaries, multi-page PDFs

3. **MalaCards Integration** (`modules/malacards/`):
   - Disease-gene association validation
   - Trait-filtered via `params.malacards_traits` (e.g., `'bmi,cad,t2d,mdd,ad,scz,ibd,breast'`)

4. **DoRothEA Correlation** (`modules/dorothea/`):
   - Transcription factor target gene enrichment correlation
   - Runs on all traits; uses `params.dorothea_scores_path`
   - Enabled/disabled by `params.run_dorothea_correlation`

**Critical Validation Dependencies:**
- OpenTargets: `params.opentargets_json_dir` must contain association JSONs
- Tissue: `params.tissue_expression_data` must have gene symbols in 'Name' column
- DoRothEA: `params.dorothea_scores_path` must point to pairwise scores CSV
- All validation requires empirical results to complete first

### Key Debugging Commands
```bash
# Check random file generation
find results/magma_random/birewire -name "*.gsa.out" | wc -l  # Should be N * num_traits

# Verify empirical results structure
head results/*/empirical_pvalues.txt  # Check for required columns

# Monitor GSA-MiXeR singularity usage
ps aux | grep singularity

# Validate GMT file structure
head -5 results/randomized_gene_sets/random_birewire/GeneSet.random1.gmt

# Check trait filtering for OpenTargets
grep -i "opentargets_supported_traits" nextflow.config

# Verify GWAS column mappings match actual files
trait="mdd"; head -1 $(jq -r '.[] | select(.trait=="'$trait'") | .gwas_file' json_files/GWAS_input_nobmi.json)
```

**Common Error Patterns:**
- "No such file" errors: Usually GMT paths in `params.gmt_dirs` incorrect, or `generate_random_gmts = false` used before GMTs exist
- "Column not found" in R: GWAS JSON column names don't match file headers
- Memory killed: GSA-MiXeR processes need 40GB, check resource allocation
- Empty empirical results: Random file count != `num_random_sets` per trait-method
- PascalX Singularity failures: Check `module load apptainer/1.3.6` is available and bind paths are correct

## Project-Specific Conventions

- **Trait filtering:** OpenTargets analyses restricted to whitelist in `params.opentargets_supported_traits`
- **Background gene sets:** Use `params.background = 'msigdbgenes'` for consistency
- **File publishing:** Real results published to tool-specific dirs, random to method-specific nested structure
- **Channel naming:** `*_by_trait_method` indicates grouped empirical results ready for downstream analysis

This pipeline's complexity stems from coordinating four independent tools with shared null models—understanding the channel grouping patterns and tool-specific requirements is essential for modifications.