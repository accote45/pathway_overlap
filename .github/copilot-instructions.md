# Pathway Overlap Pipeline - AI Copilot Instructions

## Big Picture Architecture

This is a **Nextflow DSL2 pipeline** for pathway enrichment analysis across GWAS traits with empirical p-value computation via null models. The pipeline coordinates three enrichment tools (MAGMA, PRSet, GSA-MiXeR) and validates results against biological databases.

**Core workflow pattern:**
1. **Real pathway enrichment** → Run enrichment tools on actual gene sets
2. **Null model generation** → Run same tools on randomized gene sets (BiReWire, KeepPathSize)
3. **Empirical validation** → Compare real vs null distributions to compute standardized effect sizes
4. **Biological validation** → Cross-reference top pathways with OpenTargets, MalaCards, tissue specificity

## Key Integration Points

### GWAS Input Harmonization
- Traits configured via JSON: `json_files/GWAS_input.json`
- Each trait specifies column mappings for heterogeneous GWAS formats
- Channel pattern: `trait_data` broadcasts to all enrichment workflows
- **Critical:** PRSet excludes SCZ trait (hardcoded filter in `main.nf:184`)

### Tool-Specific Data Flows
```
MAGMA: SNP → Gene → Pathway (requires PLINK bfile reference)
PRSet: GWAS → SNP deduplication → PRS → Pathway  
GSA-MiXeR: GWAS → chr-split → baseline → full model (uses Singularity)
```

### Randomization Strategy
- **BiReWire:** Preserves degree distribution, changes gene-pathway associations
- **KeepPathSize:** Maintains pathway sizes, randomly samples genes
- Both methods generate 1000 permutations per trait-method combination
- GMT files stored in `params.gmt_dirs[method]` directories

**GMT File Organization:**
```
data/randomized_gene_sets/
├── random_birewire_10k/
│   ├── GeneSet.random1.gmt
│   ├── GeneSet.random2.gmt
│   └── ... (up to GeneSet.random1000.gmt)
└── random_keeppathsize_10k/
    ├── GeneSet.random1.gmt
    └── ... (same pattern)
```

**Randomization Algorithm Details:**
- BiReWire: Uses degree-preserving network randomization - genes maintain same number of pathway memberships but in different pathways
- KeepPathSize: Pathways maintain exact size, genes randomly reassigned - pathway functional coherence is disrupted
- Both generate identical file naming: `GeneSet.random{1-1000}.gmt`

## Critical Development Patterns

### Channel Architecture
The pipeline uses **grouped channels** to coordinate real/random results. This is the most complex aspect:

```nextflow
// 1. Generate 1000 random results per trait+method combination
random_sets_inputs = gene_results_with_methods
    .combine(Channel.from(1..params.num_random_sets))
    .map { trait, gene_result, rand_method, perm -> 
        [trait, gene_result, rand_method, perm]
    }

// 2. Group random results by trait+method (collects all 1000 files)
random_results_grouped = random_sets_results.groupTuple(by: [0, 2])

// 3. Join real results with ALL random results for empirical calculation
empirical_inputs = real_results.combine(random_results_grouped, by: [0, 2])
    .map { trait, result_file, rand_method, random_files ->
        def random_dir = "${params.outdir}/magma_random/${rand_method}/${params.background}/${trait}"
        tuple(trait, "magma", result_file, random_dir)
    }
```

**Critical Channel Patterns:**
- `by: [0, 2]` groups by trait (index 0) and randomization method (index 2)
- Random results must complete BEFORE empirical calculation starts
- Each trait-method pair needs exactly 1000 random files for valid statistics

### Module Structure
- `modules/*/` contain process definitions
- `scripts/` contain R analysis code called by processes
- Each tool has separate real/random modules due to different parameters

### R Script Conventions
- Tool-agnostic empirical calculation: `scripts/calc_empirical.r`
- Tool config mapping handles column differences:
  ```r
  tool_config <- list(
    magma = list(pathway_col = "FULL_NAME", pval_col = "P", ...),
    prset = list(pathway_col = "Set", pval_col = "P", ...),
    gsamixer = list(pathway_col = "GO", pval_col = NULL, ...)
  )
  ```

**Script Parameter Patterns:**
```r
# Standard calling convention for empirical scripts
args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]        # e.g., "mdd", "cad"
full_tool <- args[2]    # e.g., "magma", "prset", "gsamixer"  
real_file <- args[3]    # Path to real enrichment results
random_dir <- args[4]   # Directory containing 1000 random files
gmt_path <- args[5]     # Optional: GMT file for validation
```

**File Format Expectations:**
- MAGMA outputs: `*.gsa.out` (tab-delimited, FULL_NAME column)
- PRSet outputs: `*.summary` (includes Competitive.P, PRS.R2 columns)
- GSA-MiXeR: `*.go_test_enrich.csv` (JSON also created but CSV used for analysis)

## Environment & Dependencies

### Scheduler Configuration
- **Default:** LSF with project codes (`clusterOptions = '-P acc_paul_oreilly'`)
- **Memory patterns:** Base processes (4-15GB), GSA-MiXeR full models (40GB, 24h)
- **To adapt:** Remove LSF blocks in `nextflow.config`, set global executor

### External Tool Requirements
- **MAGMA:** Requires `module load magma_gwas/1.10` and reference PLINK files
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
├── magma_real/{trait}/           # Real enrichment results
├── magma_random/{method}/{background}/{trait}/  # Null distributions  
├── {trait}_{tool}_empirical_pvalues.txt        # Final empirical results
└── opentargets/, tissue/, malacards/           # Validation outputs
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
   - Runs correlation analysis
   - No trait restrictions (unlike OpenTargets)

**Critical Validation Dependencies:**
- OpenTargets: `params.opentargets_json_dir` must contain association JSONs
- Tissue: `params.tissue_expression_data` must have gene symbols in 'Name' column
- All validation requires empirical results to complete first

### Key Debugging Commands
```bash
# Check random file generation
find results/magma_random/birewire -name "*.gsa.out" | wc -l  # Should be 1000 * num_traits

# Verify empirical results structure
head results/*/empirical_pvalues.txt  # Check for required columns

# Monitor GSA-MiXeR singularity usage
ps aux | grep singularity

# Validate GMT file structure
head -5 data/randomized_gene_sets/random_birewire_10k/GeneSet.random1.gmt

# Check trait filtering for OpenTargets
grep -i "opentargets_supported_traits" nextflow.config

# Verify GWAS column mappings match actual files
trait="mdd"; head -1 $(jq -r '.[] | select(.trait=="'$trait'") | .gwas_file' json_files/GWAS_input.json)
```

**Common Error Patterns:**
- "No such file" errors: Usually GMT paths in `params.gmt_dirs` incorrect
- "Column not found" in R: GWAS JSON column names don't match file headers
- Memory killed: GSA-MiXeR processes need 40GB, check resource allocation
- Empty empirical results: Random file count != 1000 per trait-method

## Project-Specific Conventions

- **Trait filtering:** OpenTargets analyses restricted to whitelist in `params.opentargets_supported_traits`
- **Background gene sets:** Use `params.background = 'msigdbgenes'` for consistency
- **File publishing:** Real results published to tool-specific dirs, random to method-specific nested structure
- **Channel naming:** `*_by_trait_method` indicates grouped empirical results ready for downstream analysis

This pipeline's complexity stems from coordinating three independent tools with shared null models—understanding the channel grouping patterns and tool-specific requirements is essential for modifications.