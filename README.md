# Pathway Overlap Pipeline

**Nextflow pipeline for pathway enrichment analysis with empirical validation**

## What This Does

Tests whether gene sets (pathways) are enriched in GWAS traits, using **null models** to calculate adjusted enrichment statistics accounting for pathway overlap and/or pathway size. Runs three complementary enrichment tools (MAGMA, PRSet, GSA-MiXeR) and validates results against biological databases.

---

## Quick Start

### 1. Install Dependencies
```bash
# Nextflow
curl -s https://get.nextflow.io | bash

# R packages
Rscript -e 'install.packages(c("tidyverse", "data.table", "GSA"))'
```

**Tool-specific requirements**: See [Detailed Setup Guide](docs/DETAILED_SETUP.md)

### 2. Clone & Configure
```bash
git clone https://github.com/accote45/pathway_overlap.git
cd pathway_overlap

# Edit GWAS trait configurations
nano json_files/GWAS_input.json

# Update paths in nextflow.config
nano nextflow.config
```

**Configuration help**: See [GWAS Input Guide](docs/DETAILED_SETUP.md#gwas-configuration)

### 3. Pre-generate Randomized Gene Sets
**⚠️ Required before running pipeline** (12h runtime)
```bash
Rscript scripts/core/generate_birewire_gmts.R \
  --input data/msigdb/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt \
  --output data/randomized_gene_sets/random_birewire_10k \
  --num_random 1000
```

### 4. Run Pipeline
```bash
# Full analysis
nextflow run main.nf -resume

# MAGMA only (for testing)
nextflow run main.nf --run_magma true --run_prset false --run_gsamixer false -resume

# Custom parameters
nextflow run main.nf \
  --num_random_sets 500 \
  --outdir /scratch/my_results \
  -resume
```

---

## Pipeline Overview

```
┌─────────────────────────────────────────────────┐
│ 1. Real Pathway Enrichment (MAGMA/PRSet/MiXeR)│
└─────────────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────────────┐
│ 2. Null Model Enrichment (1000 random sets)    │
│    - BiReWire: Preserves gene degree distribution │
│    - KeepPathSize: Preserves pathway sizes      │
└─────────────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────────────┐
│ 3. Empirical Statistics                         │
│    - empirical_pval = rank(real) / 1001        │
│    - std_effect = (real - μ_null) / σ_null     │
└─────────────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────────────┐
│ 4. Biological Validation (OpenTargets, GTEx)   │
└─────────────────────────────────────────────────┘
```

---

## Key Results

### Main Output
```
results/empirical_stats/{tool}/{trait}/
  └── {trait}_{tool}_{method}_empirical_pvalues.txt
```

**Columns**:
- `empirical_pval`: Rank-based p-value from null distribution (⭐ use this, not `p_value`)
- `std_effect_size`: Standardized effect (like Z-score)
- `p_value`: Original tool p-value (for reference only)

### Validation Outputs
```
results/
├── opentargets_correlation/{trait}/
├── tissue_correlation/{trait}/
├── malacards_correlation/{trait}/
└── dorothea_correlation/{trait}/
```

See [Output Guide](docs/DETAILED_SETUP.md#output-structure) for details.

---

## Understanding Results

### Interpreting Empirical P-values
```r
# Load results
results <- read.table("results/empirical_stats/magma/t2d/t2d_magma_birewire_empirical_pvalues.txt", 
                      header=TRUE)

# Top pathways (empirical p < 0.05)
top_pathways <- results[results$empirical_pval < 0.05 & results$std_effect_size > 0, ]
top_pathways <- top_pathways[order(top_pathways$empirical_pval), ]

# Check validation
validation <- read.table("results/opentargets_correlation/t2d/t2d_magma_rank_correlation_summary.csv",
                         header=TRUE, sep=",")
```

### Comparing Randomization Methods
- **BiReWire advantage**: Top-ranked pathways show better biological validation
- **KeepPathSize advantage**: More specific pathway size effects
- **Use BiReWire by default** unless you specifically need size-matched controls

---

## Documentation

- **[Detailed Setup Guide](docs/DETAILED_SETUP.md)**: Software/data dependencies, scheduler configuration
- **[Tool Details](docs/TOOL_DETAILS.md)**: MAGMA, PRSet, GSA-MiXeR implementation specifics
- **[Validation Guide](docs/VALIDATION.md)**: OpenTargets, GTEx, MalaCards, DoRothEA
- **[Troubleshooting](docs/TROUBLESHOOTING.md)**: Common errors and debugging strategies

---

## FAQ

**Q: How long does this take?**  
A: 5-7 days for 3 traits with all tools (1000 permutations × 2 methods). MAGMA-only: ~2 days.

**Q: Why does PRSet exclude SCZ?**  
A: Hardcoded exclusion in `main.nf` line 184 due to UKB data restrictions. Remove if you have alternative data.

**Q: What's the difference between `p_value` and `empirical_pval`?**  
A: 
- `p_value`: Original tool output (asymptotic, can be anti-conservative)
- `empirical_pval`: Rank-based value from null distribution (**use this**)

**Q: Can I use my own gene sets?**  
A: Yes! Format as GMT file and update `params.geneset_real`. See [Custom Gene Sets](docs/DETAILED_SETUP.md#custom-gene-sets).

---

## Citation

### Pipeline
```
[Add your paper citation once published]
Repository: https://github.com/accote45/pathway_overlap
```

### Methods
- **MAGMA**: de Leeuw et al. (2015) *PLOS Comput Biol*
- **PRSice**: Choi & O'Reilly (2019) *GigaScience*
- **GSA-MiXeR**: Frei et al. (2019) *Nat Commun*
- **BiReWire**: Gobbi et al. (2014) *Bioinformatics*

---

## License

MIT License - see [LICENSE](LICENSE) file

---

## Support

**Issues**: [GitHub Issues](https://github.com/accote45/pathway_overlap/issues)  
**Discussions**: [GitHub Discussions](https://github.com/accote45/pathway_overlap/discussions)
