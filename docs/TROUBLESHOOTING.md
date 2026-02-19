#### filepath: docs/TROUBLESHOOTING.md
# Troubleshooting Guide

## Common Errors

### 1. Missing Random Files
**Error**: `No random result files found`

**Solution**:
```bash
# Verify GMT files exist
find data/randomized_gene_sets -name "*.gmt" | wc -l
# Should show 2000 (1000 × 2 methods)

# Check file naming
ls data/randomized_gene_sets/random_birewire/ | head
# Should show: GeneSet.random1.gmt, GeneSet.random2.gmt, etc.
```

### 2. GWAS Column Mapping Errors
**Error**: `Column not found: SNP`

**Solution**:
```bash
# Show actual column names with numbers
head -1 /path/to/gwas.txt | tr '\t' '\n' | nl

# Update JSON to match exact column names (case-sensitive)
```

### 3. GSA-MiXeR Memory Killed
**Error**: Process killed by scheduler

**Solution**:
```groovy
// In nextflow.config
process {
  withName: gsamixer_plsa_full {
    memory = '50 GB'  // Increase from 40GB
    time = '36h'      // Increase from 24h
  }
}
```

### 4. Channel Grouping Issues
**Error**: Empirical calculation runs before random results complete

**Debug**:
```nextflow
// Add .view() to inspect channel
random_results_grouped.view { "GROUPED: ${it[0]}, files=${it[1].size()}" }
// Should show exactly 1000 files per trait-method
```

## Debugging Workflows

### Inspect Work Directory
```bash
# Find failed processes
find work/ -name ".exitcode" -exec grep -l "1" {} \;

# Check stdout/stderr
cat work/ab/cd1234.../.command.out
cat work/ab/cd1234.../.command.err

# See executed command
cat work/ab/cd1234.../.command.sh
```

### Resume Strategies
```bash
# Resume from last successful step (recommended)
nextflow run main.nf -resume

# Force re-run specific processes
rm -rf work/ab/cd1234.../
nextflow run main.nf -resume

# Nuclear option: start fresh
rm -rf work/
nextflow run main.nf
```

### Verify Randomization Quality
```r
library(GSA)

real <- GSA.read.gmt("data/msigdb/c2.all.gmt_filtered.txt")
rand <- GSA.read.gmt("data/randomized_gene_sets/random_birewire_10k/GeneSet.random1.gmt")

# Check pathway size preservation
real_sizes <- sapply(real$genesets, length)
rand_sizes <- sapply(rand$genesets, length)
identical(sort(real_sizes), sort(rand_sizes))  # Should be TRUE
```

## Performance Issues

### Pipeline Too Slow
```bash
# Run only MAGMA (fastest tool)
nextflow run main.nf --run_prset false --run_gsamixer false -resume

# Reduce randomization count for testing
nextflow run main.nf --num_random_sets 100 -resume

# Use more parallel jobs
# In nextflow.config:
executor.queueSize = 1000
```

### Out of Disk Space
```bash
# Remove intermediate files
find results/ -type d -name "*_random" -exec rm -rf {} +
find results/ -type d -name "gsamixer" -exec rm -rf {} +

# Clean work directory (only after pipeline completes)
rm -rf work/
```