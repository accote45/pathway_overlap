#### filepath: docs/DETAILED_SETUP.md
# Detailed Setup Guide

## Software Dependencies

### Required Tools
- **Nextflow** ≥ 22.10 (DSL2 support required)
- **Java** 11+ (for Nextflow execution)
- **R** ≥ 4.1

### R Package Installation
```r
install.packages(c("tidyverse", "data.table", "plyr", "jsonlite", 
                   "GSA", "MatchIt", "ggplot2", "patchwork", 
                   "gridExtra", "RColorBrewer", "pheatmap"))
```

### Tool-Specific Requirements

**MAGMA** (if `params.run_magma = true`):
- Version: v1.10
- Download: https://ctg.cncr.nl/software/magma
- Requires module load capability or PATH configuration

**PRSice/PRSet** (if `params.run_prset = true`):
- Version: v2.3+
- Download: https://www.prsice.info/
- Requires UKB validation cohort data

**GSA-MiXeR** (if `params.run_gsamixer = true`):
- Container: `gsa-mixer.sif`
- Requires Singularity/Apptainer 3.0+

## Data Dependencies

### Reference Files

1. **MSigDB gene sets**:
   ```bash
   wget http://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.2.Hs/c2.all.v2023.2.Hs.symbols.gmt
   ```

2. **Randomized gene sets** (BiReWire + KeepPathSize):
   ```
   data/randomized_gene_sets/
   ├── random_birewire/
   │   └── GeneSet.random{1-1000}.gmt
   └── random_keeppathsize/
       └── GeneSet.random{1-1000}.gmt
   ```

3. **MAGMA reference panel** (1000 Genomes EUR)

4. **PRSet UKB data**:
   - Training genotypes: `ukb18177_eur_autosomes.{bed,bim,fam}`
   - QC SNP list: `ukb18177-qc.snplist`
   - Phenotype file: `ukb_phenofile_forprset.txt`
   - Train/test splits

5. **GSA-MiXeR reference**:
   - Reference panel: `1000G.EUR.QC.@.bim`
   - Baseline annotations: `baseline_v1.1_genes_chr@.annot.gz`

### Validation Databases

See [VALIDATION.md](VALIDATION.md) for details on:
- OpenTargets JSON files
- GTEx tissue expression
- MalaCards disease genes
- DoRothEA regulatory scores

## Scheduler Configuration

### Default (LSF)
The pipeline is pre-configured for LSF batch systems. No changes needed.

### Adapting to SLURM
```groovy
// In nextflow.config, replace all `executor 'lsf'` blocks:
process {
  executor = 'slurm'
  queue = 'normal'
  
  withName: run_gene_analysis {
    cpus = 1
    memory = '50 GB'
    time = '12h'
    clusterOptions = '--account=my_account'
  }
}
```

### Local Testing
```groovy
process {
  executor = 'local'
  maxForks = 4
}
```