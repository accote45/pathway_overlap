#### filepath: docs/TOOL_DETAILS.md
# Tool-Specific Implementation Details

## MAGMA Gene-Set Analysis

### Workflow
1. SNP-to-gene mapping: 35kb upstream/downstream windows
2. Gene-level analysis: P-value aggregation via SNP-wise mean model
3. Pathway enrichment: Competitive test controlling for gene length/density

### Key Command
```bash
magma \
  --bfile ${params.bfile} \
  --pval ${gwas_file} use=SNP,P ncol=N \
  --gene-annot ${trait}.genes.annot \
  --out ${trait}_genebased
```

### Resource Requirements
- Gene annotation: 4GB RAM, 20min
- Gene analysis: 50GB RAM, 12h
- Pathway enrichment: 4GB RAM, 1h

### Common Issues
- **Missing SNPs**: Check RSID consistency between GWAS and reference
- **High memory**: Split GWAS by chromosome if needed
- **Module load failure**: Verify `magma_gwas/1.10` is available

---

## PRSet Polygenic Risk Score Analysis

### Workflow
1. SNP deduplication: Removes ALL instances of duplicate RSIDs
2. PRS calculation: Clumping with r²=0.1, 1000kb window
3. Pathway scoring: MSigDB gene sets tested for trait association
4. Validation: Uses held-out UKB test set

### ⚠️ Critical Exclusion
SCZ trait is **hardcoded to be excluded** from PRSet analysis (line 184 in `main.nf`).

### Key Parameters
```bash
PRSice_linux \
  --a1 ${effect_allele} --a2 ${other_allele} \
  --background msigdb.genes.txt:gene \
  --clump-kb 1000 --clump-p 1.0 --clump-r2 0.1 \
  --fastscore \
  --keep ukb_test_samples.txt \
  --msigdb ${pathway_gmt_file}
```

### Resource Requirements
- Deduplication: 4GB RAM, 10min
- PRS calculation: 15GB RAM, 2-4h

### Common Issues
- **Duplicate SNPs**: Pipeline removes ALL duplicates (not just highest P)
- **Phenotype misalignment**: Ensure trait names match between GWAS JSON and UKB pheno file
- **Missing genotypes**: Check `ukb18177-qc.snplist` for QC'd SNP list

---

## GSA-MiXeR Mixed-Model Enrichment

### Workflow
1. Format conversion: GWAS → GSA-MiXeR format (RSID, CHR, POS, A1, A2, Z, N, NEFF)
2. Chromosome splitting: Per-chromosome Parquet files
3. Baseline model: Universal gene set model
4. Full model: GO pathway enrichment

### Key Commands
```bash
python mixer.py plsa \
  --sumstats ${trait}.chr*.sumstats.gz \
  --go-file-base coding_genes.txt \
  --go-file-test full_gene_set.txt
```

### Resource Requirements
- Baseline model: 15GB RAM, 4h
- Full model: **40GB RAM**, 24h (critical bottleneck)
- Chromosome splitting: 1GB per chromosome

### Common Issues
- **Memory killed**: Increase to 50GB for traits with >10M SNPs
- **Missing NEFF**: Script automatically falls back to N if NEFF column absent
- **Singularity permissions**: Ensure `--home $(pwd):/home` bind mount works

### Special Pathway Handling
- `coding_genes`: Universe control (all protein-coding genes)
- `base`: Baseline model (filtered in empirical calculation)