#### filepath: docs/VALIDATION.md
# Biological Validation Guide

## OpenTargets Gene-Disease Evidence

### Purpose
Test whether top-ranked pathways are enriched for genes with known disease associations.

### Required Data
```bash
# Download associationByDatatypeDirect JSON files
opentargets_data/
├── EFO_0000311.json  # Schizophrenia
├── EFO_0001645.json  # Coronary artery disease
├── MONDO_0005148.json  # Type 2 diabetes
└── ... (7 total required)
```

### Trait Mapping
```
t2d     → MONDO_0005148
cad     → EFO_0001645
scz     → MONDO_0005090
mdd     → MONDO_0002050
ad      → MONDO_0004975
ibd     → EFO_0000555
breast  → MONDO_0007254
```

### Output Files
- `{trait}_{tool}_rank_correlation_summary.csv`
- `{trait}_{tool}_detailed_opentargets_advantage.csv`
- `{trait}_{tool}_opentargets_correlation.pdf`

### Key Metrics
- `rho_ot`: Spearman correlation between pathway ranking and OpenTargets scores
- `p_value`: Correlation significance
- `advantage_birewire`: BiReWire rank improvement over KeepPathSize

---

## GTEx Tissue Specificity

### Purpose
Test if pathway rankings correlate with tissue-specific expression patterns.

### Required Data
- CSV with 'Name' column (gene symbols) + tissue expression columns
- Download: https://gtexportal.org/home/datasets

### Output Files
- `{trait}_{tool}_tissue_correlation_summary.csv`
- `{trait}_{tool}_best_method_by_tissue.csv`
- `{trait}_{tool}_tissue_correlation.pdf`

### Key Metrics
- `rho_tissue`: Spearman correlation for each tissue
- `best_method`: BiReWire vs KeepPathSize winner per tissue
- `significant_tissues`: Count of tissues with p < 0.05

---

## MalaCards Disease-Gene Associations

### Purpose
Validate pathway rankings against disease-specific gene databases.

### Required Data
Path pattern: `${params.malacards_path}/{trait}_genes.csv`

Required traits: `bmi, cad, t2d, mdd, ad, scz, ibd, breast`

CSV format:
```
ENSEMBL,score
ENSG00000012048,8.5
ENSG00000139618,7.2
```

### Output Files
- `{trait}_{tool}_malacards_rank_correlation_summary.csv`

---

## DoRothEA Transcription Factor Regulation

### Purpose
Test pathway similarity based on shared regulatory mechanisms.

### ⚠️ Critical Requirement
This is **NOT** the standard DoRothEA TF-target database.

Requires **pre-computed pairwise pathway scores**: `dorothea_pairwise_scores.csv`

Format:
```
pathway1,pathway2,score
KEGG_GLYCOLYSIS,REACTOME_APOPTOSIS,0.456
```

### Output Files
- `{trait}_{tool}_dorothea_rank_correlation_summary.csv`
- `{trait}_{tool}_dorothea_correlation.pdf`

### Key Metrics
- `rho_dorothea`: Spearman correlation between rankings and DoRothEA scores
- `p_value`: Correlation significance