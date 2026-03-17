# references:
# https://github.com/bulik/ldsc/wiki/Partitioned-Heritability
# https://github.com/bulik/ldsc
# https://hgen471.hakyimlab.org/post/2022/02/22/partition-heritability/#partitioned-heritability

# model:
## the χ2 association statistic for a given SNP includes the effects of all SNPs tagged by this SNP. Thus, for a polygenic trait, SNPs with a high LD score will have higher χ2 statistics on average than SNPs with a low LD score16. This phenomenon might be driven either by the higher likelihood of these SNPs tagging an individual large effect or their ability to tag multiple weak effects.
## If we partition SNPs into functional categories with different contributions to heritability, then LD to a category that is enriched for heritability will increase the χ2 statistic of a SNP more than LD to a category that does not contribute to heritability.
## our method determines that a category of SNPs is enriched for heritability if SNPs with high LD to that category have higher χ2 statistics than SNPs with low LD to that category.

#### checks:
# When you implement the multipathway vs other genes annotations, the LDSC tutorial assumes you're adding your annotations to the baseline model, not replacing it. That step is critical — otherwise the enrichment can just reflect gene density or LD structure.

# Check whether mean chi2 is above ~1.02. If no, this means there is very little polygenic signal for LDSC to work with.

# A. By default, ldsc reports h2 on the observed scale. Most publications report h2 on the liability scale, since this allows for comparison across studies of the same disease with different proportions of cases and controls. You can use the --samp-prev and --pop-prev flags to convert to liability scale h2 and genetic covariance
    # observed scale = proportion of phenotypic variance in the trait due to genetic variation


#### all in /sc/arion/projects/psychgen/cotea02_prset/pathway_overlap/exploratory/sldsc

BIMDIR=/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/gsa_mixer/1000G_EUR_Phase3_plink
# Path prefix to official pre-computed baselineLD v2.2 LD scores: /sc/arion/projects/psychgen/projects/prs/sample_overlap/analysis/prepare/data/reference
# tar -xvf 1000G_Phase3_baselineLD_v2.2_ldscores.tgz -C data/ref/
BASELINE_LDSCORES=/sc/arion/projects/psychgen/projects/prs/sample_overlap/analysis/prepare/data/reference/baseline-

mkdir -p annot ldscores results data/gene_lists data/gene_coords

# --------------------------------------------------------------------------
# munge GWAS sumstats for LDSC
# --------------------------------------------------------------------------
ml ldsc/1.0.1

snplist=/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/results/ldsc/resource/eur_w_ld_chr/w_hm3.snplist

# Extract SNP-ID-only list for --print-snp (w_hm3.snplist is multi-column SNP/A1/A2;
# --print-snps requires a single SNP column with header)
awk 'NR==1{for(i=1;i<=NF;i++) if($i=="SNP") col=i; print "SNP"; next} {print $col}' \
    ${snplist} > data/hm3_snpids_only.txt
hm3_snpids=data/hm3_snpids_only.txt

munge_sumstats.py \
    --sumstats /sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/gwas/ibd_gwas.txt \
    --merge-alleles $snplist \
    --a1 Allele1 \
    --a2 Allele2 \
    --N-col N \
    --p P.value \
    --snp Rsid \
    --signed-sumstats Effect,0 \
    --chunksize 500000 \
    --out ibd

# --------------------------------------------------------------------------
# get list of multipathway genes (top 10%) and other genes
# --------------------------------------------------------------------------
ml R

Rscript - <<'EOF'
library(data.table)
gmt_file <- "/sc/arion/projects/psychgen/cotea02_prset/pathway_overlap/data/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt"

# Read GMT: col1=pathway, col2=placeholder, col3+=genes
lines <- readLines(gmt_file)
genes <- unlist(lapply(lines, function(line) {
  fields <- strsplit(line, "\t")[[1]]
  fields[-(1:2)]  # drop pathway name and placeholder
}))

# Count pathway memberships per gene
counts <- sort(table(genes))

# Top 10% threshold (90th percentile)
threshold <- quantile(counts, 0.90)
multipathway_genes <- names(counts[counts >= threshold])

cat(sprintf("Threshold: %d pathways\n", threshold))
cat(sprintf("%d multipathway genes\n", length(multipathway_genes)))

writeLines(multipathway_genes, "data/gene_lists/multipathway_genes.txt")
other_genes <- names(counts[counts < threshold])
writeLines(other_genes, "data/gene_lists/other_genes.txt")
cat(sprintf("%d other genes\n", length(other_genes)))
EOF

# --------------------------------------------------------------------------
# reformat gene coordinates file
# --------------------------------------------------------------------------
Rscript - <<'EOF'
library(data.table)

coords <- fread(
  "data/msigdbgenes.regions",
  col.names = c("ensembl_id", "CHR", "START", "END", "strand", "GENE"),
  colClasses = c("character", "character", "integer", "integer", "character", "character")
)

# Keep only autosomes
coords <- coords[CHR %in% c(as.character(1:22))]

# Write in LDSC make_annot.py expected format
fwrite(
  coords[, .(GENE = ensembl_id, CHR, START, END)],
  file = "data/gene_coords/gene_coords.txt",
  sep = "\t"
)
EOF

# --------------------------------------------------------------------------
# Generate SNP annotations using make_annot.py
# --------------------------------------------------------------------------
ml ldsc/1.0.1 bedtools

for chr in {1..22}; do
    make_annot.py \
        --gene-set-file data/gene_lists/multipathway_genes.txt \
        --gene-coord-file data/gene_coords/gene_coords.txt \
        --windowsize 35000 \
        --bimfile ${BIMDIR}/1000G.EUR.QC.${chr}.bim \
        --annot-file annot/multipathway.${chr}.annot.gz

    make_annot.py \
        --gene-set-file data/gene_lists/other_genes.txt \
        --gene-coord-file data/gene_coords/gene_coords.txt \
        --windowsize 35000 \
        --bimfile ${BIMDIR}/1000G.EUR.QC.${chr}.bim \
        --annot-file annot/othergenes.${chr}.annot.gz
done

# --------------------------------------------------------------------------
# Calculate allele frequency from 1000G EUR reference panel
# Mirrors genotype_filtering process in prset_ldsc.nf:
#   - exclude duplicate SNPs (same awk deduplication logic)
#   - apply minor allele count filter (--mac 5) and genotype missingness (--geno 0.01)
#   - compute allele frequencies with --freq
# --------------------------------------------------------------------------
ml plink

mkdir -p data/freq

for chr in {1..22}; do
    bim=${BIMDIR}/1000G.EUR.QC.${chr}.bim

    # identify duplicate SNP IDs
    awk '$2 in a {print $2; next} !($2 in a) {a[$2]=1; next}' \
        ${bim} > data/freq/dup_snps_chr${chr}.txt

    # step 1: apply filters and write a clean bfile
    plink \
        --bfile ${BIMDIR}/1000G.EUR.QC.${chr} \
        --mac 5 \
        --geno 0.01 \
        --exclude data/freq/dup_snps_chr${chr}.txt \
        --make-bed \
        --out data/freq/1000G.EUR.QC.${chr}.filtered \
        --silent

    # step 2: compute allele frequency on the filtered bfile
    plink \
        --bfile data/freq/1000G.EUR.QC.${chr}.filtered \
        --freq \
        --out data/freq/1000G.EUR.QC.${chr} \
        --silent
done

# --------------------------------------------------------------------------
# Generate LD scores for custom annotations
# --thin-annot: only output LD scores for annotation columns (not all SNPs)
# --print-snp:  restrict output SNPs to HapMap3 list for consistency with
#               the baseline LD scores (pattern from prset_ldsc.nf)
# --------------------------------------------------------------------------
ml ldsc/1.0.1

for chr in {21..22}; do
    ldsc.py \
        --l2 \
        --bfile ${BIMDIR}/1000G.EUR.QC.${chr} \
        --ld-wind-cm 1 \
        --annot annot/multipathway.${chr}.annot.gz \
        --thin-annot \
        --print-snp ${hm3_snpids} \
        --out ldscores/multipathway.${chr}

    ldsc.py \
        --l2 \
        --bfile ${BIMDIR}/1000G.EUR.QC.${chr} \
        --ld-wind-cm 1 \
        --annot annot/othergenes.${chr}.annot.gz \
        --thin-annot \
        --print-snp ${hm3_snpids} \
        --out ldscores/othergenes.${chr}
done

# --------------------------------------------------------------------------
# Run partitioned heritability
# --------------------------------------------------------------------------
ldsc.py \
    --h2 ibd.sumstats.gz \
    --ref-ld-chr ${BASELINE_LDSCORES},ldscores/multipathway.,ldscores/othergenes. \
    --frqfile-chr data/freq/1000G.EUR.QC. \
    --w-ld-chr /sc/arion/projects/psychgen/projects/prs/sample_overlap/analysis/prepare/data/reference/weight- \
    --overlap-annot \
    --print-coefficients \
    --print-delete-vals \
    --out results/ibd.multipathway_partitioned_h2