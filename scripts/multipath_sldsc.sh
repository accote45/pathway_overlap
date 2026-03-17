# references:
https://github.com/bulik/ldsc/wiki/Partitioned-Heritability?utm_source
https://github.com/bulik/ldsc?utm_source
https://hgen471.hakyimlab.org/post/2022/02/22/partition-heritability/#partitioned-heritability

# model:
## the χ2 association statistic for a given SNP includes the effects of all SNPs tagged by this SNP. Thus, for a polygenic trait, SNPs with a high LD score will have higher χ2 statistics on average than SNPs with a low LD score16. This phenomenon might be driven either by the higher likelihood of these SNPs tagging an individual large effect or their ability to tag multiple weak effects.
## If we partition SNPs into functional categories with different contributions to heritability, then LD to a category that is enriched for heritability will increase the χ2 statistic of a SNP more than LD to a category that does not contribute to heritability.
## our method determines that a category of SNPs is enriched for heritability if SNPs with high LD to that category have higher χ2 statistics than SNPs with low LD to that category.


#### checks:
# When you implement the multipathway vs other genes annotations, the LDSC tutorial assumes you’re adding your annotations to the baseline model, not replacing it. That step is critical — otherwise the enrichment can just reflect gene density or LD structure.

# Check whether mean chi2 is above ~1.02. If no, this means there is very little polygenic signal for LDSC to work with.

# A. By default, ldsc reports h2 on the observed scale. Most publications report h2 on the liability scale, since this allows for comparison across studies of the same disease with different proportions of cases and controls. You can use the --samp-prev and --pop-prev flags to convert to liability scale h2 and genetic covariance
	# observed scale = proportion of phenotypic variance in the trait due to genetic variation





#### all in /sc/arion/projects/psychgen/cotea02_prset/pathway_overlap/exploratory/sldsc

# munge GWAS sumstats for LDSC
 ml ldsc/1.0.1
 
snplist=/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/results/ldsc/resource/eur_w_ld_chr/w_hm3.snplist
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

# get list of multipathway genes
ml R

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

writeLines(multipathway_genes, "multipathway_genes.txt")
other_genes <- names(counts[counts < threshold])
writeLines(other_genes, "other_genes.txt")
cat(sprintf("%d other genes\n", length(other_genes)))


# reformat gene coordinates file
library(data.table)

coords <- fread(
  "data/msigdbgenes.regions",
  col.names = c("ensembl_id", "CHR", "START", "END", "strand", "GENE"),
  colClasses = c("character", "character", "integer", "integer", "character", "character")
)

# Keep only autosomes + sex chromosomes
coords <- coords[CHR %in% c(as.character(1:22), "X", "Y")]

# Write in LDSC make_annot.py expected format
fwrite(
  coords[, .(GENE = ensembl_id, CHR, START, END)],
  file = "data/gene_coords.txt",
  sep = "\t"
)



# multipathway annotation
ml bedtools
for chr in {1..22}; do
make_annot.py \
    --gene-set-file data/multipathway_genes.txt \
    --gene-coord-file data/gene_coords.txt \
    --windowsize 35000 \
    --bimfile /sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/gsa_mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
    --annot-file annot/multipathway.${chr}.annot.gz
done

# other genes annotation
for chr in {1..22}; do
make_annot.py \
    --gene-set-file data/gene_lists/other_genes.txt \
    --gene-coord-file data/gene_coords/gene_coords.txt \
    --windowsize 35000 \
    --bimfile /sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/gsa_mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
    --annot-file annot/othergenes.${chr}.annot.gz
done










# multipathway LD scores
for chr in {1..22}; do
ldsc.py \
    --l2 \
    --bfile /sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/gsa_mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
    --ld-wind-cm 1 \
    --annot annot/multipathway.${chr}.annot.gz \
    --thin-annot \
    --out ldscores/multipathway.${chr}
done

# other genes LD scores
for chr in {1..22}; do
  python ldsc/ldsc.py \
    --l2 \
    --bfile /sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/gsa_mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
    --ld-wind-cm 1 \
    --annot annot/othergenes.${chr}.annot.gz \
    --thin-annot \
    --out ldscores/othergenes.${chr}
done

## merge custom annotations with baseline LD 
## compute LD scores for the merged baseline LD + custom file
for chr in {1..22}; do
  python ldsc/ldsc.py \
    --l2 \
    --bfile /sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/gsa_mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
    --ld-wind-cm 1 \
    --annot annot/custom_baselineLD.${chr}.annot.gz \
    --out ldscores/custom_baselineLD.${chr}
done

#munge GWAS sumstats
python ldsc/munge_sumstats.py \
  --sumstats data/sumstats/trait_gwas.txt \
  --N 250000 \
  --out data/sumstats/trait \
  --merge-alleles data/ref/w_hm3.snplist


## run partitioned heritability
python ldsc/ldsc.py \
  --h2 data/sumstats/trait.sumstats.gz \
  --ref-ld-chr ldscores/custom_baselineLD. \
  --frqfile-chr data/ref/1000G_frq/1000G.EUR.QC. \
  --w-ld-chr data/ref/weights_hm3_no_hla/weights.hm3_noMHC. \
  --overlap-annot \
  --print-coefficients \
  --print-delete-vals \
  --out results/trait.multipathway_partitioned_h2