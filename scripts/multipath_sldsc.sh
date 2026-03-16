


# munge GWAS sumstats for LDSC

# get list of multipathway genes
# get list of all other genes (in pathway db)
# get gene coordinates for all genes
# need LDSC references



# multipathway annotation
for chr in {1..22}; do
  python ldsc/make_annot.py \
    --gene-set-file data/gene_lists/multipathway_genes.txt \
    --gene-coord-file data/gene_coords/gene_coords.txt \
    --windowsize 100000 \
    --bimfile data/ref/1kg_eur/1000G.EUR.QC.${chr}.bim \
    --annot-file annot/multipathway.${chr}.annot.gz
done

# other genes annotation
for chr in {1..22}; do
  python ldsc/make_annot.py \
    --gene-set-file data/gene_lists/other_genes.txt \
    --gene-coord-file data/gene_coords/gene_coords.txt \
    --windowsize 100000 \
    --bimfile data/ref/1kg_eur/1000G.EUR.QC.${chr}.bim \
    --annot-file annot/othergenes.${chr}.annot.gz
done

# multipathway LD scores
for chr in {1..22}; do
  python ldsc/ldsc.py \
    --l2 \
    --bfile data/ref/1kg_eur/1000G.EUR.QC.${chr} \
    --ld-wind-cm 1 \
    --annot annot/multipathway.${chr}.annot.gz \
    --thin-annot \
    --out ldscores/multipathway.${chr}
done

# other genes LD scores
for chr in {1..22}; do
  python ldsc/ldsc.py \
    --l2 \
    --bfile data/ref/1kg_eur/1000G.EUR.QC.${chr} \
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
    --bfile data/ref/1kg_eur/1000G.EUR.QC.${chr} \
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