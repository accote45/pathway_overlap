## convert GMT file to GSA-MiXer input

library(tidyverse)
library(ggplot2)
library(dplyr)
library(dorothea)
library(GSA)
library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)
gmt_path <- args[1]
gtf_path <- args[2]

# 1. Load pathway data
dat <- GSA.read.gmt(gmt_path)
path.list <- dat$genesets
names(path.list) <- dat$geneset.names

# Convert the path list to a data table
genes_long <- rbindlist(lapply(names(path.list), function(name) {
  data.table(value = path.list[[name]], name = name)
}))

# create gene files for GSA-MiXer
# 1. Full gene list
full_gene <- unique(genes_long$value)

# 2. Retrieve hg19 positions for Ensembl IDs from local GTF
gtf <- import(gtf_path)

# Filter for 'gene' features
gtf_genes <- gtf[gtf$type == "gene"]

# Extract relevant columns
gene_annot <- data.frame(
  GENE = mcols(gtf_genes)$gene_id,
  CHR = as.character(seqnames(gtf_genes)),
  FROM = start(gtf_genes),
  TO = end(gtf_genes)
)

# Filter for standard chromosomes
gene_annot <- gene_annot %>%
  dplyr::filter(CHR %in% c(as.character(1:22), "X", "Y"))

# Join with your gene list
gene_coords <- gene_annot %>%
  dplyr::filter(GENE %in% unique(genes_long$value))

# 3. Merge pathway info
setnames(genes_long, c("value", "name"), c("GENE", "GO"))
genes_merged <- genes_long %>% inner_join(gene_coords, by = "GENE")

# 4. Create baseline file (coding_genes only)
baseline <- gene_coords %>%
  dplyr::mutate(GO = "coding_genes") %>%
  dplyr::select(GO, GENE, CHR, FROM, TO)
write_tsv(baseline, "baseline.txt")

# 5. Create full,gene file (coding_genes + gene name)
full_gene <- as_tibble(gene_coords) %>%
  dplyr::select(GENE, CHR, FROM, TO)

full_gene_rows <- bind_rows(
  full_gene %>% dplyr::mutate(GO = "coding_genes"),
  full_gene %>% dplyr::mutate(GO = GENE)
) %>% dplyr::select(GO, GENE, CHR, FROM, TO)

write_tsv(full_gene_rows, "full_gene.txt")

# 6. Create full,gene set file (coding_genes + gene name + pathway name)
full_gene_set_rows <- bind_rows(
  genes_merged %>% dplyr::mutate(GO = "coding_genes"),
  genes_merged %>% dplyr::mutate(GO = GENE),
  genes_merged %>% dplyr::mutate(GO = GO)
) %>% dplyr::select(GO, GENE, CHR, FROM, TO)

write_tsv(full_gene_set_rows, "full_gene_set.txt")