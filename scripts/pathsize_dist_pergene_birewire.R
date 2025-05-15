
## determine distribution of genes across random pathways

library(data.table)
library(tidyverse)
library(GSA)

setwd('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/randomized_gene_sets/random_birewire')

# Read pathway size file (same for all GMTs)
pathway_sizes <- read.table("/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/results/pathwaydb_enrichment_OLD/msigdball_pathway_size.txt") %>% rownames_to_column("pathway_name")

# List all GMT files
gmt_files <- list.files(pattern = "\\.gmt$")

# For each GMT, get gene-pathway pairs and join with pathway size
all_gene_pathway_sizes <- rbindlist(lapply(gmt_files, function(gmt_file) {
  dat <- GSA.read.gmt(gmt_file)
  path.list <- dat$genesets
  names(path.list) <- dat$geneset.names
  gene_pathway <- rbindlist(lapply(names(path.list), function(name) {
    data.table(gene = path.list[[name]], pathway = name)
  }))
  # Join with pathway size
  merge(gene_pathway, pathway_sizes, by.x = "pathway", by.y = "pathway_name")
}), use.names = TRUE)

# Example: get the distribution of pathway sizes for a specific gene
gene_of_interest <- "GENE_SYMBOL"
gene_sizes <- all_gene_pathway_sizes[gene == gene_of_interest, pathway_size]

# Plot distribution for a specific gene
ggplot(data.frame(pathway_size = gene_sizes), aes(x = pathway_size)) +
  geom_histogram(binwidth = 1) +
  ggtitle(paste("Distribution of pathway sizes for", gene_of_interest))

# Or summarize for all genes
gene_pathway_summary <- all_gene_pathway_sizes[, .(
  mean_size = mean(pathway_size),
  median_size = median(pathway_size),
  sd_size = sd(pathway_size),
  n = .N
), by = gene]

fwrite(gene_pathway_summary, "gene_pathway_size_summary.csv")