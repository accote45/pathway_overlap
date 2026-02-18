
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
  # Extract the file identifier (assuming format like "GeneSet.random1.gmt")
  file_id <- gsub("\\.gmt$", "", gmt_file)
  
  dat <- GSA.read.gmt(gmt_file)
  path.list <- dat$genesets
  names(path.list) <- dat$geneset.names
  gene_pathway <- rbindlist(lapply(names(path.list), function(name) {
    data.table(
      gene = path.list[[name]], 
      pathway = name,
      file = file_id   # Add file identifier to each row
    )
  }))
  # Join with pathway size
  merge(gene_pathway, pathway_sizes, by.x = "pathway", by.y = "pathway_name")
}), use.names = TRUE)

# summarize for all genes
gene_pathway_summary <- all_gene_pathway_sizes[, .(
  mean_size = mean(pathway_size),
  median_size = median(pathway_size),
  sd_size = sd(pathway_size),
  n = .N
), by = gene]

fwrite(gene_pathway_summary, "gene_pathway_size_summary.csv")






# determine hub genes (those appearing in more than 100 pathways)
genefreq <- as.data.frame(table(all_gene_pathway_sizes$gene))
genefreq$Freq <- genefreq$Freq/1000
hub_genes <- genefreq[genefreq$Freq > 50, "Var1"]

# For each pathway in each file, calculate:
# 1. Total number of genes
# 2. Number of hub genes
# 3. Proportion of hub genes
hub_gene_analysis <- all_gene_pathway_sizes[, {
  # Get unique genes in this pathway-file combination
  pathway_genes <- unique(gene)
  
  # Count how many are hub genes
  hub_count <- sum(pathway_genes %in% hub_genes)
  
  # Calculate proportion
  list(
    total_genes = length(pathway_genes),
    hub_gene_count = hub_count,
    hub_gene_proportion = hub_count / length(pathway_genes)
  )
}, by = .(file, pathway)]

# Write results
fwrite(hub_gene_analysis, "hub_gene_proportion_by_pathway_file.csv")

# Calculate summary statistics across files for each pathway
hub_gene_pathway_summary <- hub_gene_analysis[, .(
  mean_hub_proportion = mean(hub_gene_proportion),
  median_hub_proportion = median(hub_gene_proportion),
  sd_hub_proportion = sd(hub_gene_proportion),
  mean_hub_count = mean(hub_gene_count),
  mean_total_genes = mean(total_genes)
), by = total_genes]

fwrite(hub_gene_pathway_summary, "hub_gene_pathway_summary.csv")



# Plot hub gene proportions across pathways

hub_gene_analysis <- fread("/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/randomized_gene_sets/random_birewire/hub_gene_proportion_by_pathway_file.csv")

# too large, subset for 100 randomized files
index <- sample(unique(hub_gene_analysis$file), 100)
temp <- hub_gene_analysis[hub_gene_analysis$file %in% index, ]

library(ggplot2)
pdf('test.pdf')
ggplot(temp, aes(x = total_genes,y=hub_gene_proportion)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(x = "Pathway size",
       y = "Proportion of Hub Genes") +
  theme_minimal()
dev.off()

library(ggplot2)
pdf('test.pdf')
ggplot(temp, aes(x = total_genes,y=hub_gene_count)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(x = "Pathway size",
       y = "Number of Hub Genes") +
  theme_minimal()
dev.off()



