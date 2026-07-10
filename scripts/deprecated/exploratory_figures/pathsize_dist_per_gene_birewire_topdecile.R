## determine distribution of genes across random pathways

library(data.table)
library(tidyverse)
library(GSA)
library(parallel)
library(ggpubr)

setwd('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/randomized_gene_sets/random_birewire')

# Read pathway size file (same for all GMTs)
pathway_sizes <- read.table("/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/results/pathwaydb_enrichment_OLD/msigdball_pathway_size.txt") %>% rownames_to_column("pathway_name")

# List all GMT files
gmt_files <- list.files(pattern = "\\.gmt$")
gmt_files <- head(gmt_files, 100)   # first 100 only, representative subset

# For each GMT, get gene-pathway pairs and join with pathway size
gene_pathway_list <- mclapply(gmt_files, function(gmt_file) {
  file_id <- gsub("\\.gmt$", "", gmt_file)
  dat <- GSA.read.gmt(gmt_file)
  names(dat$genesets) <- dat$geneset.names

  # build long table without per-pathway rbindlist
  lens <- lengths(dat$genesets)
  data.table(
    gene    = unlist(dat$genesets, use.names = FALSE),
    pathway = rep(dat$geneset.names, lens),
    file    = file_id
  )
}, mc.cores = detectCores() - 2)

all_gene_pathway_sizes <- rbindlist(gene_pathway_list, use.names = TRUE)
# join sizes once, on the full table (faster than 1000 small joins)
all_gene_pathway_sizes <- merge(
  all_gene_pathway_sizes, pathway_sizes,
  by.x = "pathway", by.y = "pathway_name"
)

# summarize for all genes
gene_pathway_summary <- all_gene_pathway_sizes[, .(
  mean_size = mean(pathway_size),
  median_size = median(pathway_size),
  sd_size = sd(pathway_size),
  n = .N
), by = gene]

fwrite(gene_pathway_summary, "gene_pathway_size_summary_top100.csv")






# determine hub genes
genefreq <- as.data.frame(table(all_gene_pathway_sizes$gene))
n_files <- length(gmt_files)
genefreq$Freq <- genefreq$Freq / n_files

# top decile = genes whose frequency is at/above the 90th percentile
cutoff <- quantile(genefreq$Freq, 0.9)
hub_genes <- genefreq[genefreq$Freq >= cutoff, "Var1"]

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
fwrite(hub_gene_analysis, "hub_gene_proportion_by_pathway_file_decile.csv")

# Calculate summary statistics across files for each pathway
hub_gene_pathway_summary <- hub_gene_analysis[, .(
  mean_hub_proportion = mean(hub_gene_proportion),
  median_hub_proportion = median(hub_gene_proportion),
  sd_hub_proportion = sd(hub_gene_proportion),
  mean_hub_count = mean(hub_gene_count),
  mean_total_genes = mean(total_genes)
), by = total_genes]

fwrite(hub_gene_pathway_summary, "hub_gene_pathway_summary_decile.csv")



# Plot hub gene proportions across pathways

hub_gene_analysis <- fread("/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/randomized_gene_sets/random_birewire/hub_gene_proportion_by_pathway_file_decile.csv")

# too large, subset for 100 randomized files
index <- sample(unique(hub_gene_analysis$file), 100)
temp <- hub_gene_analysis[hub_gene_analysis$file %in% index, ]


# first figure: proportion
png('test.png', width = 7, height = 6, units = "in", res = 300)
ggplot(hub_gene_analysis, aes(x = total_genes, y = hub_gene_proportion)) +
  geom_point(color = "black", alpha = 0.5, size = 1) +  # darker points
  geom_smooth(method = "lm", color = "blue") +
  stat_cor(method = "pearson",
           aes(label = ..r.label..),                 # r only, no p-value
           cor.coef.name = "r",                       # lowercase r notation
           size = 6,                                  # correlation text
           label.x.npc = "right", label.y.npc = "bottom",  # bottom-right
           hjust = 1) +
  labs(x = "Pathway size",
       y = "Proportion of multi-pathway genes") +
  theme_classic(base_size = 16) +                     # clean white background
  theme(axis.line  = element_line(color = "black"),   # outlined black axes
        axis.title = element_text(size = 18),         # axis titles larger than body
        axis.text  = element_text(size = 14, color = "black"))  # tick labels
dev.off()

# second figure: count
png('test2.png', width = 7, height = 6, units = "in", res = 300)
ggplot(hub_gene_analysis, aes(x = total_genes, y = hub_gene_count)) +
  geom_point(color = "black", alpha = 0.5, size = 1) +  # darker points
  geom_smooth(method = "lm", color = "blue") +
  stat_cor(method = "pearson",
           aes(label = ..r.label..),                 # r only, no p-value
           cor.coef.name = "r",                       # lowercase r notation
           size = 6,                                  # correlation text
           label.x.npc = "right", label.y.npc = "bottom",  # bottom-right
           hjust = 1) +
  labs(x = "Pathway size",
       y = "Number of multi-pathway genes") +
  theme_classic(base_size = 16) +                     # clean white background
  theme(axis.line  = element_line(color = "black"),   # outlined black axes
        axis.title = element_text(size = 18),         # axis titles larger than body
        axis.text  = element_text(size = 14, color = "black"))  # tick labels
dev.off()









# =====================================================================
# REAL GENE SET ANALYSIS (single GMT, same figures as above)
# =====================================================================

real_gmt <- "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt"

# Read the real GMT and build gene-pathway long table
dat_real <- GSA.read.gmt(real_gmt)
names(dat_real$genesets) <- dat_real$geneset.names

lens_real <- lengths(dat_real$genesets)
real_gene_pathway <- data.table(
  gene    = unlist(dat_real$genesets, use.names = FALSE),
  pathway = rep(dat_real$geneset.names, lens_real)
)

# Pathway size = number of genes in each pathway (computed directly from the GMT)
real_gene_pathway[, pathway_size := .N, by = pathway]

# summarize for all genes: distribution of pathway sizes each gene belongs to
real_gene_pathway_summary <- real_gene_pathway[, .(
  mean_size   = mean(pathway_size),
  median_size = median(pathway_size),
  sd_size     = sd(pathway_size),
  n           = .N
), by = gene]

fwrite(real_gene_pathway_summary, "real_gene_pathway_size_summary.csv")

# determine hub genes: in a single GMT, gene frequency = number of pathways it belongs to
real_genefreq <- as.data.frame(table(real_gene_pathway$gene))

# top decile = genes whose pathway membership count is at/above the 90th percentile
real_cutoff   <- quantile(real_genefreq$Freq, 0.9)
real_hub_genes <- real_genefreq[real_genefreq$Freq >= real_cutoff, "Var1"]

# For each pathway, calculate total genes, hub gene count, and hub proportion
real_hub_gene_analysis <- real_gene_pathway[, {
  pathway_genes <- unique(gene)
  hub_count     <- sum(pathway_genes %in% real_hub_genes)
  list(
    total_genes         = length(pathway_genes),
    hub_gene_count      = hub_count,
    hub_gene_proportion = hub_count / length(pathway_genes)
  )
}, by = .(pathway)]

fwrite(real_hub_gene_analysis, "real_hub_gene_proportion_by_pathway.csv")

# first figure: proportion
png('real_test.png', width = 7, height = 6, units = "in", res = 300)
ggplot(real_hub_gene_analysis, aes(x = total_genes, y = hub_gene_proportion)) +
  geom_point(color = "black", alpha = 0.5, size = 1) +
  geom_smooth(method = "lm", color = "blue") +
  stat_cor(method = "pearson",
           aes(label = ..r.label..),
           cor.coef.name = "r",
           size = 6,
           label.x.npc = "right", label.y.npc = "bottom",
           hjust = 1) +
  labs(x = "Pathway size",
       y = "Proportion of multi-pathway genes") +
  theme_classic(base_size = 16) +
  theme(axis.line  = element_line(color = "black"),
        axis.title = element_text(size = 18),
        axis.text  = element_text(size = 14, color = "black"))
dev.off()

# second figure: count
png('real_test2.png', width = 7, height = 6, units = "in", res = 300)
ggplot(real_hub_gene_analysis, aes(x = total_genes, y = hub_gene_count)) +
  geom_point(color = "black", alpha = 0.5, size = 1) +
  geom_smooth(method = "lm", color = "blue") +
  stat_cor(method = "pearson",
           aes(label = ..r.label..),
           cor.coef.name = "r",
           size = 6,
           label.x.npc = "right", label.y.npc = "bottom",
           hjust = 1) +
  labs(x = "Pathway size",
       y = "Number of multi-pathway genes") +
  theme_classic(base_size = 16) +
  theme(axis.line  = element_line(color = "black"),
        axis.title = element_text(size = 18),
        axis.text  = element_text(size = 14, color = "black"))
dev.off()
