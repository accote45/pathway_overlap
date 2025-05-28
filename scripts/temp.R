library(data.table)
library(tidyverse)
library(ggplot2)
library(GSA)
library(parallel)

setwd('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/randomized_gene_sets/random_birewire')

gmt_files <- list.files(pattern="*.gmt")

read_all_gmts_fast <- function(file_list) {
  # Use parallel processing
  n_cores <- detectCores() - 2  # Don't use all cores
  
  all_pathways <- mclapply(file_list, function(gmt) {
    tryCatch({
      dat <- GSA.read.gmt(gmt)
      pathways <- dat$genesets
      names(pathways) <- paste0(dat$geneset.names, "_", tools::file_path_sans_ext(gmt))
      return(pathways)
    }, error = function(e) {
      message(paste("Error reading", gmt, ":", e$message))
      return(NULL)
    })
  }, mc.cores = n_cores)
  
  # Flatten list and remove NULLs
  all_pathways <- unlist(all_pathways, recursive = FALSE)
  all_pathways <- all_pathways[!sapply(all_pathways, is.null)]
  
  return(all_pathways)
}

all_pathways <- read_all_gmts_fast(gmt_files)
saveRDS(all_pathways, file = "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/randomized_gene_sets/random_birewire/all_pathways.rds")



all_pathways <- readRDS("/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/randomized_gene_sets/random_birewire/all_pathways.rds")

# subset full pathway list for pathways of specific sizes (speeds up testing)
subset_pathways <- function(pathways, sizes = c(20, 50, 100, 200, 500, 1000)) {
  subsetted <- pathways[sapply(pathways, function(p) length(p) %in% sizes)]
  return(subsetted)
}

sub_pathways <- subset_pathways(all_pathways)
    # number of pathways for each size
    table(sapply(sub_pathways, length))
    
# downsample sub_pathways to same pathways per size group (take min size group)
downsample_pathways <- function(pathways) {
  size_groups <- split(pathways, sapply(pathways, length))
  n_per_size <- min(table(sapply(pathways,length)))

  downsampled <- lapply(size_groups, function(group) {
    if(length(group) > n_per_size) {
      return(sample(group, n_per_size))
    } else {
      return(group)
    }
  })
  return(unlist(downsampled, recursive = FALSE))
}
downsampled_pathways <- downsample_pathways(sub_pathways)


# read in gene frequencies
  dat <- GSA.read.gmt('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/pathway_db/msigdb/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt')
  path.list <- dat$genesets
  names(path.list) <- dat$geneset.names
  genes_long <- rbindlist(lapply(names(path.list), function(name) {
    data.table(value = path.list[[name]], name = name)
  }))
  # Calculate gene frequency
  gene_freq <- as.data.frame(table(genes_long$value))






# Modified function to use external gene frequency data
test_hub_gene_counts_by_size <- function(pathways, gene_freq_df, hub_threshold = 100) {
  
  # Use external gene frequency data
  colnames(gene_freq_df) <- c("gene", "frequency")  # Standardize column names
  hub_genes <- gene_freq_df$gene[gene_freq_df$frequency >= hub_threshold]
  
  message(paste("Found", length(hub_genes), "hub genes (frequency >=", hub_threshold, ") from external data"))
  
  # Pre-calculate pathway sizes
  pathway_sizes <- lengths(pathways)
  
  # Vectorized hub gene counting
  hub_counts <- vapply(pathways, function(p) sum(p %in% hub_genes), integer(1))
  
  pathway_stats <- data.table(
    pathway_id = seq_along(pathways),
    size = pathway_sizes,
    n_hub_genes = hub_counts,
    hub_proportion = hub_counts / pathway_sizes
  )
  
  # Correlations
  size_hubcount_cor <- cor.test(pathway_stats$size, pathway_stats$n_hub_genes)
  size_hubprop_cor <- cor.test(pathway_stats$size, pathway_stats$hub_proportion)
  
  return(list(
    stats = pathway_stats,
    hub_genes = hub_genes,
    correlations = list(
      size_count = size_hubcount_cor,
      size_proportion = size_hubprop_cor)
  ))
}

master <- test_hub_gene_counts_by_size(downsampled_pathways, gene_freq)

# Plot
pdf('test.pdf')
ggplot(master$stats, aes(x = size, y = hub_proportion)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm") +
    labs(title = "Hub Gene Count vs Pathway Size",
         x = "Pathway Size", y = "Number of Hub Genes",
         subtitle = paste("r =", round(size_hubcount_cor$estimate, 3))) +
    theme_minimal()
dev.off()





# Modified gene pairs function
test_hub_gene_cooccurrence_by_size <- function(pathways, gene_freq_df, hub_threshold = 50) {
  
  # Use external gene frequency data
  colnames(gene_freq_df) <- c("gene", "frequency")
  hub_genes <- gene_freq_df$gene[gene_freq_df$frequency >= hub_threshold]
  
  message(paste("Found", length(hub_genes), "hub genes (frequency >=", hub_threshold, ")"))
  
  # Calculate pathway-level statistics
  pathway_stats <- data.table(
    pathway_id = 1:length(pathways),
    pathway_name = names(pathways),
    size = lengths(pathways),
    n_hub_genes = vapply(pathways, function(p) sum(p %in% hub_genes), integer(1)),
    hub_proportion = vapply(pathways, function(p) mean(p %in% hub_genes), numeric(1))
  )
  
  # Hub gene pairs in each pathway (vectorized)
  pathway_stats[, hub_pairs := vapply(pathways, function(p) {
    hub_in_pathway <- sum(p %in% hub_genes)
    if(hub_in_pathway >= 2) {
      return(choose(hub_in_pathway, 2))
    } else {
      return(0L)
    }
  }, integer(1))]
  
  # Expected hub pairs using external frequencies
  total_genes <- nrow(gene_freq_df)
  hub_prob <- length(hub_genes) / total_genes
  pathway_stats[, expected_hub_pairs := vapply(pathways, function(p) {
    pathway_size <- length(p)
    if(pathway_size >= 2) {
      expected_hub_genes <- pathway_size * hub_prob
      if(expected_hub_genes >= 2) {
        return(choose(expected_hub_genes, 2))
      }
    }
    return(0)
  }, numeric(1))]
  
  # Hub pair density
  pathway_stats[, hub_pair_density := ifelse(n_hub_genes >= 2, hub_pairs / n_hub_genes, 0)]
  
  # Size bins
  pathway_stats[, size_bin := cut(size, 
                                 breaks = c(0, 30, 60, 100, 200, Inf),
                                 labels = c("Very Small", "Small", "Medium", "Large", "Very Large"))]
  
  # Statistical tests
  size_hub_pairs_cor <- cor.test(pathway_stats$size, pathway_stats$hub_pairs)
  size_hub_density_cor <- cor.test(pathway_stats$size, pathway_stats$hub_pair_density)
  size_hub_prop_cor <- cor.test(pathway_stats$size, pathway_stats$hub_proportion)
  
  # Summary by size
  size_summary <- pathway_stats[, .(
    count = .N,
    mean_size = mean(size),
    mean_hub_genes = mean(n_hub_genes),
    mean_hub_pairs = mean(hub_pairs),
    mean_hub_proportion = mean(hub_proportion),
    mean_hub_pair_density = mean(hub_pair_density),
    total_hub_pairs = sum(hub_pairs)
  ), by = size_bin]
  
  # Plots
  p1 <- ggplot(pathway_stats, aes(x = size, y = hub_pairs, color = size_bin)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = "Hub Gene Pairs vs Pathway Size (External Frequencies)",
         x = "Pathway Size", 
         y = "Number of Hub Gene Pairs",
         subtitle = paste("r =", round(size_hub_pairs_cor$estimate, 3))) +
    theme_minimal()
  
  return(list(
    pathway_stats = pathway_stats,
    size_summary = size_summary,
    hub_genes = hub_genes,
    correlations = list(
      size_hub_pairs = size_hub_pairs_cor,
      size_hub_density = size_hub_density_cor,
      size_hub_proportion = size_hub_prop_cor
    ),
    plots = list(p1)
  ))
}

# Usage with your data:
# Read in your gene frequency file (already done)
gene_freq <- as.data.frame(table(genes_long$value))

# Load your subset pathways
all_pathways <- readRDS("/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/randomized_gene_sets/random_birewire/all_pathways.rds")

# Subset to specific sizes
sub_pathways <- subset_pathways(all_pathways, sizes = c(50, 200, 500))
downsampled_pathways <- downsample_pathways(sub_pathways)

# Run analysis with external gene frequencies
simple_results <- test_hub_gene_counts_by_size(downsampled_pathways, gene_freq, hub_threshold = 50)
pairs_results <- test_hub_gene_cooccurrence_by_size(downsampled_pathways, gene_freq, hub_threshold = 50)

# Compare correlations
print("Simple count approach:")
print(paste("Size vs Hub Count:", round(simple_results$correlations$size_count$estimate, 4)))

print("Gene pairs approach:")  
print(paste("Size vs Hub Pairs:", round(pairs_results$correlations$size_hub_pairs$estimate, 4)))

# View plots
simple_results$plot
pairs_results$plots[[1]]

# Summary statistics
print("Hub gene summary by pathway size:")
print(pairs_results$size_summary)
