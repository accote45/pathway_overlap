library(data.table)
library(tidyverse)
library(ggplot2)
library(GSA)

# Function to read GMT files from a directory
read_gmt_files <- function(path) {
  setwd(path)
  files <- list.files(pattern="*gmt")
  ldf <- lapply(files, function(f) tryCatch(read.table(f, header=F), error=function(e) NULL))
  ldf <- ldf[!sapply(ldf, is.null)]  # Remove any NULL entries
  return(ldf)
}

# Complete analysis for your hypothesis
test_cooccurrence_hypothesis <- function(random_pathways, magma_results, trait = "height") {
  
  # Load data
  pathways <- read_gmt_files(random_pathways)  # Your 1000 random GMTs
  magma_genes <- magma_results[TRAIT == trait]
  
  # PART 1: Test size effect on any gene pair
  size_effect <- test_pathway_size_effect(pathways)
  
  # PART 2: Test abundance effect  
  abundance_effect <- test_gene_abundance_effect(pathways, magma_genes)
  
  # PART 3: Test combined effect (your main hypothesis)
  combined_effect <- test_combined_inflation_effect(pathways, magma_genes)
  
  # PART 4: Test gene set enrichment inflation
  gene_set_effect <- test_gene_set_inflation_effect(pathways, magma_genes)
  
  return(list(
    size_effect = size_effect,
    abundance_effect = abundance_effect, 
    combined_effect = combined_effect,
    gene_set_effect = gene_set_effect
  ))
}

# PART 1: Pathway size effect
test_pathway_size_effect <- function(pathways) {
  
  pathway_stats <- data.table(
    pathway_id = 1:length(pathways),
    size = lengths(pathways),
    n_pairs = sapply(lengths(pathways), function(x) choose(x, 2))
  )
  
  # Bin by size
  pathway_stats[, size_bin := cut(size, 
                                 breaks = c(0, 20, 50, 100, 200, Inf),
                                 labels = c("Very Small", "Small", "Medium", "Large", "Very Large"))]
  
  # Expected co-occurrence probability by size
  N_total <- length(unique(unlist(pathways)))
  pathway_stats[, expected_cooccur_prob := (size/N_total) * ((size-1)/(N_total-1))]
  
  # Test: Does larger pathways = higher co-occurrence?
  size_correlation <- cor.test(pathway_stats$size, pathway_stats$expected_cooccur_prob)
  
  # Plot
  p1 <- ggplot(pathway_stats, aes(x = size, y = expected_cooccur_prob)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm") +
    scale_y_log10() +
    labs(title = "Hypothesis Part 1: Pathway Size Effect",
         x = "Pathway Size", y = "Expected Co-occurrence Probability (log)")
  
  return(list(stats = pathway_stats, correlation = size_correlation, plot = p1))
}

# PART 2: Gene abundance effect
test_gene_abundance_effect <- function(pathways, magma_genes) {
  
  # Calculate gene frequencies
  all_genes <- unlist(pathways)
  gene_freq <- table(all_genes)
  
  # Merge with MAGMA results
  gene_stats <- data.table(
    GENE = names(gene_freq),
    frequency = as.numeric(gene_freq)
  )
  gene_stats <- merge(gene_stats, magma_genes[, .(GENE, ZSTAT)], by = "GENE")
  
  # Categorize by frequency
  gene_stats[, freq_category := cut(frequency,
                                   breaks = c(0, 10, 50, 100, 500, Inf),
                                   labels = c("Rare", "Uncommon", "Common", "Frequent", "Hub"))]
  
  # Test: Do high-frequency genes have higher associations?
  freq_correlation <- cor.test(gene_stats$frequency, abs(gene_stats$ZSTAT))
  
  # Plot
  p2 <- ggplot(gene_stats, aes(x = frequency, y = abs(ZSTAT), color = freq_category)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm") +
    scale_x_log10() +
    labs(title = "Hypothesis Part 2: Gene Abundance Effect",
         x = "Gene Frequency (log)", y = "MAGMA |Z-stat|")
  
  return(list(stats = gene_stats, correlation = freq_correlation, plot = p2))
}

# PART 3: Combined inflation effect (your main hypothesis)
test_combined_inflation_effect <- function(pathways, magma_genes) {
  
  # For each pathway, calculate:
  # 1. Size
  # 2. Number of hub genes
  # 3. Number of significant hub genes  
  # 4. Expected "inflation score"
  
  gene_freq <- table(unlist(pathways))
  hub_genes <- names(gene_freq)[gene_freq >= 50]  # Define hub threshold
  
  # Significant genes (|Z| > 1.96)
  sig_genes <- magma_genes[abs(ZSTAT) > 1.96, GENE]
  sig_hub_genes <- intersect(hub_genes, sig_genes)
  
  pathway_inflation <- data.table(
    pathway_id = 1:length(pathways),
    size = lengths(pathways),
    n_hub_genes = sapply(pathways, function(p) sum(p %in% hub_genes)),
    n_sig_hub_genes = sapply(pathways, function(p) sum(p %in% sig_hub_genes)),
    hub_proportion = sapply(pathways, function(p) mean(p %in% hub_genes)),
    sig_hub_proportion = sapply(pathways, function(p) mean(p %in% sig_hub_genes))
  )
  
  # Calculate "inflation risk score"
  pathway_inflation[, inflation_risk := size * n_sig_hub_genes * hub_proportion]
  
  # Test correlations
  size_hub_cor <- cor.test(pathway_inflation$size, pathway_inflation$n_hub_genes)
  size_sighub_cor <- cor.test(pathway_inflation$size, pathway_inflation$n_sig_hub_genes)
  
  # Plot the smoking gun
  p3 <- ggplot(pathway_inflation, aes(x = size, y = n_sig_hub_genes)) +
    geom_point(aes(color = hub_proportion), alpha = 0.7) +
    geom_smooth(method = "lm") +
    scale_color_viridis_c(name = "Hub Gene\nProportion") +
    labs(title = "Hypothesis Part 3: Combined Inflation Effect",
         x = "Pathway Size", 
         y = "Number of Significant Hub Genes",
         subtitle = "Color = Hub gene proportion")
  
  return(list(
    stats = pathway_inflation,
    correlations = list(size_hub = size_hub_cor, size_sighub = size_sighub_cor),
    plot = p3
  ))
}

# PART 4: Gene set enrichment inflation effect
test_gene_set_inflation_effect <- function(pathways, magma_genes) {
  
  gene_freq <- table(unlist(pathways))
  hub_genes <- names(gene_freq)[gene_freq >= 50]
  
  # Define significant genes by MAGMA Z-stat
  sig_genes <- magma_genes[abs(ZSTAT) > 1.96, GENE]
  strong_sig_genes <- magma_genes[abs(ZSTAT) > 2.58, GENE]  # p < 0.01
  
  # Hub genes that are also trait-associated
  sig_hub_genes <- intersect(hub_genes, sig_genes)
  strong_sig_hub_genes <- intersect(hub_genes, strong_sig_genes)
  
  pathway_inflation <- data.table(
    pathway_id = 1:length(pathways),
    size = lengths(pathways),
    n_hub_genes = sapply(pathways, function(p) sum(p %in% hub_genes)),
    n_sig_genes = sapply(pathways, function(p) sum(p %in% sig_genes)),
    n_sig_hub_genes = sapply(pathways, function(p) sum(p %in% sig_hub_genes)),
    n_strong_sig_hub_genes = sapply(pathways, function(p) sum(p %in% strong_sig_hub_genes)),
    
    # Collective signal metrics
    mean_abs_zstat = sapply(pathways, function(p) {
      pathway_genes <- intersect(p, magma_genes$GENE)
      if(length(pathway_genes) > 0) {
        mean(abs(magma_genes[GENE %in% pathway_genes, ZSTAT]), na.rm = TRUE)
      } else NA
    }),
    
    sum_abs_zstat = sapply(pathways, function(p) {
      pathway_genes <- intersect(p, magma_genes$GENE)
      if(length(pathway_genes) > 0) {
        sum(abs(magma_genes[GENE %in% pathway_genes, ZSTAT]), na.rm = TRUE)
      } else 0
    }),
    
    # Hub gene signal contribution
    hub_signal_contribution = sapply(pathways, function(p) {
      hub_in_pathway <- intersect(p, sig_hub_genes)
      if(length(hub_in_pathway) > 0) {
        sum(abs(magma_genes[GENE %in% hub_in_pathway, ZSTAT]), na.rm = TRUE)
      } else 0
    })
  )
  
  # Calculate proportion of signal from hub genes
  pathway_inflation[, hub_signal_proportion := hub_signal_contribution / sum_abs_zstat]
  pathway_inflation[hub_signal_proportion > 1, hub_signal_proportion := 1]  # Cap at 100%
  
  # Test key correlations for gene set enrichment inflation
  size_meansignal_cor <- cor.test(pathway_inflation$size, pathway_inflation$mean_abs_zstat)
  size_totalsignal_cor <- cor.test(pathway_inflation$size, pathway_inflation$sum_abs_zstat)
  size_hubsignal_cor <- cor.test(pathway_inflation$size, pathway_inflation$hub_signal_contribution)
  hub_proportion_signal_cor <- cor.test(pathway_inflation$n_hub_genes, pathway_inflation$sum_abs_zstat)
  
  # Key plot: How pathway size affects total trait signal
  p1 <- ggplot(pathway_inflation, aes(x = size, y = sum_abs_zstat)) +
    geom_point(aes(color = n_sig_hub_genes), alpha = 0.7) +
    geom_smooth(method = "lm") +
    scale_color_viridis_c(name = "Significant\nHub Genes") +
    labs(title = "Gene Set Enrichment Inflation: Total Trait Signal vs Pathway Size",
         x = "Pathway Size", 
         y = "Sum of |MAGMA Z-statistics|",
         subtitle = "Color = Number of significant hub genes")
  
  # Second plot: Hub gene signal contribution
  p2 <- ggplot(pathway_inflation, aes(x = size, y = hub_signal_proportion)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm") +
    labs(title = "Hub Gene Signal Dominance in Large Pathways",
         x = "Pathway Size",
         y = "Proportion of Signal from Hub Genes")
  
  # Third plot: Mean vs total signal (shows accumulation effect)
  p3 <- ggplot(pathway_inflation, aes(x = mean_abs_zstat, y = sum_abs_zstat, color = size)) +
    geom_point(alpha = 0.7) +
    scale_color_viridis_c(name = "Pathway\nSize") +
    labs(title = "Signal Accumulation: Mean vs Total Signal",
         x = "Mean |Z-statistic| per Gene",
         y = "Total |Z-statistic| Sum")
  
  return(list(
    stats = pathway_inflation,
    correlations = list(
      size_mean_signal = size_meansignal_cor,
      size_total_signal = size_totalsignal_cor,
      size_hub_signal = size_hubsignal_cor,
      hub_proportion_signal = hub_proportion_signal_cor
    ),
    plots = list(p1, p2, p3)
  ))
}

# Example usage:
# results <- test_cooccurrence_hypothesis(
#   random_pathways = "/path/to/random/gmt/directory",
#   magma_results = your_magma_data,
#   trait = "height"
# )