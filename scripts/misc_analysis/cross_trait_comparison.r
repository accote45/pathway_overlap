#!/usr/bin/env Rscript

# =======================================================
# Cross-Trait Pathway Comparison Script
# =======================================================
# This script analyzes how top pathways are shared across traits
# and how ranking differs between randomization methods

library(tidyverse)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(UpSetR)

# Define constants
TOP_N <- 20 # Number of top pathways to examine
BASE_RESULT_DIR <- "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results"
OUTPUT_DIR <- "cross_trait_comparison"
METHODS <- c("original", "birewire", "keeppathsize")
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# Function to load all empirical p-value files for a given method
load_empirical_results <- function(method = "birewire") {
  # Determine correct directory based on method
  result_dir <- file.path(BASE_RESULT_DIR, "empirical_pvalues", paste0("magma_", method))
  
  # Get all result files
  files <- list.files(path = result_dir, pattern = "_empirical_pvalues.txt$", full.names = TRUE, recursive=TRUE)
  
  # Initialize a list to store data for each trait
  results_by_trait <- list()
  
  # Process each file
  for (file in files) {
    # Extract trait name from filename
    trait <- basename(file) %>%
      str_replace("_magma_.*$", "")
    
    cat("Loading", trait, "results for", method, "method\n")
    
    # Read the data
    tryCatch({
      data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      
      # Store with trait and method information
      data <- data %>%
        mutate(trait = trait, 
               method = method,
               rank = rank(P)) %>%
        arrange(P)
      
      results_by_trait[[trait]] <- data
    }, error = function(e) {
      cat("Error loading", file, ":", e$message, "\n")
    })
  }
  
  # Combine all traits into a single dataframe
  all_results <- bind_rows(results_by_trait)
  return(all_results)
}

# Load results for all methods
cat("Loading empirical p-value results...\n")
original_results <- load_empirical_results("birewire") %>%
  select(FULL_NAME, P, trait) %>%
  rename(original_P = P, original_empP = P)

birewire_results <- load_empirical_results("birewire") %>%
  select(FULL_NAME, empirical_pval, trait) %>%
  rename(birewire_empP = empirical_pval)

keeppath_results <- load_empirical_results("keeppathsize") %>%
  select(FULL_NAME,empirical_pval, trait) %>%
  rename(keeppath_empP = empirical_pval)

# Merge all results
all_results <- original_results %>%
  left_join(birewire_results, by = c("FULL_NAME", "trait")) %>%
  left_join(keeppath_results, by = c("FULL_NAME", "trait"))

# Write combined results to file
write.csv(all_results, file.path(OUTPUT_DIR, "combined_pathway_results.csv"), row.names = FALSE)

# Get list of all traits
all_traits <- unique(all_results$trait)
cat("Found", length(all_traits), "traits:", paste(all_traits, collapse = ", "), "\n")

# ===============================
# 1. Get top pathways by method
# ===============================

get_top_pathways <- function(data, method = "original", n = TOP_N) {
  # Select the appropriate p-value column
  p_col <- paste0(method, "_empP")
  
  # For each trait, get top pathways
  top_sets <- data %>%
    group_by(trait) %>%
    arrange(.by_expression = !!sym(p_col)) %>%
    slice_head(n = n) %>%
    ungroup() %>%
    select(trait, FULL_NAME, !!sym(p_col))
  
  return(top_sets)
}

# Get top pathways for each method
top_original <- get_top_pathways(all_results, "original")
top_birewire <- get_top_pathways(all_results, "birewire")
top_keeppath <- get_top_pathways(all_results, "keeppath")

# ===============================
# 2. Calculate pathway overlaps
# ===============================

# For each trait, calculate overlap between methods
trait_method_overlap <- data.frame()
for (trait in all_traits) {
  cat("Calculating overlaps for trait:", trait, "\n")
  
  # Get top pathways for this trait by each method
  orig_sets <- top_original %>% filter(trait == !!trait) %>% pull(FULL_NAME)
  bw_sets <- top_birewire %>% filter(trait == !!trait) %>% pull(FULL_NAME)
  kp_sets <- top_keeppath %>% filter(trait == !!trait) %>% pull(FULL_NAME)
  
  # Calculate overlaps
  orig_bw_overlap <- length(intersect(orig_sets, bw_sets))
  orig_kp_overlap <- length(intersect(orig_sets, kp_sets))
  bw_kp_overlap <- length(intersect(bw_sets, kp_sets))
  all_overlap <- length(Reduce(intersect, list(orig_sets, bw_sets, kp_sets)))
  
  # Add to result dataframe
  trait_method_overlap <- rbind(trait_method_overlap, data.frame(
    trait = trait,
    orig_bw_overlap = orig_bw_overlap,
    orig_kp_overlap = orig_kp_overlap,
    bw_kp_overlap = bw_kp_overlap,
    all_overlap = all_overlap,
    unique_orig = length(setdiff(orig_sets, union(bw_sets, kp_sets))),
    unique_bw = length(setdiff(bw_sets, union(orig_sets, kp_sets))),
    unique_kp = length(setdiff(kp_sets, union(orig_sets, bw_sets))),
    orig_bw_pct = orig_bw_overlap / TOP_N * 100,
    orig_kp_pct = orig_kp_overlap / TOP_N * 100,
    bw_kp_pct = bw_kp_overlap / TOP_N * 100,
    all_pct = all_overlap / TOP_N * 100
  ))
}

# Save results
write.csv(trait_method_overlap, file.path(OUTPUT_DIR, "method_overlap_by_trait.csv"), row.names = FALSE)

# ===============================
# 3. Cross-trait sharing analysis
# ===============================

# Function to calculate cross-trait overlaps for a method
calc_cross_trait_overlap <- function(top_sets, method) {
  # Create a matrix to hold overlaps
  traits <- unique(top_sets$trait)
  n_traits <- length(traits)
  overlap_matrix <- matrix(0, nrow = n_traits, ncol = n_traits)
  rownames(overlap_matrix) <- traits
  colnames(overlap_matrix) <- traits
  
  # For each trait pair, calculate overlap
  for (i in 1:n_traits) {
    for (j in 1:n_traits) {
      if (i == j) {
        overlap_matrix[i, j] <- TOP_N  # Same trait has TOP_N overlap
      } else {
        sets_i <- top_sets %>% filter(trait == traits[i]) %>% pull(FULL_NAME)
        sets_j <- top_sets %>% filter(trait == traits[j]) %>% pull(FULL_NAME)
        overlap_matrix[i, j] <- length(intersect(sets_i, sets_j))
      }
    }
  }
  
  return(overlap_matrix)
}

# Calculate cross-trait overlaps for each method
cross_trait_orig <- calc_cross_trait_overlap(top_original, "original")
cross_trait_bw <- calc_cross_trait_overlap(top_birewire, "birewire")
cross_trait_kp <- calc_cross_trait_overlap(top_keeppath, "keeppath")

# Save overlap matrices
write.csv(cross_trait_orig, file.path(OUTPUT_DIR, "cross_trait_overlap_original.csv"))
write.csv(cross_trait_bw, file.path(OUTPUT_DIR, "cross_trait_overlap_birewire.csv"))
write.csv(cross_trait_kp, file.path(OUTPUT_DIR, "cross_trait_overlap_keeppath.csv"))

# ===============================
# 4. Visualizations 
# ===============================

# 1. Method Overlap Heatmap
plot_method_overlap <- function() {
  long_data <- trait_method_overlap %>%
    pivot_longer(cols = c(orig_bw_pct, orig_kp_pct, bw_kp_pct),
                 names_to = "comparison", 
                 values_to = "percent_overlap") %>%
    mutate(comparison = case_when(
      comparison == "orig_bw_pct" ~ "Original vs BireWire",
      comparison == "orig_kp_pct" ~ "Original vs KeepPath",
      comparison == "bw_kp_pct" ~ "BireWire vs KeepPath"
    ))
  
  # Create the plot
  p <- ggplot(long_data, aes(x = reorder(trait, percent_overlap), y = comparison, fill = percent_overlap)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.1f%%", percent_overlap)), color = "black", size = 3) +
    scale_fill_gradient(low = "white", high = "steelblue", name = "Overlap %") +
    labs(title = paste("Overlap in Top", TOP_N, "Pathways Between Methods"),
         x = "Trait", y = "Method Comparison") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save the plot
  ggsave(file.path(OUTPUT_DIR, "method_overlap_heatmap.pdf"), p, width = 10, height = 6)
  return(p)
}

# 2. Cross-trait Overlap Heatmaps
plot_cross_trait_heatmap <- function(overlap_matrix, method) {
  # Calculate percentage based on total possible unique pathways (TOP_N * 2)
  # For the diagonal (same trait), the denominator remains TOP_N since there can only be TOP_N unique pathways
  percent_matrix <- matrix(0, nrow = nrow(overlap_matrix), ncol = ncol(overlap_matrix))
  
  for (i in 1:nrow(overlap_matrix)) {
    for (j in 1:ncol(overlap_matrix)) {
      if (i == j) {
        # Same trait - maximum overlap is TOP_N (100%)
        percent_matrix[i, j] <- 100
      } else {
        # Different traits - calculate as percentage of possible unique pathways
        percent_matrix[i, j] <- overlap_matrix[i, j] / (TOP_N * 2) * 100
      }
    }
  }
  
  # Create annotation for diagonal (self-overlap)
  annotation_col <- data.frame(Trait = rownames(overlap_matrix))
  rownames(annotation_col) <- rownames(overlap_matrix)
  
  # Create the heatmap
  pdf(file.path(OUTPUT_DIR, paste0("cross_trait_heatmap_", method, ".pdf")), width = 10, height = 8)
  pheatmap(
    percent_matrix,
    main = paste("Cross-Trait Overlap in Top", TOP_N, "Pathways -", toupper(method)),
    color = colorRampPalette(c("white", "steelblue", "darkblue"))(50),
    display_numbers = TRUE,
    number_format = "%.1f%%",
    fontsize_number = 7,
    angle_col = 45,
    clustering_method = "ward.D2"
  )
  dev.off()
  
  # Also save the raw overlap matrix with counts
  pdf(file.path(OUTPUT_DIR, paste0("cross_trait_counts_", method, ".pdf")), width = 10, height = 8)
  pheatmap(
    overlap_matrix,
    main = paste("Number of Shared Pathways -", toupper(method)),
    color = colorRampPalette(c("white", "steelblue", "darkblue"))(50),
    display_numbers = TRUE,
    fontsize_number = 7,
    angle_col = 45,
    clustering_method = "ward.D2"
  )
  dev.off()
}

# 3. Create UpSet plots to visualize pathway sharing
plot_upset <- function(method_name) {
  # Get the data for the specified method
  if (method_name == "original") {
    top_data <- top_original
  } else if (method_name == "birewire") {
    top_data <- top_birewire
  } else if (method_name == "keeppath") {
    top_data <- top_keeppath
  } else {
    stop("Invalid method name")
  }
  
  # Get all unique pathways across traits
  all_pathways <- unique(top_data$FULL_NAME)
  
  # Create a binary matrix for UpSetR
  binary_data <- matrix(0, nrow = length(all_pathways), ncol = length(all_traits))
  colnames(binary_data) <- all_traits
  rownames(binary_data) <- all_pathways
  
  # Fill the matrix
  for (i in 1:nrow(top_data)) {
    pathway <- top_data$FULL_NAME[i]
    trait <- top_data$trait[i]
    binary_data[pathway, trait] <- 1
  }
  
  # Convert to data frame for UpSetR
  binary_df <- as.data.frame(binary_data)
  
  # Create and save UpSet plot
  pdf(file.path(OUTPUT_DIR, paste0("upset_plot_", method_name, ".pdf")), width = 12, height = 8)
  upset(binary_df, nsets = length(all_traits), 
        sets = all_traits, 
        mainbar.y.label = paste("Pathway Overlaps -", toupper(method_name)),
        sets.x.label = "Pathways per Trait",
        text.scale = 1.2,
        point.size = 3,
        line.size = 1)
  dev.off()
}

# Generate all visualizations
cat("Generating visualizations...\n")
plot_method_overlap()
plot_cross_trait_heatmap(cross_trait_orig, "original")
plot_cross_trait_heatmap(cross_trait_bw, "birewire")
plot_cross_trait_heatmap(cross_trait_kp, "keeppath")
plot_upset("original")
plot_upset("birewire")
plot_upset("keeppath")

# ===============================
# 5. Statistical Comparison
# ===============================

# Compare cross-trait overlap distributions between methods
compare_overlap_distributions <- function() {
  # Extract upper triangle values (excluding diagonal)
  get_upper_tri <- function(matrix) {
    upper <- matrix[upper.tri(matrix)]
    return(upper)
  }
  
  orig_overlaps <- get_upper_tri(cross_trait_orig)
  bw_overlaps <- get_upper_tri(cross_trait_bw)
  kp_overlaps <- get_upper_tri(cross_trait_kp)
  
  # Statistical tests
  orig_vs_bw <- wilcox.test(orig_overlaps, bw_overlaps, paired = TRUE)
  orig_vs_kp <- wilcox.test(orig_overlaps, kp_overlaps, paired = TRUE)
  bw_vs_kp <- wilcox.test(bw_overlaps, kp_overlaps, paired = TRUE)
  
  # Create summary dataframe
  summary_stats <- data.frame(
    comparison = c("Original vs BireWire", "Original vs KeepPath", "BireWire vs KeepPath"),
    p_value = c(orig_vs_bw$p.value, orig_vs_kp$p.value, bw_vs_kp$p.value),
    mean_diff = c(
      mean(orig_overlaps) - mean(bw_overlaps),
      mean(orig_overlaps) - mean(kp_overlaps),
      mean(bw_overlaps) - mean(kp_overlaps)
    ),
    median_diff = c(
      median(orig_overlaps) - median(bw_overlaps),
      median(orig_overlaps) - median(kp_overlaps),
      median(bw_overlaps) - median(kp_overlaps)
    )
  )
  
  write.csv(summary_stats, file.path(OUTPUT_DIR, "overlap_comparison_stats.csv"), row.names = FALSE)
  
  # Create boxplot
  overlap_data <- data.frame(
    method = c(rep("Original", length(orig_overlaps)), 
               rep("BireWire", length(bw_overlaps)), 
               rep("KeepPath", length(kp_overlaps))),
    overlap = c(orig_overlaps, bw_overlaps, kp_overlaps)
  )
  
  p <- ggplot(overlap_data, aes(x = method, y = overlap, fill = method)) +
    geom_boxplot() +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Distribution of Cross-Trait Pathway Sharing",
         subtitle = paste("Top", TOP_N, "pathways by method"),
         x = "Method", y = "Number of shared pathways") +
    theme_minimal() +
    theme(legend.position = "none")
  
  ggsave(file.path(OUTPUT_DIR, "cross_trait_overlap_boxplot.pdf"), p, width = 8, height = 6)
  
  return(summary_stats)
}

overlap_stats <- compare_overlap_distributions()
print(overlap_stats)

cat("Analysis complete! Results saved in", OUTPUT_DIR, "directory\n")