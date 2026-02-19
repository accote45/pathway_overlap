# ===== Tissue Specificity Correlation Analysis =====
# This script analyzes the correlation between pathway rankings
# and tissue specificity metrics

library(plyr)
library(tidyverse)
library(ggplot2)
library(data.table)
library(GSA)

# Helper function for writing results to CSV
write_results_csv <- function(data, trait, tool_base, prefix="", suffix="", verbose=TRUE) {
  filename <- paste0(trait, "_", tool_base, prefix, ".csv")
  
  write.csv(data, filename, row.names=FALSE)
  if(verbose) cat("Wrote", nrow(data), "rows to", filename, "\n")
  return(filename)
}

# Process command line arguments
args <- commandArgs(trailingOnly = TRUE)

trait <- args[1]
tool_base <- args[2]
birewire_results_file <- args[3]
keeppathsize_results_file <- args[4]
gmt_file <- args[5]
tissue_file <- args[6]

cat("======= Starting Tissue Specificity Correlation Analysis for", trait, "=======\n")

# 1. Load pathway data
cat("Loading pathway data from", gmt_file, "...\n")
dat <- GSA.read.gmt(gmt_file)
path.list <- dat$genesets
names(path.list) <- dat$geneset.names

# Convert the path list to a data table
genes_long <- rbindlist(lapply(names(path.list), function(name) {
  data.table(value = path.list[[name]], name = name)
}))

# 2. Load tissue specificity data
cat("Loading tissue specificity data from", tissue_file, "...\n")
ts <- read.csv(tissue_file)
cat("Loaded tissue data with", nrow(ts), "genes and", ncol(ts) - 1, "tissues\n")

# 3. Calculate pathway-level tissue expression scores
cat("Calculating pathway-level tissue expression scores...\n")
masterfin <- merge(genes_long, ts, by.x="value", by.y="Name")

# Get all tissue columns (excluding gene identifiers)
tissue_cols <- setdiff(colnames(masterfin), c("value", "name"))

# Calculate pathway-level tissue expression metrics (mean only, no median)
pathway_tissue_scores <- masterfin %>% 
  group_by(name) %>% 
  summarise(
    # Calculate mean for each tissue
    across(all_of(tissue_cols), 
           list(
             mean = ~mean(.x, na.rm = TRUE)
           ),
           .names = "{.col}_mean"),
    # Count genes per pathway
    num_genes = n()
  ) %>%
  as.data.frame()

cat("Generated tissue expression profiles for", nrow(pathway_tissue_scores), "pathways\n")

# 4. Load empirical results files
birewire_data <- read.table(birewire_results_file, header = TRUE,fill=TRUE)
# Remove "Base" pathway if present (not a real pathway)
birewire_data <- birewire_data[!(birewire_data$pathway_name=="Base"),]

keeppath_data <- read.table(keeppathsize_results_file, header = TRUE,fill=TRUE)
# Remove "Base" pathway if present (not a real pathway)
keeppath_data <- keeppath_data[!(keeppath_data$pathway_name=="Base"),]

# Normalize column names to ensure consistency
if(!"pathway_name" %in% colnames(birewire_data) && "FULL_NAME" %in% colnames(birewire_data)) {
  birewire_data$pathway_name <- birewire_data$FULL_NAME
}
if(!"pathway_name" %in% colnames(keeppath_data) && "FULL_NAME" %in% colnames(keeppath_data)) {
  keeppath_data$pathway_name <- keeppath_data$FULL_NAME
}

cat("Loaded", nrow(birewire_data), "pathways from BireWire results\n")
cat("Loaded", nrow(keeppath_data), "pathways from KeepPathSize results\n")

# Ensure empirical_pval is numeric
birewire_data$empirical_pval <- as.numeric(as.character(birewire_data$empirical_pval))
keeppath_data$empirical_pval <- as.numeric(as.character(keeppath_data$empirical_pval))

# 6. Start correlation analysis
cat("\n======= Performing Tissue Specificity Correlation Analysis =======\n")
rank_correlation_results <- data.frame()

# Define ranking methods for correlation analysis
ranking_methods <- list(
  # Raw rankings with tie-breaking
  list(method_name = "PvalueBeta", 
       data = birewire_data, 
       rank_col = c("p_value", "beta_value"),
       higher_better = c(FALSE, TRUE)),
  
  # Empirical p-value + std effect size
  list(method_name = "BireWire_EmpPvalStdBeta", 
       data = birewire_data, 
       rank_col = c("empirical_pval", "std_effect_size"),
       higher_better = c(FALSE, TRUE)),
  list(method_name = "KeepPathSize_EmpPvalStdBeta", 
       data = keeppath_data, 
       rank_col = c("empirical_pval", "std_effect_size"),
       higher_better = c(FALSE, TRUE))
)

# Get all tissue metrics (each individual tissue)
tissue_metrics <- list()
for(tissue_col in tissue_cols) {
  tissue_name <- gsub("_mean$", "", paste0(tissue_col, "_mean"))
  tissue_metrics[[length(tissue_metrics) + 1]] <- list(
    metric = paste0(tissue_col, "_mean"),
    name = tissue_name,
    higher_better = TRUE
  )
}

cat("Analyzing correlation for", length(tissue_metrics), "individual tissues\n")

# Process each ranking method
for(ranking in ranking_methods) {
  method_name <- ranking$method_name
  data <- ranking$data
  rank_cols <- ranking$rank_col
  higher_better <- ranking$higher_better
  
  cat(paste("\nCalculating tissue correlations for", method_name, "...\n"))
  
  # Create a clean copy of data for ranking
  ranking_data <- data
  
  # Create properly oriented ranking columns
  for(i in 1:length(rank_cols)) {
    col <- rank_cols[i]
    # If higher is better, multiply by -1 so lower values = better ranking
    if(higher_better[i]) {
      ranking_data[[paste0("ranking_", i)]] <- -ranking_data[[col]]
    } else {
      ranking_data[[paste0("ranking_", i)]] <- ranking_data[[col]]
    }
  }
  
  # Simple ordering approach
  if(length(rank_cols) == 1) {
    # If just one column, order directly by it
    ranking_data <- ranking_data[order(ranking_data$ranking_1),]
  } else if(length(rank_cols) == 2) {
    # If two columns, order by first column, then second
    ranking_data <- ranking_data[order(ranking_data$ranking_1, ranking_data$ranking_2),]
  }
  
  # Add rank as the row number after ordering
  ranking_data$pathway_rank <- 1:nrow(ranking_data)
  
  # Extract just the pathway name and rank
  ranked_paths <- ranking_data %>%
    select(pathway_name, pathway_rank) %>%
    rename(name = pathway_name)
  
  # Merge with tissue scores - this gives us all pathways with scores
  all_merged_ranks <- merge(ranked_paths, pathway_tissue_scores, by="name")
  
  # Create datasets with the top N pathways (if available)
  top50_ranked_paths <- ranked_paths %>% filter(pathway_rank <= 50)
  top100_ranked_paths <- ranked_paths %>% filter(pathway_rank <= 100)
  top250_ranked_paths <- ranked_paths %>% filter(pathway_rank <= 250)
  top500_ranked_paths <- ranked_paths %>% filter(pathway_rank <= 500)
  
  # Merge the top N pathways with tissue scores
  top50_merged_ranks <- merge(top50_ranked_paths, pathway_tissue_scores, by="name")
  top100_merged_ranks <- merge(top100_ranked_paths, pathway_tissue_scores, by="name")
  top250_merged_ranks <- merge(top250_ranked_paths, pathway_tissue_scores, by="name")
  top500_merged_ranks <- merge(top500_ranked_paths, pathway_tissue_scores, by="name")
  
  # Calculate correlations for all metrics and subsets
  for(subset_name in c("All Pathways", "Top 500 Pathways", "Top 250 Pathways", "Top 100 Pathways", "Top 50 Pathways")) {
    # Select the appropriate dataset
    if(subset_name == "All Pathways") {
      merged_ranks <- all_merged_ranks
    } else if(subset_name == "Top 500 Pathways") {
      merged_ranks <- top500_merged_ranks
    } else if(subset_name == "Top 250 Pathways") {
      merged_ranks <- top250_merged_ranks
    } else if(subset_name == "Top 100 Pathways") {
      merged_ranks <- top100_merged_ranks
    } else {
      merged_ranks <- top50_merged_ranks
    }
    
    # Skip if no data
    if(nrow(merged_ranks) == 0) {
      cat("  No matching pathways found for", method_name, "- skipping", subset_name, "\n")
      next
    }
    
    cat(paste("  Found", nrow(merged_ranks), "pathways for", subset_name, "\n"))
    
    # Calculate correlations for each tissue metric without generating plots
    for(tissue_metric in tissue_metrics) {
      metric <- tissue_metric$metric
      metric_name <- tissue_metric$name
      
      # Skip if metric not available
      if(!metric %in% colnames(merged_ranks)) {
        cat("  Metric", metric, "not found in data - skipping\n")
        next
      }
      
      # Calculate Spearman correlation only
      spearman_cor <- suppressWarnings(cor.test(merged_ranks$pathway_rank,
                                                merged_ranks[[metric]],
                                                method = "spearman"))

      # Add to results dataframe
      method_result <- data.frame(
        trait = trait,
        tool_base = tool_base,
        method = method_name,
        subset = subset_name,
        tissue_metric = metric_name,
        n_pathways = nrow(merged_ranks),
        spearman_rho = unname(spearman_cor$estimate),
        correlation_pvalue = spearman_cor$p.value
      )
      
      rank_correlation_results <- rbind(rank_correlation_results, method_result)
    }
  }
}

# Write correlation results to CSV
if(nrow(rank_correlation_results) > 0) {
  write_results_csv(rank_correlation_results, trait, tool_base, "_tissue_correlation_summary")
}

# Create a table of best performing methods per tissue
best_methods_by_tissue <- rank_correlation_results %>%
  group_by(tissue_metric, subset) %>%
  slice_min(correlation_pvalue, n = 1) %>%
  select(tissue_metric, subset, method, spearman_rho, correlation_pvalue) %>%
  arrange(tissue_metric, subset)

write_results_csv(best_methods_by_tissue, trait, tool_base, "_best_method_by_tissue")

# Create a summary selecting the strongest (lowest p-value) tissue per method and subset
best_tissue_overall_by_method_subset <- rank_correlation_results %>%
  group_by(method, subset) %>%
  slice_min(correlation_pvalue, n = 1) %>%
  select(method, subset, tissue_metric, spearman_rho, correlation_pvalue, n_pathways) %>%
  arrange(subset, method)

write_results_csv(best_tissue_overall_by_method_subset, trait, tool_base, "_best_tissue_overall_by_method_subset")

cat("\nAnalysis complete. Summary files written to current directory.\n")
cat("======= Tissue Specificity Correlation Analysis Complete =======\n")



