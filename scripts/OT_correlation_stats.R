library(plyr)
library(tidyverse)
library(data.table)
library(jsonlite)
library(GSA)

write_results_csv <- function(data, trait, tool_base, prefix="", suffix="", verbose=TRUE) {
  filename <- paste0(trait, "_", tool_base, prefix, ".csv")
  
  write.csv(data, filename, row.names=FALSE)
  if(verbose) cat("Wrote", nrow(data), "rows to", filename, "\n")
  return(filename)
}

# Process command line arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 5) {
  stop("Usage: Rscript size_matched_OT_stats_optimized.R <trait> <tool_base> <birewire_results_file> <keeppathsize_results_file> <gmt_file> [top_n_values]")
}

trait <- args[1]
tool_base <- args[2]
birewire_results_file <- args[3]
keeppathsize_results_file <- args[4]
gmt_file <- args[5]


# Optional: parse comma-separated n_values or use default
n_values <- c(10, 20, 50, 100)  # Default
if(length(args) >= 6) {
  n_values <- as.numeric(unlist(strsplit(args[6], ",")))
}

# Define trait mapping
trait_mapping <- list(
  "t2d" = "MONDO_0005148",
  "cad" = "EFO_0001645", 
  "ad" = "MONDO_0004975",
  "mdd" = "MONDO_0002050",
  "scz" = "MONDO_0005090",
  "ibd" = "EFO_0000555",
  "breast" = "MONDO_0007254",
  "HDL_cholesterol" = "EFO_0004612",
  "Lymphocyte_count" = "EFO_0004587",
  "Platelet_crit" = "EFO_0007985",
  "Alkaline_phosphatase" = "EFO_0004533"
)

# Get disease ID from mapping
disease_id <- trait_mapping[[trait]]
if(is.null(disease_id)) {
  stop(paste("No mapping found for trait:", trait))
}

cat("======= Starting Size-Matched Analysis for", trait, "=======\n")
cat("Using disease ID:", disease_id, "\n")

# 1. Load pathway data
cat("Loading pathway data from", gmt_file, "...\n")
dat <- GSA.read.gmt(gmt_file)
path.list <- dat$genesets
names(path.list) <- dat$geneset.names

# Convert the path list to a data table
genes_long <- rbindlist(lapply(names(path.list), function(name) {
  data.table(value = path.list[[name]], name = name)
}))
genes_long$targetId <- genes_long$value

# 2. Process OpenTargets JSON files for disease
# Fix redundant file loading
files <- list.files(path = "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/results/drugtarget_test/associationByDatatypeDirect", 
                    pattern = "*.json", full.names=TRUE)
if(length(files) == 0) {
  stop("No JSON files found in current directory")
}

cat("Reading JSON files for disease", disease_id, "...\n")
dat <- lapply(files, function(file) {
  tryCatch({
    json_data <- stream_in(file(file), verbose = FALSE)
    json_data <- as.data.frame(json_data)
    json_data %>%
      filter(diseaseId == disease_id & datatypeId %in% c('animal_model', 'known_drug', 'literature'))
  }, error = function(e) {
    cat("Error reading file:", file, "\n")
    return(NULL)
  })
})

disease_data <- do.call(base::rbind, dat[!sapply(dat, is.null)])
if(nrow(disease_data) == 0) {
  stop(paste("No OpenTargets data found for disease ID:", disease_id))
}

# Filter and format the gene-disease associations
gene_disease_associations <- disease_data %>%
  select(targetId, diseaseId, datatypeId, score, evidenceCount) %>%
  arrange(desc(evidenceCount))

# Rank Correlation Analysis
cat("\n======= Performing Rank Correlation Analysis =======\n")
rank_correlation_results <- data.frame()

# Get all unique pathways that have OpenTargets scores
all_paths_with_scores <- pathway_scores %>% 
  filter(!is.na(mean_score)) %>%
  select(name, mean_score, evidence_density)


# Save to CSV
write_results_csv(gene_disease_associations, trait, tool_base, "_gene_disease_associations")

# Select max evidence count for each target - directly from disease data
disease_targets <- disease_data %>%
  group_by(targetId) %>%
  slice_max(evidenceCount, with_ties = FALSE) %>%
  ungroup()

# 3. Calculate pathway scores
cat("Calculating pathway-level scores...\n")
masterfin <- merge(genes_long, disease_targets, by="targetId", all.x=TRUE)
masterfin[is.na(masterfin)] <- 0

pathway_scores <- masterfin %>% 
  group_by(name) %>% 
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    num_genes = n(),
    num_with_evidence = sum(score > 0, na.rm = TRUE),
    evidence_density = num_with_evidence / num_genes,
    max_score = max(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE)
  ) %>% 
  as.data.frame()

# 4. Load empirical results files
birewire_data <- read.table(birewire_results_file, header=T)
keeppath_data <- read.table(keeppathsize_results_file, header=T)

cat("Loaded", nrow(birewire_data), "pathways from BireWire results\n")
cat("Loaded", nrow(keeppath_data), "pathways from KeepPathSize results\n")

# Ensure empirical_pval is numeric
birewire_data$empirical_pval <- as.numeric(as.character(birewire_data$empirical_pval))
keeppath_data$empirical_pval <- as.numeric(as.character(keeppath_data$empirical_pval))

# Define ranking methods for correlation analysis
ranking_methods <- list(
  list(method_name = "BireWire_empP", 
       data = birewire_data, 
       rank_col = "empirical_pval",
       sig_col = "empirical_pval",
       higher_better = FALSE),
  list(method_name = "KeepPathSize_empP", 
       data = keeppath_data, 
       rank_col = "empirical_pval",
       sig_col = "empirical_pval",
       higher_better = FALSE),
  list(method_name = "RawP", 
       data = birewire_data, 
       rank_col = "p_value",
       sig_col = "p_value",
       higher_better = FALSE),
  list(method_name = "PvalueBeta", 
       data = birewire_data, 
       rank_col = c("p_value", "beta_value"),
       sig_col = "p_value",
       higher_better = c(FALSE, TRUE)),
  list(method_name = "BireWire_EmpPvalStdBeta", 
       data = birewire_data, 
       rank_col = c("empirical_pval", "std_effect_size"),
       sig_col = "empirical_pval",
       higher_better = c(FALSE, TRUE)),
  list(method_name = "KeepPathSize_EmpPvalStdBeta", 
       data = keeppath_data, 
       rank_col = c("empirical_pval", "std_effect_size"),
       sig_col = "empirical_pval",
       higher_better = c(FALSE, TRUE))
)

# Process each ranking method
for(ranking in ranking_methods) {
  method_name <- ranking$method_name
  data <- ranking$data
  rank_cols <- ranking$rank_col
  higher_better <- ranking$higher_better
  
  cat(paste("\nCalculating rank correlation for", method_name, "...\n"))
  
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
  
  # Merge with OpenTargets scores - this gives us all pathways with scores
  all_merged_ranks <- merge(ranked_paths, all_paths_with_scores, by="name")
  
  # Also create a dataset with only the top 500 pathways (if available)
  top500_ranked_paths <- ranked_paths %>%
    filter(pathway_rank <= 500)
  
  # Merge the top 500 with OpenTargets scores
  top500_merged_ranks <- merge(top500_ranked_paths, all_paths_with_scores, by="name")
  
  # Calculate correlations for all pathways
  if(nrow(all_merged_ranks) > 0) {
    cat(paste("  Found", nrow(all_merged_ranks), "pathways with OpenTargets scores\n"))
    
    # Spearman correlation with mean_score
    all_spearman_cor <- cor.test(all_merged_ranks$pathway_rank, 
                                all_merged_ranks$mean_score, 
                                method = "spearman")
    
    # Spearman correlation with evidence_density
    all_spearman_density <- cor.test(all_merged_ranks$pathway_rank, 
                                    all_merged_ranks$evidence_density, 
                                    method = "spearman")
    
    # Add to results
    all_method_result <- data.frame(
      method = method_name,
      subset = "All Pathways",
      n_pathways = nrow(all_merged_ranks),
      rank_mean_score_correlation = all_spearman_cor$estimate,
      rank_mean_score_pvalue = all_spearman_cor$p.value,
      rank_evidence_density_correlation = all_spearman_density$estimate,
      rank_evidence_density_pvalue = all_spearman_density$p.value
    )
    
    rank_correlation_results <- rbind(rank_correlation_results, all_method_result)
    
    # Create scatterplots for all pathways
    pdf(paste0(trait, "_", tool_base, "_", gsub("[[:punct:]]", "_", method_name), "_all_pathways_correlation.pdf"), width=12, height=6)
    
    # Set up a 1x2 plot layout
    par(mfrow=c(1,2))
    
    # Plot 1: All Pathway Ranks vs Mean Score
    plot(all_merged_ranks$pathway_rank, all_merged_ranks$mean_score,
         main=paste(method_name, "All Pathway Ranks vs OpenTargets Score"),
         xlab="Pathway Rank (lower is better)",
         ylab="OpenTargets Mean Score",
         pch=19, col=rgb(0,0,1,0.3))
    abline(lm(mean_score ~ pathway_rank, data=all_merged_ranks), col="red", lwd=2)
    legend("topright", 
           legend=c(paste("Spearman rho =", round(all_spearman_cor$estimate, 3)),
                   paste("p-value =", format.pval(all_spearman_cor$p.value, digits=3)),
                   paste("n =", nrow(all_merged_ranks), "pathways")),
           bty="n")
    
    # Plot 2: All Pathway Ranks vs Evidence Density
    plot(all_merged_ranks$pathway_rank, all_merged_ranks$evidence_density,
         main=paste(method_name, "All Pathway Ranks vs Evidence Density"),
         xlab="Pathway Rank (lower is better)",
         ylab="OpenTargets Evidence Density",
         pch=19, col=rgb(0,0.7,0,0.3))
    abline(lm(evidence_density ~ pathway_rank, data=all_merged_ranks), col="red", lwd=2)
    legend("topright", 
           legend=c(paste("Spearman rho =", round(all_spearman_density$estimate, 3)),
                   paste("p-value =", format.pval(all_spearman_density$p.value, digits=3)),
                   paste("n =", nrow(all_merged_ranks), "pathways")),
           bty="n")
    
    dev.off()
    cat("  Created correlation plots for all pathways using", method_name, "\n")
  } else {
    cat("  No matching pathways found for", method_name, "- skipping correlation analysis\n")
  }
  
  # Calculate correlations for top 500 pathways
  if(nrow(top500_merged_ranks) > 0) {
    cat(paste("  Found", nrow(top500_merged_ranks), "of top 500 pathways with OpenTargets scores\n"))
    
    # Spearman correlation with mean_score
    top500_spearman_cor <- cor.test(top500_merged_ranks$pathway_rank, 
                               top500_merged_ranks$mean_score, 
                               method = "spearman")
    
    # Spearman correlation with evidence_density
    top500_spearman_density <- cor.test(top500_merged_ranks$pathway_rank, 
                                   top500_merged_ranks$evidence_density, 
                                   method = "spearman")
    
    # Add to results
    top500_method_result <- data.frame(
      method = method_name,
      subset = "Top 500 Pathways",
      n_pathways = nrow(top500_merged_ranks),
      rank_mean_score_correlation = top500_spearman_cor$estimate,
      rank_mean_score_pvalue = top500_spearman_cor$p.value,
      rank_evidence_density_correlation = top500_spearman_density$estimate,
      rank_evidence_density_pvalue = top500_spearman_density$p.value
    )
    
    rank_correlation_results <- rbind(rank_correlation_results, top500_method_result)
    
    # Create scatterplots for top 500 pathways
    pdf(paste0(trait, "_", tool_base, "_", gsub("[[:punct:]]", "_", method_name), "_top500_correlation.pdf"), width=12, height=6)
    
    # Set up a 1x2 plot layout
    par(mfrow=c(1,2))
    
    # Plot 1: Top 500 Pathway Ranks vs Mean Score
    plot(top500_merged_ranks$pathway_rank, top500_merged_ranks$mean_score,
         main=paste(method_name, "Top 500 Pathway Ranks vs OpenTargets Score"),
         xlab="Pathway Rank (lower is better)",
         ylab="OpenTargets Mean Score",
         pch=19, col=rgb(0,0,1,0.3))
    abline(lm(mean_score ~ pathway_rank, data=top500_merged_ranks), col="red", lwd=2)
    legend("topright", 
           legend=c(paste("Spearman rho =", round(top500_spearman_cor$estimate, 3)),
                   paste("p-value =", format.pval(top500_spearman_cor$p.value, digits=3)),
                   paste("n =", nrow(top500_merged_ranks), "pathways")),
           bty="n")
    
    # Plot 2: Top 500 Pathway Ranks vs Evidence Density
    plot(top500_merged_ranks$pathway_rank, top500_merged_ranks$evidence_density,
         main=paste(method_name, "Top 500 Pathway Ranks vs Evidence Density"),
         xlab="Pathway Rank (lower is better)",
         ylab="OpenTargets Evidence Density",
         pch=19, col=rgb(0,0.7,0,0.3))
    abline(lm(evidence_density ~ pathway_rank, data=top500_merged_ranks), col="red", lwd=2)
    legend("topright", 
           legend=c(paste("Spearman rho =", round(top500_spearman_density$estimate, 3)),
                   paste("p-value =", format.pval(top500_spearman_density$p.value, digits=3)),
                   paste("n =", nrow(top500_merged_ranks), "pathways")),
           bty="n")
    
    dev.off()
    cat("  Created correlation plots for top 500 pathways using", method_name, "\n")
  } else {
    cat("  No matching top 500 pathways found for", method_name, "- skipping top 500 correlation analysis\n")
  }
}

# Write correlation results to CSV
if(nrow(rank_correlation_results) > 0) {
  write_results_csv(rank_correlation_results, trait, tool_base, "_rank_correlation_summary")
  
  # Print summary of results
  cat("\nRank Correlation Analysis Summary:\n")
  
  # Sort by subset first, then by mean score correlation p-value
  rank_correlation_results <- rank_correlation_results %>%
    arrange(subset, rank_mean_score_pvalue)
  
  # Print all pathways results
  cat("\nALL PATHWAYS:\n")
  all_pathways_results <- rank_correlation_results %>% filter(subset == "All Pathways")
  for(i in 1:nrow(all_pathways_results)) {
    method <- all_pathways_results$method[i]
    cor_score <- all_pathways_results$rank_mean_score_correlation[i]
    pval_score <- all_pathways_results$rank_mean_score_pvalue[i]
    cor_density <- all_pathways_results$rank_evidence_density_correlation[i]
    pval_density <- all_pathways_results$rank_evidence_density_pvalue[i]
    n_paths <- all_pathways_results$n_pathways[i]
    
    cat(paste0(method, " (", n_paths, " pathways): \n"))
    cat(paste0("  Mean Score Correlation: ", round(cor_score, 3), 
              " (p = ", format.pval(pval_score, digits=3), ")\n"))
    cat(paste0("  Evidence Density Correlation: ", round(cor_density, 3), 
              " (p = ", format.pval(pval_density, digits=3), ")\n"))
  }
  
  # Print top 500 pathways results
  cat("\nTOP 500 PATHWAYS:\n")
  top500_results <- rank_correlation_results %>% filter(subset == "Top 500 Pathways")
  for(i in 1:nrow(top500_results)) {
    method <- top500_results$method[i]
    cor_score <- top500_results$rank_mean_score_correlation[i]
    pval_score <- top500_results$rank_mean_score_pvalue[i]
    cor_density <- top500_results$rank_evidence_density_correlation[i]
    pval_density <- top500_results$rank_evidence_density_pvalue[i]
    n_paths <- top500_results$n_pathways[i]
    
    cat(paste0(method, " (", n_paths, " of top 500 pathways): \n"))
    cat(paste0("  Mean Score Correlation: ", round(cor_score, 3), 
              " (p = ", format.pval(pval_score, digits=3), ")\n"))
    cat(paste0("  Evidence Density Correlation: ", round(cor_density, 3), 
              " (p = ", format.pval(pval_density, digits=3), ")\n"))
  }
}

cat("======= Rank Correlation Analysis Complete =======\n")











