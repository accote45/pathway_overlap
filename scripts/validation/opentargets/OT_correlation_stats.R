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

trait <- args[1]
tool_base <- args[2]
birewire_results_file <- args[3]
keeppathsize_results_file <- args[4]
gmt_file <- args[5]

# Optional: comma-separated list of Top-N cutoffs (default: 100,250,500)
# Example: Rscript OT_correlation_stats.R ... "100,250,500"
top_ns <- c(100, 250, 500)

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
birewire_data <- read.table(birewire_results_file, header = TRUE,fill=TRUE)
# Remove "Base" pathway if present (not a real pathway)
birewire_data <- birewire_data[!(birewire_data$pathway_name=="Base"),]

keeppath_data <- read.table(keeppathsize_results_file, header = TRUE,fill=TRUE)
# Remove "Base" pathway if present (not a real pathway)
keeppath_data <- keeppath_data[!(keeppath_data$pathway_name=="Base"),]

cat("Loaded", nrow(birewire_data), "pathways from BireWire results\n")
cat("Loaded", nrow(keeppath_data), "pathways from KeepPathSize results\n")

# Ensure empirical_pval is numeric
birewire_data$empirical_pval <- as.numeric(as.character(birewire_data$empirical_pval))
keeppath_data$empirical_pval <- as.numeric(as.character(keeppath_data$empirical_pval))

# Rank Correlation Analysis
cat("\n======= Performing Rank Correlation Analysis =======\n")
rank_correlation_results <- data.frame()

# Get all unique pathways that have OpenTargets scores
all_paths_with_scores <- pathway_scores %>% 
  filter(!is.na(mean_score)) %>%
  select(name, mean_score, evidence_density)

# Define ranking methods for correlation analysis
ranking_methods <- list(
  # Raw rankings (p-value + beta tie-break)
  list(method_name = "PvalueBeta", 
       data = birewire_data, 
       rank_col = c("p_value", "beta_value"),
       sig_col = "p_value",
       higher_better = c(FALSE, TRUE)),
  
  # REMOVED: RawP method (p-value only ranking)
  
  # Empirical p-value + std effect size (both randomization methods)
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

  # Build subset list dynamically: "All Pathways" + requested Top-N sets
  subset_list <- list("All Pathways" = all_merged_ranks)
  for (N in sort(unique(top_ns))) {
    topN_ranked_paths <- ranked_paths %>% filter(pathway_rank <= N)
    subset_list[[sprintf("Top %d Pathways", N)]] <- merge(topN_ranked_paths, all_paths_with_scores, by = "name")
  }

  # Helper to compute Spearman correlation for mean_score only
  add_corr_row <- function(df, subset_label) {
    if (nrow(df) <= 1) return(NULL)
    sp_mean <- suppressWarnings(cor.test(df$pathway_rank, df$mean_score, method = "spearman"))
    data.frame(
      trait = trait,
      tool_base = tool_base,
      method = method_name,
      subset = subset_label,
      n_pathways = nrow(df),
      spearman_rho = -unname(sp_mean$estimate),  # Flip sign: positive = concordance
      spearman_pvalue = sp_mean$p.value,
      stringsAsFactors = FALSE
    )
  }

  # Compute correlations for all subsets and append to results
  for (lbl in names(subset_list)) {
    df_sub <- subset_list[[lbl]]
    if (nrow(df_sub) > 0) {
      cat(sprintf("  Found %d pathways for %s with OpenTargets scores\n", nrow(df_sub), lbl))
      row <- add_corr_row(df_sub, lbl)
      if (!is.null(row)) rank_correlation_results <- rbind(rank_correlation_results, row)
    } else {
      cat("  No matching pathways found for", method_name, "- skipping", lbl, "correlation analysis\n")
    }
  }
}

# Write correlation results to CSV (now includes All + all requested Top-N subsets)
if(nrow(rank_correlation_results) > 0) {
  write_results_csv(rank_correlation_results, trait, tool_base, "_rank_correlation_summary")
}


cat("======= Rank Correlation Analysis Complete =======\n")











