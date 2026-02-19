library(plyr)
  library(tidyverse)
  library(data.table)
  library(GSA)

write_results_csv <- function(data, trait, tool_base, prefix="", suffix="", verbose=TRUE) {
  filename <- paste0(trait, "_", tool_base, prefix, ".csv")
  write.csv(data, filename, row.names=FALSE)
  if (verbose) cat("Wrote", nrow(data), "rows to", filename, "\n")
  return(filename)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript dorothea_correlation.R <trait> <tool_base> <dorothea_path> <birewire_results> <keeppathsize_results>")
}

trait <- args[1]
tool_base <- args[2]
dorothea_path <- args[3]
birewire_results_file <- args[4]
keeppathsize_results_file <- args[5]

# Top-N cutoffs
top_ns <- c(100, 250, 500)

cat("======= Starting Dorothea correlation for", trait, "(", tool_base, ") =======\n")


# Read dorothea csv
pathway_scores <- read.csv(dorothea_path)

# 4) Load empirical results files (BireWire, KeepPathSize)
birewire_data <- read.table(birewire_results_file, header = TRUE,fill=TRUE)
# Remove "Base" pathway if present (not a real pathway)
birewire_data <- birewire_data[!(birewire_data$pathway_name=="Base"),]
keeppath_data <- read.table(keeppathsize_results_file, header = TRUE,fill=TRUE)
# Remove "Base" pathway if present (not a real pathway)
keeppath_data <- keeppath_data[!(keeppath_data$pathway_name=="Base"),]

cat("Loaded", nrow(birewire_data), "pathways from BireWire results\n")
cat("Loaded", nrow(keeppath_data), "pathways from KeepPathSize results\n")

# Ensure numeric fields
num_cols <- c("empirical_pval", "p_value", "beta_value", "std_effect_size")
for (nm in num_cols) {
  if (nm %in% names(birewire_data)) birewire_data[[nm]] <- as.numeric(as.character(birewire_data[[nm]]))
  if (nm %in% names(keeppath_data)) keeppath_data[[nm]] <- as.numeric(as.character(keeppath_data[[nm]]))
}

# 5) Rank Correlation Analysis (style kept close to OT script)
cat("\n======= Performing Rank Correlation Analysis (MalaCards) =======\n")
rank_correlation_results <- data.frame()

all_paths_with_scores <- pathway_scores %>% 
  filter(!is.na(score))
all_paths_with_scores$name <- all_paths_with_scores$pathway

ranking_methods <- list(
  # Raw rankings (p-value + beta tie-break)
  list(method_name = "PvalueBeta", 
       data = birewire_data, 
       rank_col = c("p_value", "beta_value"),
       sig_col = "p_value",
       higher_better = c(FALSE, TRUE)),
  
  # REMOVED: RawP method (p-value only ranking)
  
  # Empirical p-value + std effect size
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

# 5) Rank Correlation Analysis (style kept close to OT script)
cat("\n======= Performing Rank Correlation Analysis (MalaCards) =======\n")
rank_correlation_results <- data.frame()

all_paths_with_scores <- pathway_scores %>% 
  filter(!is.na(score))
all_paths_with_scores$name <- all_paths_with_scores$pathway

for (ranking in ranking_methods) {
  method_name <- ranking$method_name
  data <- ranking$data
  rank_cols <- ranking$rank_col
  higher_better <- ranking$higher_better
  
  cat(paste("\nCalculating rank correlation for", method_name, "...\n"))
  
  ranking_data <- data
  
  # Build rank keys (invert if higher is better so smaller is better)
  for (i in seq_along(rank_cols)) {
    col <- rank_cols[i]
    if (!col %in% names(ranking_data)) next
    ranking_data[[paste0("ranking_", i)]] <- if (isTRUE(higher_better[i])) -ranking_data[[col]] else ranking_data[[col]]
  }
  
  # Order by 1 or 2 keys (as in OT script)
  if (length(rank_cols) == 1 && "ranking_1" %in% names(ranking_data)) {
    ranking_data <- ranking_data[order(ranking_data$ranking_1), ]
  } else if (length(rank_cols) >= 2 && all(c("ranking_1", "ranking_2") %in% names(ranking_data))) {
    ranking_data <- ranking_data[order(ranking_data$ranking_1, ranking_data$ranking_2), ]
  } else {
    # If required columns missing, skip
    cat("  Missing columns for", method_name, "- skipping.\n")
    next
  }
  
  ranking_data$pathway_rank <- seq_len(nrow(ranking_data))
  ranked_paths <- ranking_data %>% select(pathway_name, pathway_rank) %>% rename(name = pathway_name)
  
  # Merge with MalaCards pathway scores
  all_merged_ranks <- merge(ranked_paths, all_paths_with_scores, by = "name")
  
  # Subsets: All + Top-N
  subset_list <- list("All Pathways" = all_merged_ranks)
  for (N in sort(unique(top_ns))) {
    topN_ranked_paths <- ranked_paths %>% filter(pathway_rank <= N)
    subset_list[[sprintf("Top %d Pathways", N)]] <- merge(topN_ranked_paths, all_paths_with_scores, by = "name")
  }
  
  add_corr_row <- function(df, subset_label) {
    if (nrow(df) <= 1) return(NULL)
    sp_mean <- suppressWarnings(cor.test(df$pathway_rank, df$score, method = "spearman"))
    data.frame(
      trait = trait,
      tool_base = tool_base,
      method = method_name,
      subset = subset_label,
      n_pathways = nrow(df),
      rank_mean_score_correlation = unname(sp_mean$estimate),
      rank_mean_score_pvalue = sp_mean$p.value,
      stringsAsFactors = FALSE
    )
  }
  
  for (lbl in names(subset_list)) {
    df_sub <- subset_list[[lbl]]
    if (nrow(df_sub) > 0) {
      cat(sprintf("  Found %d pathways for %s with Dorothea scores\n", nrow(df_sub), lbl))
      row <- add_corr_row(df_sub, lbl)
      if (!is.null(row)) rank_correlation_results <- rbind(rank_correlation_results, row)
    } else {
      cat("  No matching pathways found for", method_name, "- skipping", lbl, "\n")
    }
  }
}

# 6) Write correlation results
if (nrow(rank_correlation_results) > 0) {
  write_results_csv(rank_correlation_results, trait, tool_base, "_dorothea_rank_correlation_summary")
  cat("======= Rank Correlation Analysis Complete =======\n")
} else {
  cat("No results to write.\n")
}


