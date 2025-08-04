library(plyr)
library(tidyverse)
library(data.table)
library(jsonlite)
library(GSA)
library(MatchIt)

# ===== Helper Functions =====

# Normalize column names in a data frame
normalize_column <- function(data, target_col, possible_cols) {
  if(!(target_col %in% colnames(data))) {
    for(col in possible_cols) {
      if(col %in% colnames(data)) {
        cat("Using", col, "as", target_col, "column\n")
        data[[target_col]] <- data[[col]]
        return(data)
      }
    }
    cat("Warning: Could not find column for", target_col, "\n")
  }
  return(data)
}

# Get target column name based on method name
get_target_col <- function(method_name) {
  switch(method_name,
         "BireWire" = "in_birewire_top",
         "KeepPathSize" = "in_keeppath_top",
         "RawP" = "in_raw_p_top",
         "SigBeta" = "in_sig_beta_top",
         "PvalueBeta" = "in_p_beta_top",  # New method
         "EmpPvalEffect" = "in_emp_effect_top",  # New method
         # Default fallback
         paste0("in_", tolower(method_name), "_top"))
}

# Get file suffix based on ranking type
get_file_suffix <- function(ranking_type) {
  switch(ranking_type,
         "empirical_p" = "",
         "raw_p" = "_raw_p",
         "sig_beta" = "_sig_beta",
         "p_beta" = "_p_beta",  # New type
         "emp_effect" = "_emp_effect",  # New type
         # Default
         paste0("_", ranking_type))
}

# Format method name for file output
format_method_name <- function(method_name, ranking_type) {
  if(ranking_type == "empirical_p") {
    return(tolower(gsub("PathSize", "pathsize", method_name)))
  } else {
    return(tolower(method_name))
  }
}

# Calculate comparison metrics between two groups
calculate_metrics <- function(target_group, control_group, metrics) {
  comparison <- data.frame()
  
  for(metric in metrics) {
    if(metric %in% colnames(target_group) && metric %in% colnames(control_group)) {
      target_value <- mean(target_group[[metric]], na.rm=TRUE)
      control_value <- mean(control_group[[metric]], na.rm=TRUE)
      
      # Run t-test and extract statistics
      t_test_result <- t.test(target_group[[metric]], control_group[[metric]])
      
      comparison <- rbind(comparison, data.frame(
        metric = metric,
        target_value = target_value,
        control_value = control_value,
        difference = target_value - control_value,
        percent_diff = 100 * (target_value - control_value) / ifelse(control_value == 0, 1, control_value),
        t_statistic = t_test_result$statistic,
        df = t_test_result$parameter,
        p_value = t_test_result$p.value,
        significant = t_test_result$p.value < 0.05
      ))
    }
  }
  
  return(comparison)
}

# Write results to CSV with consistent naming
write_results_csv <- function(data, trait, tool_base, prefix="", suffix="", verbose=TRUE) {
  filename <- paste0(trait, "_", tool_base, prefix, suffix, ".csv")
  write.csv(data, filename, row.names=FALSE)
  if(verbose) cat("Wrote", nrow(data), "rows to", filename, "\n")
  return(filename)
}

# ===== Main Script =====

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

# 6. Run analysis for each N value
all_results <- list()

for(n in n_values) {
  cat(paste("\nAnalyzing top", n, "pathways...\n"))
  
  # Ensure n is numeric
  n <- as.numeric(n)
  
  # First, identify all unique pathways from both methods
  all_pathways <- unique(c(birewire_data$pathway_name, keeppath_data$pathway_name))
  
  # Create a dataset with method assignment
  matching_data <- data.frame(
    name = all_pathways,
    in_birewire_empP_top = all_pathways %in% head(birewire_data[order(birewire_data$empirical_pval),]$pathway_name, n),
    in_keeppath_empP_top = all_pathways %in% head(keeppath_data[order(keeppath_data$empirical_pval),]$pathway_name, n)
  )
  
  # Add rankings for raw P and beta
  matching_data$in_raw_p_top <- all_pathways %in% head(birewire_data[order(birewire_data$p_value),]$pathway_name, n)
  
	# Add ranking by p_value (lowest) and then beta_value (highest positive)
  matching_data$in_p_beta_top <- all_pathways %in% head(birewire_data[order(birewire_data$p_value, -birewire_data$beta_value), ]$pathway_name,n)
  
  # Add ranking by empirical_pval (lowest) and then std_effect_size (highest positive)
  matching_data$in_birewire_empP_stdbeta_top <- all_pathways %in% head(
    birewire_data[order(birewire_data$empirical_pval, -birewire_data$std_effect_size), ]$pathway_name,n)
  matching_data$in_keeppath_empP_stdbeta_top <- all_pathways %in% head(
    keeppath_data[order(keeppath_data$empirical_pval, -keeppath_data$std_effect_size), ]$pathway_name,n)

  # Add pathway size information - More efficiently with joins
  sizes_bw <- birewire_data %>% select(pathway_name, ngenes) %>% rename(name=pathway_name, size=ngenes)
  sizes_kp <- keeppath_data %>% select(pathway_name, ngenes) %>% rename(name=pathway_name, size=ngenes)
  sizes <- rbind(sizes_bw, sizes_kp) %>% distinct(name, .keep_all=TRUE)
  matching_data <- left_join(matching_data, sizes, by="name")
  
  # Remove pathways with missing size
  matching_data <- matching_data[!is.na(matching_data$size),]
  
  # Add OpenTargets evidence information
  matching_data <- merge(matching_data, pathway_scores, by="name", all.x=TRUE)
  
  # Define ranking methods consistently
  ranking_methods <- list(
    list(method_name = "BireWire_empP", target_col = "in_birewire_empP_top", ranking_type = "empirical_p"),
    list(method_name = "KeepPathSize_empP", target_col = "in_keeppath_empP_top", ranking_type = "empirical_p"),
    list(method_name = "RawP", target_col = "in_raw_p_top", ranking_type = "raw_p"),
    list(method_name = "PvalueBeta", target_col = "in_p_beta_top", ranking_type = "p_beta"),
    list(method_name = "BireWire_EmpPvalStdBeta", target_col = "in_birewire_empP_stdbeta_top", ranking_type = "p_beta"),
    list(method_name = "KeepPathSize_EmpPvalStdBeta", target_col = "in_keeppath_empP_stdbeta_top", ranking_type = "p_beta")
  )
  
  results_list <- list()
  
  # Process each ranking method
  for(ranking in ranking_methods) {
    method_name <- ranking$method_name
    target_col <- ranking$target_col
    ranking_type <- ranking$ranking_type
    
    cat(paste("\nProcessing", method_name, "ranking...\n"))
    
    # Count how many pathways match this ranking
    top_count <- sum(matching_data[[target_col]], na.rm = TRUE)
    if(top_count < n/2) {
      cat("Warning: Only", top_count, "pathways found for", method_name, "ranking. Expected", n, "\n")
      if(top_count < 5) {
        cat("Skipping", method_name, "ranking due to insufficient pathways\n")
        next
      }
    }
    
    # Perform matching
    tryCatch({
      formula <- as.formula(paste(target_col, "~ size"))
      m.out <- matchit(formula, data = matching_data, method = "nearest", ratio = 1)
      matched_data <- match.data(m.out)
      
      if(is.null(matched_data)) {
        cat("Matching failed for", method_name, "top", n, "pathways\n")
        next
      }
      
      # Add metadata
      matched_data$N <- n
      matched_data$method <- method_name
      matched_data$ranking_type <- ranking_type
      
      # Extract groups and calculate metrics
      target_group <- matched_data[matched_data[[target_col]] == TRUE, ]
      control_group <- matched_data[matched_data[[target_col]] == FALSE, ]
      
      metrics <- c("mean_score", "evidence_density", "max_score", "median_score")
      comparison <- calculate_metrics(target_group, control_group, metrics)
      
      # Size verification
      target_size_avg <- mean(target_group$size)
      control_size_avg <- mean(control_group$size)
      size_diff_pct <- 100 * abs(target_size_avg - control_size_avg) / control_size_avg
      
      # Display results
      cat("\nSize-matched comparison (", method_name, " top", n, "vs. size-matched pathways):\n")
      cat("- Average size in", method_name, "top:", round(target_size_avg, 1), 
          "vs. matched pathways:", round(control_size_avg, 1), 
          "(", round(size_diff_pct, 1), "% diff)\n")
      
      for(i in 1:nrow(comparison)) {
        metric_name <- comparison$metric[i]
        target_val <- comparison$target_value[i]
        control_val <- comparison$control_value[i]
        diff <- comparison$difference[i]
        p_val <- comparison$p_value[i]
        
        cat("- Average", metric_name, "in", method_name, "top", n, ":", round(target_val, 4), 
            "vs. matched pathways:", round(control_val, 4), 
            "(diff:", round(diff, 4), ",", 
            ifelse(diff > 0, method_name, "Control"), "better", 
            ", p =", format.pval(p_val, digits=3), ")\n")
      }
      
      # Save matched data files - use helper functions for consistent naming
      file_suffix <- get_file_suffix(ranking_type)
      method_suffix <- format_method_name(method_name, ranking_type)
      
      output_file <- write_results_csv(
        matched_data, 
        trait, 
        tool_base, 
        paste0("_n", n, "_", method_suffix), 
        file_suffix, 
        verbose = FALSE
      )
      cat("Saved matched data to", output_file, "\n")
      
      # Store results
      results_list[[method_name]] <- list(
        matched_data = matched_data,
        metrics = comparison,
        size_balance = list(
          target_size = target_size_avg,
          control_size = control_size_avg,
          difference_pct = size_diff_pct
        )
      )
    }, error = function(e) {
      cat("Error in matching:", e$message, "\n")
      next
    })
  }
  
  # Store results for this N value
  all_results[[paste0("n_", n)]] <- results_list
}

# 7. Calculate advantage metrics
advantage_data <- data.frame()
method_labels <- list(
  "BireWire" = "Fix gene frequency and pathway size",
  "KeepPathSize" = "Fix pathway size",
  "RawP" = "Standard p-value",
  "SigBeta" = "Significant beta"
)

for(n in n_values) {
  result_name <- paste0("n_", n)
  
  if(!is.null(all_results[[result_name]])) {
    n_results <- all_results[[result_name]]
    
    # Process each method's results
    for(method_name in names(n_results)) {
      if(!is.null(n_results[[method_name]])) {
        method_data <- n_results[[method_name]]
        matched_data <- method_data$matched_data
        metrics <- method_data$metrics
        
        # Get the target column consistently using helper function
        target_col <- get_target_col(method_name)
        cat("Using target column", target_col, "for method", method_name, "\n")
        
        # Extract method-specific data
        top_data <- matched_data[matched_data[[target_col]] == TRUE, ]
        control_data <- matched_data[matched_data[[target_col]] == FALSE, ]
        
        # Create advantage data
        method_advantage <- data.frame(
          N = factor(n, levels=n_values),
          method = method_name,
          ranking_type = matched_data$ranking_type[1],
          mean_score_top = mean(top_data$mean_score, na.rm=TRUE),
          mean_score_control = mean(control_data$mean_score, na.rm=TRUE),
          mean_score_advantage = mean(top_data$mean_score, na.rm=TRUE) - mean(control_data$mean_score, na.rm=TRUE),
          evidence_density_top = mean(top_data$evidence_density, na.rm=TRUE),
          evidence_density_control = mean(control_data$evidence_density, na.rm=TRUE),
          evidence_density_advantage = mean(top_data$evidence_density, na.rm=TRUE) - mean(control_data$evidence_density, na.rm=TRUE)
        )
        
        # Add p-values and statistics
        for(metric_name in c("mean_score", "evidence_density")) {
          metric_idx <- which(metrics$metric == metric_name)
          if(length(metric_idx) > 0) {
            method_advantage[[paste0("p_value_", metric_name)]] <- metrics$p_value[metric_idx]
            method_advantage[[paste0("t_statistic_", metric_name)]] <- metrics$t_statistic[metric_idx]
            method_advantage[[paste0("df_", metric_name)]] <- metrics$df[metric_idx]
          }
        }
        
        advantage_data <- rbind(advantage_data, method_advantage)
      }
    }
  }
}

# Add significance indicators
advantage_data$mean_score_sig <- ifelse(advantage_data$p_value_mean_score < 0.05, "*", "")
advantage_data$evidence_density_sig <- ifelse(advantage_data$p_value_evidence_density < 0.05, "*", "")

# 8. Write detailed advantage data
write_results_csv(advantage_data, trait, tool_base, "_detailed_advantage")

# 9. Create summary tables by ranking type
for(ranking_type in unique(advantage_data$ranking_type)) {
  # Filter data for this ranking type
  type_data <- advantage_data[advantage_data$ranking_type == ranking_type,]
  file_suffix <- get_file_suffix(ranking_type)
  
  # Check if we have enough data for this ranking type
  if(nrow(type_data) > 0) {
    tryCatch({
      # Create and save summary for this ranking type
      summary_data <- type_data %>%
        select(N, method, 
               mean_score_advantage, t_statistic_mean_score, df_mean_score, mean_score_sig, 
               evidence_density_advantage, t_statistic_evidence_density, df_evidence_density, evidence_density_sig)
      
      # Check if we have multiple methods before using pivot_wider
      unique_methods <- unique(summary_data$method)
      
      if(length(unique_methods) > 1) {
        # We can use pivot_wider when we have multiple methods
        summary_data <- summary_data %>%
          pivot_wider(
            names_from = method,
            values_from = c(mean_score_advantage, t_statistic_mean_score, df_mean_score, mean_score_sig, 
                           evidence_density_advantage, t_statistic_evidence_density, df_evidence_density, evidence_density_sig)
          )
      }
      
      # Write summary data to CSV
      write_results_csv(summary_data, trait, tool_base, "_advantage_summary", file_suffix)
      
    }, error = function(e) {
      cat("Error creating summary for ranking type", ranking_type, ":", e$message, "\n")
      
      # Create a simple summary without pivot_wider
      simple_summary <- type_data %>%
        select(N, method, ranking_type, 
               mean_score_advantage, mean_score_sig, 
               evidence_density_advantage, evidence_density_sig)
      
      write_results_csv(simple_summary, trait, tool_base, "_advantage_summary", paste0(file_suffix, "_simple"))
    })
  } else {
    cat("No data for ranking type:", ranking_type, "- skipping summary file\n")
  }
}

cat("======= Analysis Complete =======\n")