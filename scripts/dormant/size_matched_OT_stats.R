# ===== Top Pathway Prioritization Analysis - Statistics =====
# This script performs statistical analysis of top-ranked pathways
# identified by different randomization methods

library(plyr)
library(tidyverse)
library(data.table)
library(jsonlite)
library(GSA)
library(MatchIt)

# Process command line arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 5) {
  stop("Usage: Rscript size_matched_OT_stats.R <trait> <tool_base> <birewire_results_file> <keeppathsize_results_file> <gmt_file> [top_n_values]")
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

# === 1. Data Loading Functions ===

load_pathway_data <- function(gmt_file) {
  cat("Loading pathway data from", gmt_file, "...\n")
  dat <- GSA.read.gmt(gmt_file)
  path.list <- dat$genesets
  names(path.list) <- dat$geneset.names
  
  # Convert the path list to a data table
  genes_long <- rbindlist(lapply(names(path.list), function(name) {
    data.table(value = path.list[[name]], name = name)
  }))
  genes_long$targetId <- genes_long$value
  
  return(genes_long)
}

process_disease_data <- function(files, disease_id) {
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
  
  combined_dat <- do.call(base::rbind, dat[!sapply(dat, is.null)])
  
  if(nrow(combined_dat) == 0) {
    stop(paste("No OpenTargets data found for disease ID:", disease_id))
  }
  return(combined_dat)
}

calculate_pathway_scores <- function(genes_long, disease_targets) {
  cat("Calculating pathway-level scores...\n")
  
  # Merge OpenTargets data with pathway genes
  masterfin <- merge(genes_long, disease_targets, by="targetId", all.x=TRUE)
  masterfin[is.na(masterfin)] <- 0
  
  # Calculate pathway-level scores
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
  
  return(pathway_scores)
}

# === 2. Matching and Analysis Functions ===

prepare_matching_data <- function(birewire_results, keeppath_results, pathway_scores, n) {
  cat(paste("\nPreparing data for top", n, "pathways...\n"))
  
  # Ensure n is numeric
  n <- as.numeric(n)
  
  # Print column names to debug
  cat("BireWire columns:", paste(colnames(birewire_results), collapse=", "), "\n")
  cat("KeepPath columns:", paste(colnames(keeppath_results), collapse=", "), "\n")
  
  # First, identify all unique pathways from both methods
  all_pathways <- unique(c(birewire_results$FULL_NAME, keeppath_results$FULL_NAME))
  
  # Check columns and normalize names - needed columns are:
  # - FULL_NAME: pathway name
  # - empirical_pval: empirical p-value (method-specific)
  # - P: raw p-value (same for both methods)
  # - BETA: effect size (same for both methods)
  # - NGENES: gene count (same for both methods)
  
  # Check and normalize column names for empirical p-values
  prepare_empirical_cols(birewire_results, keeppath_results)
  
  # Check and normalize column names for raw p-values and beta 
  # (only need to check one file since these are identical)
  prepare_raw_stats_cols(birewire_results)
  
  # Ensure empirical_pval is numeric
  birewire_results$empirical_pval <- as.numeric(as.character(birewire_results$empirical_pval))
  keeppath_results$empirical_pval <- as.numeric(as.character(keeppath_results$empirical_pval))
  
  # Create a dataset with method assignment and pathway size
  matching_data <- data.frame(
    name = all_pathways,
    in_birewire_top = all_pathways %in% head(birewire_results[order(birewire_results$empirical_pval),]$FULL_NAME, n),
    in_keeppath_top = all_pathways %in% head(keeppath_results[order(keeppath_results$empirical_pval),]$FULL_NAME, n)
  )
  
  # Add rankings for raw P and beta (these are the same for both methods, so only need one column each)
  matching_data$in_raw_p_top <- all_pathways %in% head(birewire_results[order(birewire_results$P),]$FULL_NAME, n)
  
  # Add rankings for significant beta
  sig_results <- birewire_results[birewire_results$P < 0.05,]
  if(nrow(sig_results) > 0) {
    matching_data$in_sig_beta_top <- all_pathways %in% head(sig_results[order(-sig_results$BETA),]$FULL_NAME, n)
  } else {
    matching_data$in_sig_beta_top <- FALSE
  }
  
  # Add pathway size information
  matching_data$size <- NA
  for(i in 1:nrow(matching_data)) {
    idx <- which(birewire_results$FULL_NAME == matching_data$name[i])
    if(length(idx) > 0) {
      matching_data$size[i] <- birewire_results$NGENES[idx[1]]
    } else {
      idx <- which(keeppath_results$FULL_NAME == matching_data$name[i])
      if(length(idx) > 0) {
        matching_data$size[i] <- keeppath_results$NGENES[idx[1]]
      }
    }
  }
  
  # Remove pathways with missing size
  matching_data <- matching_data[!is.na(matching_data$size),]
  
  # Add OpenTargets evidence information
  matching_data <- merge(matching_data, pathway_scores, by="name", all.x=TRUE)
  
  return(matching_data)
}

# Helper functions to check and normalize column names
prepare_empirical_cols <- function(birewire_results, keeppath_results) {
  # Check for empirical_pval column and normalize if needed
  if(!"empirical_pval" %in% colnames(birewire_results)) {
    possible_columns <- c("EMPIRICAL_P", "P", "P.value", "pvalue", "empirical_pvalue")
    for(col in possible_columns) {
      if(col %in% colnames(birewire_results)) {
        cat("Using", col, "as empirical p-value column for BireWire\n")
        birewire_results$empirical_pval <- birewire_results[[col]]
        break
      }
    }
  }
  
  if(!"empirical_pval" %in% colnames(keeppath_results)) {
    possible_columns <- c("EMPIRICAL_P", "P", "P.value", "pvalue", "empirical_pvalue")
    for(col in possible_columns) {
      if(col %in% colnames(keeppath_results)) {
        cat("Using", col, "as empirical p-value column for KeepPathSize\n")
        keeppath_results$empirical_pval <- keeppath_results[[col]]
        break
      }
    }
  }
  
  # Check for FULL_NAME column and normalize if needed
  if(!"FULL_NAME" %in% colnames(birewire_results)) {
    possible_columns <- c("PATHWAY", "name", "pathway", "SET")
    for(col in possible_columns) {
      if(col %in% colnames(birewire_results)) {
        cat("Using", col, "as pathway name column for BireWire\n")
        birewire_results$FULL_NAME <- birewire_results[[col]]
        break
      }
    }
  }
  
  if(!"FULL_NAME" %in% colnames(keeppath_results)) {
    possible_columns <- c("PATHWAY", "name", "pathway", "SET")
    for(col in possible_columns) {
      if(col %in% colnames(keeppath_results)) {
        cat("Using", col, "as pathway name column for KeepPathSize\n")
        keeppath_results$FULL_NAME <- keeppath_results[[col]]
        break
      }
    }
  }
}

prepare_raw_stats_cols <- function(results) {
  # Check for raw P-value column
  if(!"P" %in% colnames(results)) {
    possible_columns <- c("P", "p", "p.value", "pval", "P.value")
    for(col in possible_columns) {
      if(col %in% colnames(results)) {
        cat("Using", col, "as raw p-value column\n")
        results$P <- results[[col]]
        break
      }
    }
  }
  
  # Check for BETA column
  if(!"BETA" %in% colnames(results)) {
    possible_columns <- c("BETA", "beta", "effect", "EFFECT")
    for(col in possible_columns) {
      if(col %in% colnames(results)) {
        cat("Using", col, "as beta/effect size column\n")
        results$BETA <- results[[col]]
        break
      }
    }
  }
}

perform_matching <- function(matching_data, target_col) {
  # Perform size matching for target column (either in_birewire_top or in_keeppath_top)
  tryCatch({
    formula <- as.formula(paste(target_col, "~ size"))
    m.out <- matchit(formula, data = matching_data, method = "nearest", ratio = 1)
    matched_data <- match.data(m.out)
    return(matched_data)
  }, error = function(e) {
    cat("Error in matching:", e$message, "\n")
    return(NULL)
  })
}

calculate_comparison_metrics <- function(matched_data, target_col, metrics=c("mean_score", "evidence_density", "max_score", "median_score")) {
  # Extract group of interest and control group based on target column
  target_group <- matched_data[matched_data[[target_col]] == TRUE, ]
  control_group <- matched_data[matched_data[[target_col]] == FALSE, ]
  
  # Compare metrics
  comparison <- data.frame()
  
  for(metric in metrics) {
    if(metric %in% colnames(matched_data)) {
      target_value <- mean(target_group[[metric]], na.rm=TRUE)
      control_value <- mean(control_group[[metric]], na.rm=TRUE)
      
      # Run t-test and extract all statistics
      t_test_result <- t.test(target_group[[metric]], control_group[[metric]])
      
      comparison <- rbind(comparison, data.frame(
        metric = metric,
        target_value = target_value,
        control_value = control_value,
        difference = target_value - control_value,
        percent_diff = 100 * (target_value - control_value) / ifelse(control_value == 0, 1, control_value),
        t_statistic = t_test_result$statistic,
        df = t_test_result$parameter,  # degrees of freedom
        p_value = t_test_result$p.value,
        significant = t_test_result$p.value < 0.05
      ))
    }
  }
  
  # Size verification
  target_size_avg <- mean(target_group$size)
  control_size_avg <- mean(control_group$size)
  size_diff_pct <- 100 * abs(target_size_avg - control_size_avg) / control_size_avg
  
  attr(comparison, "size_balance") <- list(
    target_size = target_size_avg,
    control_size = control_size_avg,
    difference_pct = size_diff_pct
  )
  
  return(comparison)
}

# === 3. Analysis Functions ===

analyze_method_at_n <- function(birewire_results, keeppath_results, pathway_scores, n, disease_name) {
  cat(paste("\nAnalyzing top", n, "pathways...\n"))
  
  # 1. Prepare data for matching
  matching_data <- prepare_matching_data(birewire_results, keeppath_results, pathway_scores, n)
  
  # List of ranking methods to analyze (modified to reflect identical raw p-value and beta)
  ranking_methods <- list(
    list(method_name = "BireWire", target_col = "in_birewire_top", ranking_type = "empirical_p"),
    list(method_name = "KeepPathSize", target_col = "in_keeppath_top", ranking_type = "empirical_p"),
    list(method_name = "RawP", target_col = "in_raw_p_top", ranking_type = "raw_p"),
    list(method_name = "SigBeta", target_col = "in_sig_beta_top", ranking_type = "sig_beta")
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
    matched_data <- perform_matching(matching_data, target_col)
    if(is.null(matched_data)) {
      cat("Matching failed for", method_name, "top", n, "pathways\n")
      next
    }
    
    # Add metadata
    matched_data$N <- n
    matched_data$method <- method_name
    matched_data$ranking_type <- ranking_type
    
    # Calculate metrics
    comparison <- calculate_comparison_metrics(matched_data, target_col)
    size_balance <- attr(comparison, "size_balance")
    
    # Display results
    cat("\nSize-matched comparison (", method_name, " top", n, "vs. size-matched pathways):\n")
    cat("- Average size in", method_name, "top:", round(size_balance$target_size, 1), 
        "vs. matched pathways:", round(size_balance$control_size, 1), 
        "(", round(size_balance$difference_pct, 1), "% diff)\n")
    
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
    
    # Save matched data files
    file_suffix <- switch(ranking_type,
                        "empirical_p" = "",
                        "raw_p" = "_raw_p",
                        "sig_beta" = "_sig_beta")
    
    method_suffix <- if(ranking_type == "empirical_p") {
      tolower(gsub("PathSize", "pathsize", method_name))
    } else {
      tolower(method_name)
    }
    
    output_file <- paste0(disease_name, "_", tool_base, "_n", n, "_", 
                         method_suffix, file_suffix, "_matched.csv")
    write.csv(matched_data, file = output_file)
    
    # Store results
    results_list[[method_name]] <- list(
      matched_data = matched_data,
      metrics = comparison,
      size_balance = size_balance
    )
  }
  
  return(results_list)
}

calculate_method_advantages <- function(all_results, n_values, new_labels) {
  advantage_data <- data.frame()
  
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
          
          # Determine target column name
          target_col <- names(matched_data)[grep("^in_", names(matched_data))][1]
          
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
  
  return(advantage_data)
}

# === Main Execution ===
main <- function() {
  # Get disease ID from mapping
  disease_id <- trait_mapping[[trait]]
  if(is.null(disease_id)) {
    stop(paste("No mapping found for trait:", trait))
  }
  
  cat("======= Starting Size-Matched Analysis for", trait, "=======\n")
  cat("Using disease ID:", disease_id, "\n")
  
  # 1. Load data
  genes_long <- load_pathway_data(gmt_file)
  
  # Process OpenTargets JSON files for disease
  files.temp <- list.files(path = "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/results/drugtarget_test/associationByDatatypeDirect", pattern = "*.json",full.names=TRUE)
  files <- c(files.temp,list.files(path = "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/results/drugtarget_test/associationByDatatypeDirect", pattern = "*.json",full.names = TRUE))
  if(length(files) == 0) {
    stop("No JSON files found in current directory")
  }
  
  disease_data <- process_disease_data(files, disease_id)
  
    # Filter and format the gene-disease associations
  gene_disease_associations <- disease_data %>%
    select(targetId, diseaseId, datatypeId, score, evidenceCount) %>%
    arrange(desc(evidenceCount))
  
  # Create output file name
  gene_disease_file <- paste0(trait, "_", tool_base, "_gene_disease_associations.csv")
  
  # Save to CSV
  write.csv(gene_disease_associations, gene_disease_file, row.names = FALSE)
  cat("Saved", nrow(gene_disease_associations), "gene-disease associations to", gene_disease_file, "\n")
  
  # Select max evidence count for each target
  disease_targets <- disease_data %>%
    group_by(targetId) %>%
    slice_max(evidenceCount, with_ties = FALSE) %>%
    ungroup()
  
  # Load results files based on tool type
  birewire_data <- read.table(birewire_results_file,header=T)
  keeppath_data <- read.table(keeppathsize_results_file,header=T)

  cat("Loaded", nrow(birewire_data), "pathways from BireWire results\n")
  cat("Loaded", nrow(keeppath_data), "pathways from KeepPathSize results\n")
  
  # Calculate pathway scores
  pathway_scores <- calculate_pathway_scores(genes_long, disease_targets)
  
  # 2. Run analysis for each N value
  all_results <- list()
  
  for(n in n_values) {
    result_key <- paste0("n_", n)
    all_results[[result_key]] <- analyze_method_at_n(
      birewire_data, 
      keeppath_data, 
      pathway_scores, 
      n, 
      trait
    )
  }
  
  # 3. Create advantage calculations
  method_labels <- list(
    birewire = "Fix gene frequency and pathway size",
    keeppath = "Fix pathway size"
  )
  
  advantage_data <- calculate_method_advantages(all_results, n_values, method_labels)
  
  # Write detailed advantage data
  write.csv(advantage_data, paste0(trait, "_detailed_advantage.csv"), row.names=FALSE)
  
  # Create summary tables by ranking type
  for(ranking_type in unique(advantage_data$ranking_type)) {
    # Filter data for this ranking type
    type_data <- advantage_data[advantage_data$ranking_type == ranking_type,]
    
    # Get suffix for filename
    file_suffix <- switch(ranking_type,
                         "empirical_p" = "",
                         "raw_p" = "_raw_p",
                         "sig_beta" = "_sig_beta")
    
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
        write.csv(summary_data, paste0(trait, "_advantage_summary", file_suffix, ".csv"), row.names=FALSE)
        cat("Created summary file for ranking type:", ranking_type, "\n")
        
      }, error = function(e) {
        cat("Error creating summary for ranking type", ranking_type, ":", e$message, "\n")
        
        # Create a simple summary without pivot_wider
        simple_summary <- type_data %>%
          select(N, method, ranking_type, 
                 mean_score_advantage, mean_score_sig, 
                 evidence_density_advantage, evidence_density_sig)
        
        write.csv(simple_summary, paste0(trait, "_advantage_summary", file_suffix, "_simple.csv"), row.names=FALSE)
        cat("Created simplified summary file for ranking type:", ranking_type, "\n")
      })
    } else {
      cat("No data for ranking type:", ranking_type, "- skipping summary file\n")
    }
  }
  
  # 4. Write out the advantage data files (no visualization)
  write.csv(advantage_data, paste0(trait, "_detailed_advantage.csv"), row.names=FALSE)
  
  # Create summary table
    summary_data <- advantage_data %>%
    select(N, method, 
           mean_score_advantage,  t_statistic_mean_score, df_mean_score, mean_score_sig, 
           evidence_density_advantage, t_statistic_evidence_density, df_evidence_density, evidence_density_sig) %>%
    pivot_wider(
      names_from = method,
      values_from = c(mean_score_advantage, t_statistic_mean_score, df_mean_score, mean_score_sig, 
                     evidence_density_advantage, t_statistic_evidence_density, df_evidence_density, evidence_density_sig)
    )
  
  
  # Write summary data to CSV
  write.csv(summary_data, paste0(trait, "_advantage_summary.csv"), row.names=FALSE)
  
  # 5. Create a summary file for Nextflow tracking
  summary_file <- paste0(tolower(trait), "_", tolower(tool_base), "_size_matched_analysis_summary.tsv")
  write.table(
    data.frame(
      trait = trait,
      disease_id = disease_id,
      tool_base = tool_base,
      n_values = paste(n_values, collapse = ","),
      num_birewire_pathways = nrow(birewire_data),
      num_keeppath_pathways = nrow(keeppath_data)
    ),
    file = summary_file,
    sep = "\t", 
    row.names = FALSE, 
    quote = FALSE
  )
  
  cat("Summary written to", summary_file, "\n")
  cat("======= Analysis Complete =======\n")
}

# Execute the main function
main()