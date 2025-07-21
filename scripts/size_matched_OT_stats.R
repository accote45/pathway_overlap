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
  
  # Ensure we have the necessary columns
  if(!"empirical_pval" %in% colnames(birewire_results)) {
    stop("Could not find empirical p-value column in BireWire results")
  }
  if(!"FULL_NAME" %in% colnames(birewire_results)) {
    stop("Could not find pathway name column in BireWire results")
  }
  if(!"empirical_pval" %in% colnames(keeppath_results)) {
    stop("Could not find empirical p-value column in KeepPathSize results")
  }
  if(!"FULL_NAME" %in% colnames(keeppath_results)) {
    stop("Could not find pathway name column in KeepPathSize results")
  }
  
  # Ensure empirical_pval is numeric
  birewire_results$empirical_pval <- as.numeric(as.character(birewire_results$empirical_pval))
  keeppath_results$empirical_pval <- as.numeric(as.character(keeppath_results$empirical_pval))
  
  # First, identify all unique pathways from both methods
  all_pathways <- unique(c(birewire_results$FULL_NAME, keeppath_results$FULL_NAME))
  
  # Create a dataset with method assignment and pathway size
  matching_data <- data.frame(
    name = all_pathways,
    in_birewire_top = all_pathways %in% head(birewire_results[order(birewire_results$empirical_pval),]$FULL_NAME, n),
    in_keeppath_top = all_pathways %in% head(keeppath_results[order(keeppath_results$empirical_pval),]$FULL_NAME, n)
  )
  
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
  
  # 2. Perform matching for BireWire
  bw_matched <- perform_matching(matching_data, "in_birewire_top")
  if(is.null(bw_matched)) {
    cat("Matching failed for BireWire top", n, "pathways\n")
    return(NULL)
  }
  bw_matched$N <- n
  bw_matched$method <- "BireWire"
  
  # 3. Calculate metrics for BireWire
  bw_comparison <- calculate_comparison_metrics(bw_matched, "in_birewire_top")
  size_balance_bw <- attr(bw_comparison, "size_balance")
  
  # 4. Perform matching for KeepPathSize
  kp_matched <- perform_matching(matching_data, "in_keeppath_top")
  if(is.null(kp_matched)) {
    cat("Matching failed for KeepPathSize top", n, "pathways\n")
    return(NULL)
  }
  kp_matched$N <- n
  kp_matched$method <- "KeepPathSize"
  
  # 5. Calculate metrics for KeepPathSize
  kp_comparison <- calculate_comparison_metrics(kp_matched, "in_keeppath_top")
  size_balance_kp <- attr(kp_comparison, "size_balance")
  
  # 6. Display results
  cat("\nSize-matched comparison (BireWire top", n, "vs. size-matched pathways):\n")
  cat("- Average size in BireWire top:", round(size_balance_bw$target_size, 1), 
      "vs. matched pathways:", round(size_balance_bw$control_size, 1), 
      "(", round(size_balance_bw$difference_pct, 1), "% diff)\n")
  
  for(i in 1:nrow(bw_comparison)) {
    metric_name <- bw_comparison$metric[i]
    target_val <- bw_comparison$target_value[i]
    control_val <- bw_comparison$control_value[i]
    diff <- bw_comparison$difference[i]
    p_val <- bw_comparison$p_value[i]
    
    cat("- Average", metric_name, "in BireWire top", n, ":", round(target_val, 4), 
        "vs. matched pathways:", round(control_val, 4), 
        "(diff:", round(diff, 4), ",", 
        ifelse(diff > 0, "BireWire better", "Control better"),
        ", p =", format.pval(p_val, digits=3), ")\n")
  }
  
  cat("\nSize-matched comparison (KeepPathSize top", n, "vs. size-matched pathways):\n")
  cat("- Average size in KeepPathSize top:", round(size_balance_kp$target_size, 1), 
      "vs. matched pathways:", round(size_balance_kp$control_size, 1),
      "(", round(size_balance_kp$difference_pct, 1), "% diff)\n")
  
  for(i in 1:nrow(kp_comparison)) {
    metric_name <- kp_comparison$metric[i]
    target_val <- kp_comparison$target_value[i]
    control_val <- kp_comparison$control_value[i]
    diff <- kp_comparison$difference[i]
    p_val <- kp_comparison$p_value[i]
    
    cat("- Average", metric_name, "in KeepPathSize top", n, ":", round(target_val, 4), 
        "vs. matched pathways:", round(control_val, 4), 
        "(diff:", round(diff, 4), ",", 
        ifelse(diff > 0, "KeepPathSize better", "Control better"),
        ", p =", format.pval(p_val, digits=3), ")\n")
  }
  
  # 7. Save matched data files
  write.csv(bw_matched, file = paste0(disease_name, "_", tool_base, "_n", n, "_birewire_matched.csv"))
  write.csv(kp_matched, file = paste0(disease_name, "_", tool_base, "_n", n, "_keeppath_matched.csv"))
  
  # 8. Return results for this N value
  return(list(
    birewire_matched = bw_matched,
    keeppath_matched = kp_matched,
    birewire_metrics = bw_comparison,
    keeppath_metrics = kp_comparison
  ))
}

calculate_method_advantages <- function(all_results, n_values, new_labels) {
  # Convert method names
  new_bw_label <- new_labels$birewire
  new_kp_label <- new_labels$keeppath
  
  advantage_data <- data.frame()
  
  for(n in n_values) {
    result_name <- paste0("n_", n)
    
    if(!is.null(all_results[[result_name]])) {
      bw_data <- all_results[[result_name]]$birewire_matched
      kp_data <- all_results[[result_name]]$keeppath_matched
      
      # Calculate BireWire advantage
      bw_top <- bw_data[bw_data$in_birewire_top == TRUE, ]
      bw_control <- bw_data[bw_data$in_birewire_top == FALSE, ]
      bw_advantage <- data.frame(
        N = factor(n, levels=n_values),
        method = new_bw_label,
        mean_score_top = mean(bw_top$mean_score, na.rm=TRUE),
        mean_score_control = mean(bw_control$mean_score, na.rm=TRUE),
        mean_score_advantage = mean(bw_top$mean_score, na.rm=TRUE) - mean(bw_control$mean_score, na.rm=TRUE),
        evidence_density_top = mean(bw_top$evidence_density, na.rm=TRUE),
        evidence_density_control = mean(bw_control$evidence_density, na.rm=TRUE),
        evidence_density_advantage = mean(bw_top$evidence_density, na.rm=TRUE) - mean(bw_control$evidence_density, na.rm=TRUE),
        p_value_score = all_results[[result_name]]$birewire_metrics$p_value[all_results[[result_name]]$birewire_metrics$metric == "mean_score"],
        t_statistic_score = all_results[[result_name]]$birewire_metrics$t_statistic[all_results[[result_name]]$birewire_metrics$metric == "mean_score"],
        df_score = all_results[[result_name]]$birewire_metrics$df[all_results[[result_name]]$birewire_metrics$metric == "mean_score"],
        p_value_density = all_results[[result_name]]$birewire_metrics$p_value[all_results[[result_name]]$birewire_metrics$metric == "evidence_density"],
        t_statistic_density = all_results[[result_name]]$birewire_metrics$t_statistic[all_results[[result_name]]$birewire_metrics$metric == "evidence_density"],
        df_density = all_results[[result_name]]$birewire_metrics$df[all_results[[result_name]]$birewire_metrics$metric == "evidence_density"]
      )
      
      # Calculate KeepPathSize advantage
      kp_top <- kp_data[kp_data$in_keeppath_top == TRUE, ]
      kp_control <- kp_data[kp_data$in_keeppath_top == FALSE, ]
      kp_advantage <- data.frame(
        N = factor(n, levels=n_values),
        method = new_kp_label,
        mean_score_top = mean(kp_top$mean_score, na.rm=TRUE),
        mean_score_control = mean(kp_control$mean_score, na.rm=TRUE),
        mean_score_advantage = mean(kp_top$mean_score, na.rm=TRUE) - mean(kp_control$mean_score, na.rm=TRUE),
        evidence_density_top = mean(kp_top$evidence_density, na.rm=TRUE),
        evidence_density_control = mean(kp_control$evidence_density, na.rm=TRUE),
        evidence_density_advantage = mean(kp_top$evidence_density, na.rm=TRUE) - mean(kp_control$evidence_density, na.rm=TRUE),
        p_value_score = all_results[[result_name]]$keeppath_metrics$p_value[all_results[[result_name]]$keeppath_metrics$metric == "mean_score"],
        t_statistic_score = all_results[[result_name]]$keeppath_metrics$t_statistic[all_results[[result_name]]$keeppath_metrics$metric == "mean_score"],
        df_score = all_results[[result_name]]$keeppath_metrics$df[all_results[[result_name]]$keeppath_metrics$metric == "mean_score"],
        p_value_density = all_results[[result_name]]$keeppath_metrics$p_value[all_results[[result_name]]$keeppath_metrics$metric == "evidence_density"],
        t_statistic_density = all_results[[result_name]]$keeppath_metrics$t_statistic[all_results[[result_name]]$keeppath_metrics$metric == "evidence_density"],
        df_density = all_results[[result_name]]$keeppath_metrics$df[all_results[[result_name]]$keeppath_metrics$metric == "evidence_density"]
      )
      
      advantage_data <- rbind(advantage_data, bw_advantage, kp_advantage)
    }
  }
  
  # Add significance indicators
  advantage_data$mean_score_sig <- ifelse(advantage_data$p_value_score < 0.05, "*", "")
  advantage_data$evidence_density_sig <- ifelse(advantage_data$p_value_density < 0.05, "*", "")
  
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
  
  # 4. Write out the advantage data files (no visualization)
  write.csv(advantage_data, paste0(trait, "_detailed_advantage.csv"), row.names=FALSE)
  
  # Create summary table
  summary_data <- advantage_data %>%
    select(N, method, 
           mean_score_advantage, t_statistic_score, df_score, mean_score_sig, 
           evidence_density_advantage, t_statistic_density, df_density, evidence_density_sig) %>%
    pivot_wider(
      names_from = method,
      values_from = c(mean_score_advantage, t_statistic_score, df_score, mean_score_sig, 
                     evidence_density_advantage, t_statistic_density, df_density, evidence_density_sig)
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