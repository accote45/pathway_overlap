# ===== Top Pathway Prioritization Analysis =====
# ===== Top Pathway Prioritization Analysis =====
# This script compares the biological relevance of top-ranked pathways 
# identified by different randomization methods

library(plyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(data.table)
library(jsonlite)
library(GSA)
library(MatchIt)

# Process command line arguments
args <- commandArgs(trailingOnly = TRUE)

trait <- args[1]
tool_base <- args[2]
birewire_results_file <- args[3]
keeppathsize_results_file <- args[4]
gmt_file <- args[5]

# Fix this parsing to ensure n_values becomes a numeric vector
n_values <- as.numeric(unlist(strsplit(args[6], ",")))
ts <- read.csv(args[7])

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

calculate_pathway_tissue_scores <- function(genes_long, tissue_expression_data) {
  cat("Calculating pathway-level tissue expression scores...\n")
  
  # Merge tissue specificity data with pathway genes
  # we exclude genes with no expression information
  masterfin <- merge(genes_long, tissue_expression_data, by.x="value", by.y="Name")
  
  # Get all tissue columns (excluding gene identifiers)
  tissue_cols <- setdiff(colnames(masterfin), c("value", "name", "targetId"))
  
  # Calculate pathway-level tissue expression metrics
  pathway_tissue_scores <- masterfin %>% 
    group_by(name) %>% 
    summarise(
      # Calculate metrics for each tissue
      across(all_of(tissue_cols), 
             list(
               mean = ~mean(.x, na.rm = TRUE),
               median = ~median(.x, na.rm = TRUE)
             ),
             .names = "{.col}_{.fn}"),
      # Count genes per pathway
      num_genes = n()
    ) %>%
    as.data.frame()
  
  # Add summary across all tissues for each pathway
  pathway_tissue_scores <- pathway_tissue_scores %>%
    rowwise() %>%
    mutate(
      # Find the tissue with maximum expression
      max_tissue = names(which.max(across(ends_with("_mean")))),
    ) %>%
    ungroup()
  
  cat("Generated expression profiles for", nrow(pathway_tissue_scores), "pathways across", 
      length(tissue_cols), "tissues\n")
  
  return(pathway_tissue_scores)
}

# === 2. Matching and Analysis Functions ===

prepare_matching_data <- function(birewire_results, keeppath_results, pathway_scores, n) {
  cat(paste("\nPreparing data for top", n, "pathways...\n"))
  
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

calculate_comparison_metrics <- function(matched_data, target_col) {
  # Extract group of interest and control group based on target column
  target_group <- matched_data[matched_data[[target_col]] == TRUE, ]
  control_group <- matched_data[matched_data[[target_col]] == FALSE, ]
  
    comparison <- data.frame()

    # Find all tissue mean columns (format: TissueName_mean)
    tissue_means <- grep("_mean$", colnames(matched_data), value=TRUE)
    tissue_medians <- grep("_median$", colnames(matched_data), value=TRUE)
    
    # Process tissue means
    for(metric in tissue_means) {
      if(metric %in% colnames(matched_data)) {
        # Extract tissue name from column name
        tissue_name <- sub("_mean$", "", metric)
        
        target_value <- mean(target_group[[metric]], na.rm=TRUE)
        control_value <- mean(control_group[[metric]], na.rm=TRUE)
        
        # Run t-test only if there's sufficient data
        if(sum(!is.na(target_group[[metric]])) >= 3 && sum(!is.na(control_group[[metric]])) >= 3) {
          t_test_result <- tryCatch({
            t.test(target_group[[metric]], control_group[[metric]])
          }, error = function(e) {
            list(statistic = NA, parameter = NA, p.value = NA)
          })
          
          comparison <- rbind(comparison, data.frame(
            metric = paste0(tissue_name, "_mean"),
            target_value = target_value,
            control_value = control_value,
            difference = target_value - control_value,
            percent_diff = 100 * (target_value - control_value) / ifelse(control_value == 0, 1, control_value),
            t_statistic = ifelse(is.na(t_test_result$statistic), NA, t_test_result$statistic),
            df = ifelse(is.na(t_test_result$parameter), NA, t_test_result$parameter),
            p_value = ifelse(is.na(t_test_result$p.value), NA, t_test_result$p.value),
            significant = ifelse(is.na(t_test_result$p.value), FALSE, t_test_result$p.value < 0.05),
            tissue = tissue_name
          ))
        }
      }
    }
    
    # Process tissue medians
    for(metric in tissue_medians) {
      if(metric %in% colnames(matched_data)) {
        # Extract tissue name from column name
        tissue_name <- sub("_median$", "", metric)
        
        target_value <- mean(target_group[[metric]], na.rm=TRUE)
        control_value <- mean(control_group[[metric]], na.rm=TRUE)
        
        # Run t-test only if there's sufficient data
        if(sum(!is.na(target_group[[metric]])) >= 3 && sum(!is.na(control_group[[metric]])) >= 3) {
          t_test_result <- tryCatch({
            t.test(target_group[[metric]], control_group[[metric]])
          }, error = function(e) {
            list(statistic = NA, parameter = NA, p.value = NA)
          })
          
          comparison <- rbind(comparison, data.frame(
            metric = paste0(tissue_name, "_median"),
            target_value = target_value,
            control_value = control_value,
            difference = target_value - control_value,
            percent_diff = 100 * (target_value - control_value) / ifelse(control_value == 0, 1, control_value),
            t_statistic = ifelse(is.na(t_test_result$statistic), NA, t_test_result$statistic),
            df = ifelse(is.na(t_test_result$parameter), NA, t_test_result$parameter),
            p_value = ifelse(is.na(t_test_result$p.value), NA, t_test_result$p.value),
            significant = ifelse(is.na(t_test_result$p.value), FALSE, t_test_result$p.value < 0.05),
            tissue = tissue_name
          ))
        }
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


# === 4. Individual Analysis Functions ===

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
    # 8. Return results for this N value
  return(list(
    birewire_matched = bw_matched,
    keeppath_matched = kp_matched,
    birewire_metrics = bw_comparison,
    keeppath_metrics = kp_comparison
  ))
}  

write_detailed_results <- function(all_results, trait, tool_base) {
  # Initialize an empty dataframe to store combined results
  all_metrics <- data.frame()
  
  # Loop through all N values
  for (n in c(10, 20, 50, 100)) {
    result_key <- paste0("n_", n)
    
    if (!is.null(all_results[[result_key]])) {
      # Get BireWire metrics and add N value
      if (!is.null(all_results[[result_key]]$birewire_metrics)) {
        bw_metrics <- all_results[[result_key]]$birewire_metrics
        bw_metrics$N <- n
        bw_metrics$method <- "birewire"
        
        # Add to combined results
        all_metrics <- rbind(all_metrics, bw_metrics)
      }
      
      # Get KeepPath metrics and add N value
      if (!is.null(all_results[[result_key]]$keeppath_metrics)) {
        kp_metrics <- all_results[[result_key]]$keeppath_metrics
        kp_metrics$N <- n
        kp_metrics$method <- "keeppath"
        
        # Add to combined results
        all_metrics <- rbind(all_metrics, kp_metrics)
      }
    }
  }
  
  # Write combined results to file
  combined_filename <- paste0(trait, "_", tool_base, "_all_detailed_metrics.csv")
  write.csv(all_metrics, combined_filename, row.names = FALSE)
  cat("All detailed metrics written to", combined_filename, "\n")
  
  return(all_metrics)
}

# === Main Execution ===
main <- function() {
  cat("======= Starting Size-Matched Analysis for", trait, "=======\n")
  
  # 1. Load data
  genes_long <- load_pathway_data(gmt_file)
    
  # Load results files based on tool type
  birewire_data <- read.table(birewire_results_file,header=T)
  keeppath_data <- read.table(keeppathsize_results_file,header=T)

  cat("Loaded", nrow(birewire_data), "pathways from BireWire results\n")
  cat("Loaded", nrow(keeppath_data), "pathways from KeepPathSize results\n")
  
  # Calculate pathway scores
  pathway_scores <- calculate_pathway_tissue_scores(genes_long, ts)
  
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
  
  # 3. write results to file
  detailed_metrics <- write_detailed_results(all_results, trait, tool_base)
  
  # 4. Create a summary file for Nextflow tracking
  summary_file <- paste0(trait, "_", tool_base, "_size_matched_analysis_summary.tsv")
  write.table(
    data.frame(
      trait = trait,
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