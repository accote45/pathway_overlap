process calc_empirical_pvalues {
    executor 'lsf'
    tag "${trait}_${tool}_empP"
    
    input:
    tuple val(trait),
          val(tool),
          path(real_results),
          path("random_results/*")
    
    output:
    tuple val(trait),
          val(tool), 
          path("${trait}_${tool}_empirical_pvalues.txt")
    
    publishDir "${params.outdir}/empirical_pvalues/${tool}/${trait}", mode: 'copy', overwrite: true
    
    script:
    """
    #!/bin/bash
    module load R
    
    Rscript - <<EOF
    library(tidyverse)
    library(data.table)
    
    # Tool-specific column mapping
    tool_config <- list(
      magma = list(
        pathway_col = "FULL_NAME",
        pval_col = "P",
        required_cols = c("FULL_NAME", "P", "NGENES", "BETA", "BETA_STD", "SE")
      ),
      prset = list(
        pathway_col = "Set",
        pval_col = "P",
        required_cols = c("Set", "P", "Coefficient", "R2", "P.adj")
      )
    )
    
    tool <- "${tool}"
    trait <- "${trait}"
    
    cat("Processing empirical p-values for", tool, "results from", trait, "\\n")
    
    # Get tool configuration
    if (!tool %in% names(tool_config)) {
      stop("Unsupported tool: ", tool)
    }
    
    config <- tool_config[[tool]]
    pathway_col <- config\$pathway_col
    pval_col <- config\$pval_col
    required_cols <- config\$required_cols
    
    # Read real results
    cat("Reading real results file: ${real_results}\\n")
    
    # Handle different file formats
    if (grepl("\\\\.gsa\\\\.out\$", "${real_results}")) {
      # MAGMA format - tab separated, has header
      real_data <- read.table("${real_results}", header = TRUE, sep = "\\t", stringsAsFactors = FALSE)
    } else if (grepl("\\\\.summary\$", "${real_results}")) {
      # PRSet format - space separated
      real_data <- read.table("${real_results}", header = TRUE, stringsAsFactors = FALSE)
    } else {
      # Generic format - try to auto-detect
      real_data <- fread("${real_results}", header = TRUE)
    }
    
    # Check required columns exist
    missing_cols <- required_cols[!required_cols %in% colnames(real_data)]
    if (length(missing_cols) > 0) {
      cat("Warning: Missing columns in real data:", paste(missing_cols, collapse = ", "), "\\n")
      required_cols <- required_cols[required_cols %in% colnames(real_data)]
    }
    
    cat("Found", nrow(real_data), "pathways in real results\\n")
    
    # Read random results
    random_files <- list.files("random_results", full.names = TRUE)
    cat("Found", length(random_files), "random result files\\n")
    
    if (length(random_files) == 0) {
      stop("No random result files found")
    }
    
    # Read and combine random data - only need pathway and p-value columns
    random_data_list <- list()
    for (i in seq_along(random_files)) {
      tryCatch({
        if (grepl("\\\\.gsa\\\\.out\$", random_files[i])) {
          temp_data <- read.table(random_files[i], header = TRUE, sep = "\\t", stringsAsFactors = FALSE)
        } else if (grepl("\\\\.summary\$", random_files[i])) {
          temp_data <- read.table(random_files[i], header = TRUE, stringsAsFactors = FALSE)
        } else {
          temp_data <- fread(random_files[i], header = TRUE)
        }
        
        if (all(c(pathway_col, pval_col) %in% colnames(temp_data))) {
          temp_data\$random_iter <- i
          random_data_list[[i]] <- temp_data[, c(pathway_col, pval_col, "random_iter")]
        }
      }, error = function(e) {
        cat("Warning: Could not read", random_files[i], "\\n")
      })
    }
    
    if (length(random_data_list) == 0) {
      stop("No valid random data files found")
    }
    
    # Combine random data
    random_data <- do.call(rbind, random_data_list)
    colnames(random_data)[1:2] <- c("pathway_name", "p_value")  # Standardize column names
    
    cat("Combined random data:", nrow(random_data), "rows from", length(unique(random_data\$random_iter)), "iterations\\n")
    
    # Calculate empirical p-values per pathway
    real_data\$pathway_name_std <- real_data[[pathway_col]]
    real_data\$p_value_std <- real_data[[pval_col]]
    
    empirical_results <- real_data %>%
      rowwise() %>%
      mutate(
        n_more_extreme = sum(random_data\$p_value[random_data\$pathway_name == pathway_name_std] <= p_value_std, na.rm = TRUE),
        n_total_random = sum(random_data\$pathway_name == pathway_name_std, na.rm = TRUE),
        empirical_pval = ifelse(n_total_random > 0, 
                               (n_more_extreme + 1) / (n_total_random + 1), 
                               NA),
        FPR = empirical_pval  # For backward compatibility
      ) %>%
      ungroup() %>%
      filter(!is.na(empirical_pval)) %>%
      arrange(empirical_pval) %>%
      mutate(
        trait = trait,
        tool = tool
      ) %>%
      select(-pathway_name_std, -p_value_std)  # Remove temporary columns
    
    # Write results
    fwrite(empirical_results, "${trait}_${tool}_empirical_pvalues.txt", sep = "\\t")
    
    cat("Calculated empirical p-values for", nrow(empirical_results), "pathways\\n")
    cat("Results written to ${trait}_${tool}_empirical_pvalues.txt\\n")
    
    # Print summary
    cat("Summary statistics:\\n")
    cat("Mean empirical p-value:", round(mean(empirical_results\$empirical_pval, na.rm = TRUE), 4), "\\n")
    cat("Median empirical p-value:", round(median(empirical_results\$empirical_pval, na.rm = TRUE), 4), "\\n")
    cat("FPR < 0.05:", sum(empirical_results\$FPR < 0.05, na.rm = TRUE), "pathways\\n")
    EOF
    """
}

process combine_empirical_results {
    executor 'local'
    
    input:
    path("results/*")
    
    output:
    path("all_empirical_pvalues_combined.txt")
    
    publishDir "${params.outdir}/empirical_pvalues/combined", mode: 'copy', overwrite: true
    
    script:
    """
    #!/bin/bash
    module load R
    
    Rscript - <<EOF
    library(data.table)
    
    # Find all empirical p-value files
    files <- list.files("results", pattern = "*_empirical_pvalues.txt", full.names = TRUE, recursive = TRUE)
    cat("Found", length(files), "empirical p-value files\\n")
    
    if (length(files) == 0) {
      stop("No empirical p-value files found")
    }
    
    # Read and combine all files
    all_results <- lapply(files, fread)
    combined <- rbindlist(all_results, fill = TRUE)
    
    # Write combined results
    fwrite(combined, "all_empirical_pvalues_combined.txt", sep = "\\t")
    
    cat("Combined results:", nrow(combined), "rows\\n")
    cat("Traits:", length(unique(combined\$trait)), "\\n")
    cat("Tools:", paste(unique(combined\$tool), collapse = ", "), "\\n")
    EOF
    """
}