process calc_empirical_pvalues {
    executor 'lsf'
    tag "${trait}_${tool}_empP"
    
    input:
    tuple val(trait),
          val(tool),
          path(real_results),
          path(random_dir)
    
    output:
    tuple val(trait),
          val(tool), 
          path("${trait}_${tool}_empirical_pvalues.txt")
    
    publishDir "${params.outdir}/empirical_pvalues/${tool}/${trait}", mode: 'copy', overwrite: true
    
    script:
    """
    #!/bin/bash
    module load R
    
    # Create R script
    cat > calc_empirical.R << 'RSCRIPT'
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
      ),
      gsea = list(
        pathway_col = "pathway",
        pval_col = "pval", 
        required_cols = c("pathway", "pval", "ES", "NES", "size")
      )
    )
    
    # Get command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    tool <- args[1]
    trait <- args[2]
    real_results_file <- args[3]
    random_dir <- args[4]
    
    cat("Processing empirical p-values for", tool, "results from", trait, "\\n")
    
    # Get tool configuration
    if (!tool %in% names(tool_config)) {
      stop("Unsupported tool: ", tool)
    }
    
    config <- tool_config[[tool]]
    pathway_col <- config[["pathway_col"]]
    pval_col <- config[["pval_col"]]
    required_cols <- config[["required_cols"]]
    
    # Read real results - determine file type by extension
    cat("Reading real results file:", real_results_file, "\\n")
    
    # Check file extension
    is_magma <- endsWith(real_results_file, ".gsa.out")
    is_prset <- endsWith(real_results_file, ".summary")
    
    if (is_magma) {
      real_data <- read.table(real_results_file, header = TRUE, stringsAsFactors = FALSE)
    } else if (is_prset) {
      real_data <- read.table(real_results_file, header = TRUE, stringsAsFactors = FALSE)
    } else {
      real_data <- fread(real_results_file, header = TRUE)
    }
    
    # Check required columns exist
    missing_cols <- required_cols[!required_cols %in% colnames(real_data)]
    if (length(missing_cols) > 0) {
      cat("Warning: Missing columns in real data:", paste(missing_cols, collapse = ", "), "\\n")
      required_cols <- required_cols[required_cols %in% colnames(real_data)]
    }
    
    cat("Found", nrow(real_data), "pathways in real results\\n")
    
    # Read random results
    cat("Reading random results from directory:", random_dir, "\\n")

    random_files <- list.files(random_dir, full.names = TRUE)
    cat("Found", length(random_files), "random result files\\n")
    
    if (length(random_files) == 0) {
      stop("No random result files found")
    }
    
    # Read and combine random data
    random_data_list <- list()
    for (i in seq_along(random_files)) {
      tryCatch({
        # Determine file type
        current_file <- random_files[i]
        is_magma_rand <- endsWith(current_file, ".gsa.out")
        is_prset_rand <- endsWith(current_file, ".summary")
        
        if (is_magma_rand) {
          temp_data <- read.table(current_file, header = TRUE, stringsAsFactors = FALSE)
        } else (is_prset_rand) {
          temp_data <- read.table(current_file, header = TRUE, stringsAsFactors = FALSE)
        }
        
        if (all(c(pathway_col, pval_col) %in% colnames(temp_data))) {
          temp_data[["random_iter"]] <- i
          random_data_list[[i]] <- temp_data[, c(pathway_col, pval_col, "random_iter")]
        }
      }, error = function(e) {
        cat("Warning: Could not read", current_file, "\\n")
      })
    }
    
    if (length(random_data_list) == 0) {
      stop("No valid random data files found")
    }
    
    # Combine random data
    random_data <- do.call(rbind, random_data_list)
    colnames(random_data)[1:2] <- c("pathway_name", "p_value")
    
    cat("Combined random data:", nrow(random_data), "rows\\n")
    
    # Calculate empirical p-values per pathway
    real_data[["pathway_name_std"]] <- real_data[[pathway_col]]
    real_data[["p_value_std"]] <- real_data[[pval_col]]
    
    empirical_results <- real_data %>%
      rowwise() %>%
      mutate(
        n_more_extreme = sum(random_data[["p_value"]][random_data[["pathway_name"]] == pathway_name_std] <= p_value_std, na.rm = TRUE),
        n_total_random = sum(random_data[["pathway_name"]] == pathway_name_std, na.rm = TRUE),
        empirical_pval = ifelse(n_total_random > 0, 
                               (n_more_extreme + 1) / (n_total_random + 1), 
                               NA),
        FPR = empirical_pval
      ) %>%
      ungroup() %>%
      filter(!is.na(empirical_pval)) %>%
      arrange(empirical_pval) %>%
      mutate(
        trait = trait,
        tool = tool
      ) %>%
      select(-pathway_name_std, -p_value_std)
    
    # Write results
    output_file <- paste0(trait, "_", tool, "_empirical_pvalues.txt")
    fwrite(empirical_results, output_file, sep = "\\t")
    
    cat("Calculated empirical p-values for", nrow(empirical_results), "pathways\\n")
    cat("Results written to", output_file, "\\n")
    
    # Print summary
    cat("Summary statistics:\\n")
    cat("Mean empirical p-value:", round(mean(empirical_results[["empirical_pval"]], na.rm = TRUE), 4), "\\n")
    cat("Median empirical p-value:", round(median(empirical_results[["empirical_pval"]], na.rm = TRUE), 4), "\\n")
    cat("FPR < 0.05:", sum(empirical_results[["FPR"]] < 0.05, na.rm = TRUE), "pathways\\n")
    
RSCRIPT
    
    # Run the R script
    Rscript calc_empirical.R "${tool}" "${trait}" "${real_results}" "${random_dir}"
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