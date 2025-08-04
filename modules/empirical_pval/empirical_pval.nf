process calc_empirical_pvalues {
    executor 'lsf'
    tag "${trait}_${tool}_empP"
    
    input:
    tuple val(trait),
          val(tool),
          path(real_results),
          val(random_dir)
    
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
    cat > calc_empirical.R << 'EOF'
    library(tidyverse)
    library(data.table)
    library(parallel)

    # Tool-specific column mapping
    tool_config <- list(
      magma = list(
        pathway_col = "FULL_NAME",
        pval_col = "P",
        beta_col = "BETA", 
        required_cols = c("FULL_NAME", "P", "NGENES", "BETA", "BETA_STD", "SE")
      ),
      prset = list(
        pathway_col = "Set",
        pval_col = "P",
        beta_col = "Coefficient",
        required_cols = c("Set", "P", "Coefficient", "R2", "P.adj")
      )
    )

    # Get command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    trait <- args[1]
    full_tool <- args[2] 
    real_results_file <- args[3]
    random_dir <- args[4]

    # Extract base tool name within R - much clearer!
    base_tool <- ifelse(grepl("_", full_tool), strsplit(full_tool, "_")[[1]][1], full_tool)

    cat("Processing empirical p-values and standardized effect sizes for", full_tool, "results from", trait, "\\n")

    # Get tool configuration using base tool name
    if (!base_tool %in% names(tool_config)) {
      stop("Unsupported base tool: ", base_tool)
    }

    config <- tool_config[[base_tool]]
    pathway_col <- config[["pathway_col"]]
    pval_col <- config[["pval_col"]]
    beta_col <- config[["beta_col"]]
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
    
    # Verify beta column exists
    if (!beta_col %in% colnames(real_data)) {
      stop("Required beta column '", beta_col, "' not found in real data")
    }

    cat("Found", nrow(real_data), "pathways in real results\\n")

    # Read random results
    cat("Reading random results from directory:", random_dir, "\\n")

    random_files <- list.files(random_dir, full.names = TRUE,pattern="*gsa.out")
    cat("Found", length(random_files), "random result files\\n")

    if (length(random_files) == 0) {
      stop("No random result files found")
    }

    # Read and combine random data - now including beta values
    random_data_list <- list()
    valid_file_count <- 0

    # Parallel processing setup
    num_cores <- min(detectCores() - 1, 8) # Use all but one core, max 8

    process_random_files <- function(file_batch, path_col, p_col, beta_col) {
      require(tidyverse)
      require(data.table)
      
      result_list <- list()
      count <- 0
      
      # Debug info at function start
      cat("Processing", length(file_batch), "files in batch\n")
      cat("Looking for columns:", path_col, p_col, beta_col, "\n")
      
      for (file_idx in seq_along(file_batch)) {
        file <- file_batch[file_idx]
        
        # Debug every 50 files
        if (file_idx %% 50 == 1) {
          cat("Processing file", file_idx, ":", basename(file), "\n")
        }
        
        # Check if file exists
        if (!file.exists(file)) {
          cat("File doesn't exist:", file, "\n")
          next
        }
        
        # Check file size
        file_size <- file.info(file)$size
        if (file_size == 0) {
          cat("Empty file:", file, "\n")
          next
        }
        
        tryCatch({
          is_magma <- endsWith(file, ".gsa.out")
          is_prset <- endsWith(file, ".summary")
          
          # Debug first file in detail
          if (file_idx == 1) {
            header <- readLines(file, n=1)
            cat("First file header:", header, "\n")
          }
          
          # Read the file
          if (is_magma || is_prset) {
            # Try fread first (much faster)
            temp_data <- tryCatch({
              fread(file, header = TRUE)
            }, error = function(e) {
              cat("fread failed for", basename(file), ":", e$message, "\n")
              # Fall back to read.table
              read.table(file, header = TRUE, stringsAsFactors = FALSE)
            })
            
            # Debug column names
            if (file_idx == 1) {
              cat("Columns in first file:", paste(colnames(temp_data), collapse=", "), "\n")
            }
            
            # Check if required columns exist
            if (all(c(path_col, p_col, beta_col) %in% colnames(temp_data))) {
              # Extract only the columns we need
              subset_data <- as.data.frame(temp_data)[, c(path_col, p_col, beta_col)]
              
              # Add iteration column
              count <- count + 1
              subset_data$random_iter <- count
              
              # Add to results
              result_list[[count]] <- subset_data
              
              if (count %% 100 == 0) {
                cat("Successfully processed", count, "files\n")
              }
            } else {
              missing <- setdiff(c(path_col, p_col, beta_col), colnames(temp_data))
              if (file_idx == 1) {
                cat("Missing columns in file", basename(file), ":", paste(missing, collapse=", "), "\n")
              }
            }
          }
        }, error = function(e) {
          cat("Error processing", basename(file), ":", e$message, "\n")
        })
      }
      
      cat("Batch finished, processed", count, "files successfully\n")
      
      # Return results
      if (length(result_list) > 0) {
        # Convert to data.table for efficiency
        combined <- rbindlist(result_list, fill=TRUE)
        
        # Rename columns to standardized names
        col_names <- names(combined)
        col_names[col_names == path_col] <- "pathway_name"
        col_names[col_names == p_col] <- "p_value" 
        col_names[col_names == beta_col] <- "beta_value"
        setnames(combined, col_names)
        
        cat("Returning", nrow(combined), "rows from", count, "files\n")
        return(combined)
      } else {
        cat("WARNING: No valid data found in this batch\n")
        return(NULL)
      }
    }

    # Use parallel processing for file reading
    if (length(random_files) > 100) {
      cat("Using parallel processing with", num_cores, "cores\n")
      chunks <- split(random_files, cut(seq_along(random_files), num_cores))
      random_chunks <- mclapply(chunks, process_random_files, 
                              path_col = pathway_col, p_col = pval_col, beta_col = beta_col,
                              mc.cores = num_cores)
      random_data <- rbindlist(random_chunks[!sapply(random_chunks, is.null)])
      
      # Ensure correct column names
      if (ncol(random_data) == 4) {
        setnames(random_data, names(random_data)[1:3], c("pathway_name", "p_value", "beta_value"))
      }
    } else {
      # Original sequential code for small file counts
      random_data_list <- list()
      valid_file_count <- 0
      for (i in seq_along(random_files)) {
        tryCatch({
          # Determine file type
          current_file <- random_files[i]
          is_magma_rand <- endsWith(current_file, ".gsa.out")
          is_prset_rand <- endsWith(current_file, ".summary")
          
          # Only process files that match our criteria
          if (is_magma_rand || is_prset_rand) {
            if (is_magma_rand) {
              temp_data <- read.table(current_file, header = TRUE, stringsAsFactors = FALSE)
            } else if (is_prset_rand) {
              temp_data <- read.table(current_file, header = TRUE, stringsAsFactors = FALSE)
            }
            
            # Make sure we have pathway, p-value, and beta columns
            if (all(c(pathway_col, pval_col, beta_col) %in% colnames(temp_data))) {
              valid_file_count <- valid_file_count + 1
              temp_data[["random_iter"]] <- valid_file_count  # Use sequential numbering for valid files
              # Add beta column to extracted data
              random_data_list[[valid_file_count]] <- temp_data[, c(pathway_col, pval_col, beta_col, "random_iter")]
            }
          }
        }, error = function(e) {
          cat("Warning: Could not read", current_file, "\\n")
        })
        
        # Print progress every 500 files
        if (i %% 500 == 0) {
          cat("Processed", i, "random files\\n")
        }
      }

      cat("Successfully read", length(random_data_list), "valid random files\\n")

      if (length(random_data_list) == 0) {
        stop("No valid random data files found")
      }

      # Combine random data
      random_data <- do.call(rbind, random_data_list)
      colnames(random_data)[1:3] <- c("pathway_name", "p_value", "beta_value")

      cat("Combined random data:", nrow(random_data), "rows\\n")
    }

    # Convert to data.table for efficiency
    random_dt <- as.data.table(random_data)
    real_dt <- as.data.table(real_data)

    # Add standardized column names to real data
    real_dt[, `:=`(
      pathway_name_std = get(pathway_col),
      p_value_std = get(pval_col),
      beta_value_std = get(beta_col)
    )]

    # Set keys for optimized joins
    setkey(random_dt, pathway_name)
    setkey(real_dt, pathway_name_std)

    # Calculate pathway counts in one go
    pathway_counts <- random_dt[, .(.N), by = pathway_name]
    setnames(pathway_counts, c("pathway_name", "n_total_random"))

    # Calculate beta statistics in one go
    beta_stats <- random_dt[, .(
      mean_beta_perm = mean(beta_value, na.rm = TRUE),
      sd_beta_perm = sd(beta_value, na.rm = TRUE),
      n_perms = .N  # Add number of permutations per pathway
    ), by = pathway_name]

    # Join counts and beta stats to real data efficiently
    results_dt <- merge(real_dt, pathway_counts, by.x = "pathway_name_std", by.y = "pathway_name", all.x = TRUE)
    results_dt <- merge(results_dt, beta_stats, by.x = "pathway_name_std", by.y = "pathway_name", all.x = TRUE)

    # Calculate p-values and effect sizes in one pass
    results_dt[, `:=`(
      n_more_extreme = {
        # For each pathway in results_dt
        mapply(function(path, p_val) {
          # Find matching permutations and count more extreme p-values
          sum(random_dt[.(path), p_value <= p_val], na.rm = TRUE)
        }, pathway_name_std, p_value_std)
      },
      empirical_pval = {
        # Will be calculated after n_more_extreme is available
        NA_real_
      },
      std_effect_size = {
        # Calculate standardized effect size
        ifelse(!is.na(sd_beta_perm) & sd_beta_perm > 0,
               (beta_value_std - mean_beta_perm) / sd_beta_perm,
               NA_real_)
      }
    )]

    # Calculate empirical p-value
    results_dt[, empirical_pval := (n_more_extreme + 1) / (n_total_random + 1)]

    # Only create final tibble once
    empirical_results <- results_dt[!is.na(empirical_pval),
                                   .(trait = trait, 
                                     tool = full_tool,
                                     pathway_name = get(pathway_col),
                                     p_value = get(pval_col),
                                     beta_value = get(beta_col),
                                     empirical_pval,
                                     std_effect_size,
                                     mean_beta_perm,
                                     sd_beta_perm,
                                     n_perms,
                                     n_more_extreme,
                                     n_total_random)]

    # Write results
    output_file <- paste0(trait, "_", full_tool, "_empirical_pvalues.txt")
    fwrite(empirical_results, output_file, sep = "\\t")

    cat("Calculated empirical p-values and standardized effect sizes for", nrow(empirical_results), "pathways\\n")
    cat("Results written to", output_file, "\\n")

    # Print summary
    cat("Summary statistics:\\n")
    cat("Mean empirical p-value:", round(mean(empirical_results$empirical_pval, na.rm = TRUE), 4), "\\n")
    cat("Median empirical p-value:", round(median(empirical_results$empirical_pval, na.rm = TRUE), 4), "\\n")
    cat("Mean standardized effect size:", round(mean(empirical_results$std_effect_size, na.rm = TRUE), 4), "\\n")
    
    # Print top 10 pathways by standardized effect size
    cat("\\nTop 10 pathways by standardized effect size:\\n")
    empirical_results %>%
      arrange(desc(std_effect_size)) %>%
      head(10) %>%
      select(all_of(c(pathway_col, pval_col, beta_col, "empirical_pval", "std_effect_size"))) %>%
      print(n = 10)
    EOF

    # Run the R script with both the full tool name and base tool name
    Rscript calc_empirical.R "${trait}" "${tool}" "${real_results}" "${random_dir}"
    """
}