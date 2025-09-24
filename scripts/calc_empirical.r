library(tidyverse)
    library(data.table)
    library(parallel)

    # Tool-specific column mapping
    tool_config <- list(
      magma = list(
        pathway_col = "FULL_NAME",
        pval_col = "P",
        beta_col = "BETA",
        ngenes_col = "NGENES",
        required_cols = c("FULL_NAME", "P", "NGENES", "BETA", "BETA_STD", "SE"),
        calc_pvalue = TRUE
      ),
      prset = list(
        pathway_col = "Set",
        pval_col = "P",
        beta_col = "Coefficient",
        ngenes_col = "nGenes",
        required_cols = c("Set", "P", "Coefficient", "PRS.R2", "Competitive.P"),
        calc_pvalue = TRUE
      ),
      gsamixer = list(
        pathway_col = "GO",
        pval_col = NULL,  # No p-value column for GSA-Mixer
        beta_col = "enrich",  # Use enrichment score instead
        ngenes_col = "NGENES",
        required_cols = c("GO", "enrich", "NGENES", "h2", "se_enrich"),
        calc_pvalue = FALSE  # Don't calculate empirical p-values for GSA-Mixer
      )
    )

    # Get command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    trait <- args[1]
    full_tool <- args[2] 
    real_results_file <- args[3]
    random_dir <- args[4]
    
    # Add GMT file path as 5th argument
    gmt_path <- if (length(args) >= 5) args[5] else NULL

    # Extract base tool name within R - much clearer!
    base_tool <- ifelse(grepl("_", full_tool), strsplit(full_tool, "_")[[1]][1], full_tool)

    cat("Processing", ifelse(base_tool == "gsamixer", "standardized effect sizes", "empirical p-values and standardized effect sizes"), 
        "for", full_tool, "results from", trait, "\n")

    # Get tool configuration using base tool name
    if (!base_tool %in% names(tool_config)) {
      stop("Unsupported base tool: ", base_tool)
    }

    config <- tool_config[[base_tool]]
    pathway_col <- config[["pathway_col"]]
    pval_col <- config[["pval_col"]]
    beta_col <- config[["beta_col"]]
    ngenes_col <- config[["ngenes_col"]]
    required_cols <- config[["required_cols"]]
    calc_pvalue <- config[["calc_pvalue"]]

    # Read GMT file to get valid pathway names when dealing with GSA-Mixer
    valid_pathways <- NULL
    if (base_tool == "gsamixer") {
      # For GSA-Mixer, GMT file is required
      if (is.null(gmt_path) || !file.exists(gmt_path)) {
        stop("Error: GMT file is required for GSA-Mixer analysis but was not found at path: ", gmt_path)
      }
      
      # Read the GMT file
      cat("Reading pathway definitions from GMT file:", gmt_path, "\n")
      gmt_content <- readLines(gmt_path)
      valid_pathways <- sapply(strsplit(gmt_content, "\t"), function(x) x[1])
      cat("Found", length(valid_pathways), "valid pathways in GMT file\n")
    }

    # Read real results - determine file type by extension
    cat("Reading real results file:", real_results_file, "\n")

    # Check file extension
    is_magma <- endsWith(real_results_file, ".gsa.out")
    is_prset <- endsWith(real_results_file, ".summary")
    is_gsamixer <- endsWith(real_results_file, ".go_test_enrich.csv")

    if (is_magma) {
      real_data <- read.table(real_results_file, header = TRUE, stringsAsFactors = FALSE)
    } else if (is_prset) {
      real_data <- read.table(real_results_file, header = TRUE, stringsAsFactors = FALSE)
    } else if (is_gsamixer) {
      real_data <- fread(real_results_file, header = TRUE)
      # For gsamixer, filter out rows without valid pathway names (skip base and coding_genes)
      real_data <- real_data[!is.na(real_data[[pathway_col]]) & 
                            real_data[[pathway_col]] != "base" & 
                            real_data[[pathway_col]] != "coding_genes"]
      
      # Filter pathways based on GMT file - strictly required for GSA-Mixer
      original_count <- nrow(real_data)
      real_data <- real_data[real_data[[pathway_col]] %in% valid_pathways]
      cat("Filtered GSA-Mixer results using GMT pathways:", 
          original_count, "->", nrow(real_data), "pathways\n")
      
      # If no pathways remain after filtering, stop the script
      if (nrow(real_data) == 0) {
        stop("No valid pathways found after filtering GSA-Mixer results with GMT file")
      }
    } else {
      real_data <- fread(real_results_file, header = TRUE)
    }

    # Check required columns exist
    missing_cols <- required_cols[!required_cols %in% colnames(real_data)]
    if (length(missing_cols) > 0) {
      cat("Warning: Missing columns in real data:", paste(missing_cols, collapse = ", "), "\n")
      required_cols <- required_cols[required_cols %in% colnames(real_data)]
    }
    
    # Verify beta column exists (enrichment for gsamixer)
    if (!beta_col %in% colnames(real_data)) {
      stop("Required beta/enrichment column '", beta_col, "' not found in real data")
    }

    cat("Found", nrow(real_data), "pathways in real results\n")

    # Read random results
    cat("Reading random results from directory:", random_dir, "\n")

    # Set pattern based on tool type
    file_pattern <- if(is_gsamixer) {
      "*go_test_enrich.csv"
    } else if(is_prset) {
      "*.summary"
    } else {
      "*gsa.out"
    }
    random_files <- list.files(random_dir, full.names = TRUE, pattern = file_pattern)
    cat("Found", length(random_files), "random result files\n")

    if (length(random_files) == 0) {
      stop("No random result files found")
    }

    # Read and combine random data - now including beta values
    random_data_list <- list()
    valid_file_count <- 0

    # Parallel processing setup
    num_cores <- min(detectCores() - 1, 8) # Use all but one core, max 8

    process_random_files <- function(file_batch, path_col, p_col, beta_col, valid_pathways) {
      require(tidyverse)
      require(data.table)
      
      result_list <- list()
      count <- 0
      
      # Debug info at function start
      cat("Processing", length(file_batch), "files in batch\n")
      if (!is.null(p_col)) {
        cat("Looking for columns:", path_col, p_col, beta_col, "\n")
      } else {
        cat("Looking for columns:", path_col, beta_col, "(no p-value column)\n")
      }
      
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
          is_gsamixer <- endsWith(file, ".go_test_enrich.csv")
          
          # Debug first file in detail
          if (file_idx == 1) {
            header <- readLines(file, n=1)
            cat("First file header:", header, "\n")
          }
          
          # Read the file
          if (is_magma || is_prset || is_gsamixer) {
            # Try fread first (much faster)
            temp_data <- tryCatch({
              if (is_magma) {
                # Always use read.table for MAGMA files
                read.table(file, header = TRUE, stringsAsFactors = FALSE)
              } else {
                # Try fread for non-MAGMA files
                fread(file, header = TRUE)
              }
            }, error = function(e) {
              cat("File reading failed for", basename(file), ":", e$message, "\n")
              # Fall back to read.table as a last resort
              read.table(file, header = TRUE, stringsAsFactors = FALSE)
            })
            
            # Debug column names
            if (file_idx == 1) {
              cat("Columns in first file:", paste(colnames(temp_data), collapse=", "), "\n")
            }
            
            # For GSA-Mixer, filter out base and coding_genes rows and only keep valid pathways from GMT
            if (is_gsamixer) {
              temp_data <- temp_data[temp_data[[path_col]] != "base" & 
                                   temp_data[[path_col]] != "coding_genes" & 
                                   temp_data[[path_col]] %in% valid_pathways, ]
            }
            
            # Check if required columns exist
            if (is.null(p_col)) {
              # For GSA-Mixer, we only need pathway and beta columns
              cols_exist <- all(c(path_col, beta_col) %in% colnames(temp_data))
            } else {
              cols_exist <- all(c(path_col, p_col, beta_col) %in% colnames(temp_data))
            }
            
            if (cols_exist) {
              # Extract only the columns we need
              if (is.null(p_col)) {
                # For GSA-Mixer, we don't need p-values
                subset_data <- as.data.frame(temp_data)[, c(path_col, beta_col)]
              } else {
                subset_data <- as.data.frame(temp_data)[, c(path_col, p_col, beta_col)]
              }
              
              # Add iteration column
              count <- count + 1
              subset_data$random_iter <- count
              
              # Add to results
              result_list[[count]] <- subset_data
              
              if (count %% 100 == 0) {
                cat("Successfully processed", count, "files\n")
              }
            } else {
              missing <- setdiff(c(path_col, if(is.null(p_col)) c() else p_col, beta_col), colnames(temp_data))
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
        if (!is.null(p_col) && p_col %in% col_names) {
          col_names[col_names == p_col] <- "p_value"
        }
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
                              valid_pathways = valid_pathways,
                              mc.cores = num_cores)
      random_data <- rbindlist(random_chunks[!sapply(random_chunks, is.null)])
      
      # Ensure correct column names for GSA-Mixer
      if (is.null(pval_col)) {
        if (!"p_value" %in% names(random_data)) {
          random_data$p_value <- NA_real_  # Add dummy column for consistency
        }
      } else if (ncol(random_data) >= 3) {
        setnames(random_data, names(random_data)[1:min(3,ncol(random_data))], 
                c("pathway_name", "p_value", "beta_value")[1:min(3,ncol(random_data))])
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
          is_gsamixer_rand <- endsWith(current_file, ".go_test_enrich.csv")
          
          # Only process files that match our criteria
          if (is_magma_rand || is_prset_rand || is_gsamixer_rand) {
            if (is_gsamixer_rand) {
              temp_data <- fread(current_file, header = TRUE)
              # Filter GSA-Mixer data to only include valid pathways from GMT
              temp_data <- temp_data[temp_data[[pathway_col]] %in% valid_pathways]
            } else {
              temp_data <- read.table(current_file, header = TRUE, stringsAsFactors = FALSE)
            }
            
            # Make sure we have pathway and beta columns (p-value optional for GSA-Mixer)
            if (is.null(pval_col)) {
              # For GSA-Mixer, we only need pathway and beta columns
              if (all(c(pathway_col, beta_col) %in% colnames(temp_data))) {
                valid_file_count <- valid_file_count + 1
                temp_data[["random_iter"]] <- valid_file_count
                random_data_list[[valid_file_count]] <- temp_data[, c(pathway_col, beta_col, "random_iter")]
              }
            } else if (all(c(pathway_col, pval_col, beta_col) %in% colnames(temp_data))) {
              valid_file_count <- valid_file_count + 1
              temp_data[["random_iter"]] <- valid_file_count
              random_data_list[[valid_file_count]] <- temp_data[, c(pathway_col, pval_col, beta_col, "random_iter")]
            }
          }
        }, error = function(e) {
          cat("Warning: Could not read", current_file, "\n")
        })
        
        # Print progress every 500 files
        if (i %% 500 == 0) {
          cat("Processed", i, "random files\n")
        }
      }

      cat("Successfully read", length(random_data_list), "valid random files\n")

      if (length(random_data_list) == 0) {
        stop("No valid random data files found")
      }

      # Combine random data - handling GSA-Mixer case separately
      random_data <- do.call(rbind, random_data_list)
      if (is.null(pval_col)) {
        # For GSA-Mixer, only two columns (pathway and beta)
        colnames(random_data)[1:2] <- c("pathway_name", "beta_value")
        # Add dummy p_value column for consistency
        random_data$p_value <- NA_real_
      } else {
        colnames(random_data)[1:3] <- c("pathway_name", "p_value", "beta_value")
      }

      cat("Combined random data:", nrow(random_data), "rows\n")
    }

    # Convert to data.table for efficiency
    random_dt <- as.data.table(random_data)
    real_dt <- as.data.table(real_data)

    # Add standardized column names to real data
    real_dt[, `:=`(
      pathway_name_std = get(pathway_col),
      p_value_std = if (!is.null(pval_col)) get(pval_col) else NA_real_,
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

    # Calculate effect sizes and p-values (if applicable) in one pass
    results_dt[, `:=`(
      std_effect_size = {
        # Calculate standardized effect size
        ifelse(!is.na(sd_beta_perm) & sd_beta_perm > 0,
               (beta_value_std - mean_beta_perm) / sd_beta_perm,
               NA_real_)
      }
    )]
    
    # Only calculate empirical p-values for tools that need them
    if (calc_pvalue) {
      results_dt[, `:=`(
        n_more_extreme = {
          # For each pathway in results_dt
          mapply(function(path, p_val) {
            # Find matching permutations and count more extreme p-values
            sum(random_dt[.(path), p_value <= p_val], na.rm = TRUE)
          }, pathway_name_std, p_value_std)
        },
        empirical_pval = NA_real_  # Will be calculated next
      )]
      
      # Calculate empirical p-value
      results_dt[, empirical_pval := (n_more_extreme + 1) / (n_total_random + 1)]
    }

    # Create final output - handle GSA-Mixer specially
    if (base_tool == "gsamixer") {
      empirical_results <- results_dt[!is.na(std_effect_size),
                                    .(trait = trait, 
                                      tool = full_tool,
                                      pathway_name = get(pathway_col),
                                      ngenes = get(ngenes_col),
                                      enrich = get(beta_col),  # Original enrichment score
                                      std_effect_size,
                                      mean_beta_perm,
                                      sd_beta_perm,
                                      n_perms,
                                      n_total_random)]
    } else {
      empirical_results <- results_dt[!is.na(empirical_pval),
                                    .(trait = trait, 
                                      tool = full_tool,
                                      pathway_name = get(pathway_col),
                                      ngenes = get(ngenes_col),
                                      p_value = get(pval_col),
                                      beta_value = get(beta_col),
                                      empirical_pval,
                                      std_effect_size,
                                      mean_beta_perm,
                                      sd_beta_perm,
                                      n_perms,
                                      n_more_extreme,
                                      n_total_random)]
    }

    # Write results
    output_file <- paste0(trait, "_", full_tool, "_", 
                         if(base_tool == "gsamixer") "standardized_effects.txt" else "empirical_pvalues.txt")
    fwrite(empirical_results, output_file, sep = "\t")

    cat("Results calculated for", nrow(empirical_results), "pathways\n")
    cat("Results written to", output_file, "\n")