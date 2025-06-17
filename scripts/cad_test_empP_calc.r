## calculate empP, Zscore/n, Zscore/sqrt(n)

library(data.table)
library(GSA)
library(tidyverse)
library(ggplot2)

cad <- read.table('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_real/cad/cad_real_set.gsa.out',header=T)

read_files <- function(path) {
  setwd(path)
  files <- list.files(pattern="*gsa.out")
  ldf <- lapply(files, function(f) tryCatch(read.table(f, header=T), error=function(e) NULL))
  ldf <- ldf[!sapply(ldf, is.null)]  # Remove any NULL entries
  return(ldf)
}

calc_empirical_pvalues_fast <- function(real_results, background, random_method) {
  # Convert real results to data.table
  real_dt <- as.data.table(real_results)

  # Base path for random results
  base_path <- paste0('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_random/',
                     random_method,'/',background,'/cad')
  
  if(!dir.exists(base_path)) {
    stop("Directory not found: ", base_path)
  }
  
  # Step 1: Read all random results at once
  message("Reading all random results files...")
  setwd(base_path)
  random_files <- list.files(pattern="*gsa.out")
  total_files <- length(random_files)
  message("Found ", total_files, " random results files")
  
  # Create a master table of all random results
  message("Creating master random results table...")
  all_random_results <- rbindlist(
    lapply(random_files, function(file) {
      tryCatch({
        dt <- fread(file, select = c("FULL_NAME", "P"))
        dt[, file_id := file]  # Add file identifier
        return(dt)
      }, error = function(e) {
        message("Error reading file ", file, ": ", e$message)
        return(NULL)
      })
    }),
    fill = TRUE,
    use.names = TRUE
  )
  
  # Step 2: Calculate empirical p-values directly
  message("Calculating empirical p-values...")
  
  results <- real_dt %>% group_by(FULL_NAME) %>% mutate(
        n_more_extreme = sum(random_data$P <= P, na.rm = TRUE),
        empirical_pval = (n_more_extreme + 1) / (nrow(random_data) + 1),
        dataset = trait,
        rand = rand
      ) %>%
      ungroup()
    

  # Count occurrences per pathway where random P ≤ real P
  counts <- real_dt[, {
    pathway <- FULL_NAME
    real_p <- P
    
    # Extract relevant random results for this pathway
    random_subset <- all_random_results[FULL_NAME == pathway]
    
    # Count how many are more significant than real
    more_significant_count <- sum(random_subset$P <= real_p, na.rm = TRUE)
    
    # Count unique files this pathway appears in
    unique_files <- uniqueN(random_subset$file_id)
    
    # Calculate empirical p-value
    list(
      real_pvalue = real_p,
      more_significant_count = more_significant_count,
      total_random = unique_files,
      empirical_pvalue = (more_significant_count + 1) / (unique_files + 1)
    )
  }, by = FULL_NAME]
  
  message("Done! Processed ", uniqueN(real_pathways_dt$FULL_NAME), " pathways.")
  
  return(as.data.frame(counts))
}


# Calculate empirical p-values for different randomization methods
empirical_keeppathsize <- calc_empirical_pvalues_fast(cad, "msigdbgenes", "keeppathsize")
empirical_birewire <- calc_empirical_pvalues_fast(cad, "msigdbgenes", "birewire")

# Add method information
empirical_keeppathsize$method <- "keeppathsize"
empirical_birewire$method <- "birewire"

# Combine results
empirical_combined <- rbind(empirical_keeppathsize, empirical_birewire)

# Save results
write.csv(empirical_combined, "cad_empirical_pvalues.csv", row.names = FALSE)








calc_empP <- function(background, rand, trait, output_name = NULL) {
  
  cat("Processing", trait, "...\n")
  
  # Read real data - KEEP ALL COLUMNS
  real_file <- paste0('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_real/',trait,'/',trait,'_real_set.gsa.out')
    
  if (!file.exists(real_file)) {
    cat("Warning: Real data file not found for", trait, "\n")
    return(NULL)
  }
  
  real_data <- read.table(real_file, header = TRUE)
  
  # Read random data - NEED BOTH P and FULL_NAME for pathway matching
  random_files <- list.files(
    path = paste0('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_random/',rand,'/',background,'/',trait),
    pattern = "\\.gsa\\.out$", 
    full.names = TRUE
  )
  
  if (length(random_files) == 0) {
    cat("Warning: No random files found for", trait, "\n")
    return(NULL)
  }
  
  # Read and combine random data - keep FULL_NAME and P for pathway matching
  random_data_list <- list()
  for (i in seq_along(random_files)) {
    tryCatch({
      temp_data <- read.table(random_files[i], header = TRUE)
      if (all(c("P", "FULL_NAME") %in% colnames(temp_data))) {
        temp_data$random_iter <- i  # Track which random iteration
        random_data_list[[i]] <- temp_data[, c("FULL_NAME", "P", "random_iter")]
      }
    }, error = function(e) {
      cat("Warning: Could not read", random_files[i], "\n")
    })
  }
  
  if (length(random_data_list) == 0) {
    cat("Warning: No valid random data for", trait, "\n")
    return(NULL)
  }
  
  # Combine all random data
  random_data <- do.call(rbind, random_data_list)
  
  # Calculate empirical p-values PER PATHWAY
  trait_results <- real_data %>%
    rowwise() %>%
    mutate(
      # For each pathway, count how many random iterations had P <= real P
      n_more_extreme = sum(random_data$P[random_data$FULL_NAME == FULL_NAME] <= P, na.rm = TRUE),
      n_total_random = sum(random_data$FULL_NAME == FULL_NAME, na.rm = TRUE),
      empirical_pval = ifelse(n_total_random > 0, 
                             (n_more_extreme + 1) / (n_total_random + 1), 
                             NA)
    ) %>%
    ungroup() %>%
    # Add metadata columns
    mutate(
      dataset = trait,
      rand = rand,
      FPR = empirical_pval  # For backward compatibility
    ) %>%
    # Rename FULL_NAME to Var1 for compatibility with existing code
    rename(Var1 = FULL_NAME) %>%
    # Remove pathways where empirical p-value couldn't be calculated
    filter(!is.na(empirical_pval)) %>%
    # Sort by empirical p-value
    arrange(empirical_pval)
  
  # Assign to global environment if output_name is provided
  if (!is.null(output_name)) {
    assign(output_name, trait_results, envir = .GlobalEnv)
    cat("Results assigned to", output_name, "with", nrow(trait_results), "rows\n")
  }
  
  # Print summary
  cat("Completed", trait, "- found", nrow(trait_results), "pathways with empirical p-values\n")
  cat("Columns in output:\n")
  print(colnames(trait_results))
  
  # Print pathway matching summary
  cat("Random tests per pathway: mean =", round(mean(trait_results$n_total_random, na.rm = TRUE), 1),
      ", range =", range(trait_results$n_total_random, na.rm = TRUE), "\n")
  
  return(trait_results)
}

# Process individual traits
df_cad_birewire <- calc_empP("msigdbgenes", "birewire", "cad", 
                           output_name = "df_cad_birewire")

df_cad_keeppathsize <- calc_empP("msigdbgenes", "keeppathsize", "cad", 
                           output_name = "df_cad_keeppathsize")


# Combine results from both methods
df_cad_birewire$method <- "birewire"
df_cad_keeppathsize$method <- "keeppathsize"
df_cad_combined <- rbind(df_cad_birewire, df_cad_keeppathsize)

# calculate Zscore/n, Zscore/sqrt(n)
df_cad_combined$Zscore <- df_cad_combined$BETA/df_cad_combined$SE
df_cad_combined$Zscore_n <- df_cad_combined$Zscore / df_cad_combined$NGENES
df_cad_combined$Zscore_sqrt_n <- df_cad_combined$Zscore / sqrt(df_cad_combined$NGENES)

# Save combined results
write.csv(df_cad_combined, "cad_empirical_pvalues_combined.csv", row.names = FALSE)






