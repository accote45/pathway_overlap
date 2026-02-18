library(tidyverse)
library(data.table)

# Process command line arguments
args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]
tool_base <- args[2]
rand_method <- args[3]
random_dir <- args[4]

cat("Calculating FPR for trait:", trait, "tool:", tool_base, "method:", rand_method, "\n")
cat("Random directory:", random_dir, "\n")

# Tool-specific file pattern and significance column
get_tool_config <- function(tool_base) {
  switch(tool_base,
    "magma" = list(
      pattern = "*.gsa.out",
      pathway_col = "FULL_NAME",
      pval_col = "P",
      sig_threshold = 0.05
    ),
    "prset" = list(
      pattern = "*.summary", 
      pathway_col = "Set",
      pval_col = "Competitive.P",
      sig_threshold = 0.05
    ),
    "gsamixer" = list(
      pattern = "*go_test_enrich.csv",
      pathway_col = "GO", 
      pval_col = "p_value",  # This would need to be confirmed
      sig_threshold = 0.05
    ),
    stop("Unsupported tool: ", tool_base)
  )
}

config <- get_tool_config(tool_base)

# Function to read files based on tool type
read_random_files <- function(random_dir, pattern) {
  if (!dir.exists(random_dir)) {
    stop("Random directory does not exist: ", random_dir)
  }
  
  files <- list.files(random_dir, pattern = glob2rx(pattern), full.names = TRUE)
  
  if (length(files) == 0) {
    stop("No files found matching pattern '", pattern, "' in ", random_dir)
  }
  
  cat("Found", length(files), "random result files\n")
  
  # Read all files
  file_data <- list()
  for (i in seq_along(files)) {
    tryCatch({
      if (tool_base == "gsamixer") {
        data <- fread(files[i], header = TRUE)
      } else {
        data <- read.table(files[i], header = TRUE, stringsAsFactors = FALSE)
      }
      
      # Basic validation
      if (nrow(data) > 0 && config$pathway_col %in% colnames(data) && config$pval_col %in% colnames(data)) {
        file_data[[i]] <- data
      } else {
        cat("Warning: Skipping file", basename(files[i]), "- missing required columns\n")
      }
    }, error = function(e) {
      cat("Warning: Error reading file", basename(files[i]), ":", e$message, "\n")
    })
  }
  
  # Remove NULL entries
  file_data <- file_data[!sapply(file_data, is.null)]
  
  if (length(file_data) == 0) {
    stop("No valid files could be read")
  }
  
  cat("Successfully read", length(file_data), "files\n")
  return(file_data)
}

# Calculate FPR
calculate_pathway_fpr <- function(file_data, config) {
  # Combine all data
  all_data <- do.call(rbind, file_data)
  
  # Remove "Base" pathway if present (common in pathway analysis)
  if ("Base" %in% all_data[[config$pathway_col]]) {
    all_data <- all_data[all_data[[config$pathway_col]] != "Base", ]
    cat("Removed 'Base' pathway from analysis\n")
  }
  
  # Find significant pathways across all permutations
  significant_pathways <- all_data[[config$pathway_col]][all_data[[config$pval_col]] < config$sig_threshold]
  
  # Count frequency of each pathway being significant
  pathway_counts <- table(significant_pathways)
  
  # Calculate FPR (frequency / total number of files)
  total_files <- length(file_data)
  fpr_results <- data.frame(
    pathway_name = names(pathway_counts),
    significant_count = as.vector(pathway_counts),
    fpr = as.vector(pathway_counts) / total_files,
    total_permutations = total_files,
    trait = trait,
    tool = tool_base,
    randomization_method = rand_method,
    stringsAsFactors = FALSE
  )
  
  return(fpr_results)
}

# Main execution
tryCatch({
  # Read random files
  file_data <- read_random_files(random_dir, config$pattern)
  
  # Calculate FPR
  fpr_results <- calculate_pathway_fpr(file_data, config)
  
  # Write results
  output_file <- paste0(trait, "_", tool_base, "_", rand_method, "_fpr_results.csv")
  write.csv(fpr_results, output_file, row.names = FALSE, quote = FALSE)
  
  cat("FPR calculated for", nrow(fpr_results), "pathways\n")
  cat("Results written to:", output_file, "\n")
  
}, error = function(e) {
  cat("Error in FPR calculation:", e$message, "\n")
  quit(status = 1)
})
