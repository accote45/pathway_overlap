#!/usr/bin/env Rscript

# Generate randomized GMT files using BiRewire
# Usage: Rscript generate_birewire_gmts.R <input_gmt> <output_dir> <num_random_sets>

library(BiRewire)
library(data.table)
library(GSA)
library(slam)  # For reading BiRewire's sparse matrix format

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript generate_birewire_gmts.R <input_gmt> <output_dir> <num_random_sets>")
}

input_gmt <- args[1]
output_dir <- args[2]
num_random_sets <- as.integer(args[3])

cat("========================================\n")
cat("BiRewire GMT Generation\n")
cat("========================================\n")
cat("Input GMT file:", input_gmt, "\n")
cat("Output directory:", output_dir, "\n")
cat("Number of random sets:", num_random_sets, "\n")
cat("========================================\n\n")

# Create output directory and temp directory for BiRewire output
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
birewire_temp_dir <- file.path(output_dir, "birewire_temp")
dir.create(birewire_temp_dir, recursive = TRUE, showWarnings = FALSE)

# Read GMT file and convert to long format
cat("Reading GMT file...\n")
dat <- GSA.read.gmt(input_gmt)
path.list <- dat$genesets
names(path.list) <- dat$geneset.names

real_genes_long <- rbindlist(lapply(names(path.list), function(name) {
  data.table(value = path.list[[name]], name = name)
}))

cat("Loaded", length(path.list), "pathways with", length(unique(real_genes_long$value)), "unique genes\n\n")

# Generate binary matrix (rows=genes, columns=pathways)
cat("Generating binary matrix...\n")
bin_matrix <- table(real_genes_long$value, real_genes_long$name)
cat("Binary matrix dimensions:", nrow(bin_matrix), "genes x", ncol(bin_matrix), "pathways\n\n")

# Run BiRewire sampling - saves to disk in numbered batch directories
cat("Running BiRewire sampling with K =", num_random_sets, "\n")
cat("This may take several hours for large gene sets.\n")
cat("Parameters: max.iter='n', accuracy=0.00005, verbose=TRUE\n")
cat("BiRewire will create numbered batch directories (1, 2, 3, ...) in:", birewire_temp_dir, "\n\n")

birewire.sampler.bipartite(
  bin_matrix,
  max.iter = "n",
  accuracy = 0.00005,
  verbose = TRUE,
  path = birewire_temp_dir,
  K = num_random_sets
)

cat("\n\nBiRewire sampling complete!\n")
cat("Loading randomized networks from batch directories and converting to GMT format...\n\n")

# BiRewire creates numbered directories (1, 2, 3, ...) with network files in sparse format
batch_dirs <- list.dirs(birewire_temp_dir, recursive = FALSE, full.names = TRUE)
batch_dirs <- batch_dirs[grepl("^[0-9]+$", basename(batch_dirs))]  # Only numeric directory names
batch_dirs <- batch_dirs[order(as.integer(basename(batch_dirs)))]  # Sort numerically

if (length(batch_dirs) == 0) {
  stop("No batch directories found in BiRewire output. Expected numbered directories (1, 2, 3, ...)")
}

cat("Found", length(batch_dirs), "BiRewire batch directories\n\n")

# Function to convert a sparse matrix to GMT format
write_sparse_to_gmt <- function(sparse_mat, output_file) {
  # Convert sparse matrix to dense matrix
  dense_mat <- as.matrix(sparse_mat)
  
  # Get labels
  pathway_names <- colnames(dense_mat)
  gene_names <- rownames(dense_mat)
  description_placeholder <- "PLACEHOLDER"
  
  # Open GMT file for writing
  gmt_con <- file(output_file, "w")
  
  # Loop over each pathway (column)
  for (j in 1:ncol(dense_mat)) {
    # Get genes (rows) that are non-zero
    gene_indices <- which(dense_mat[, j] != 0)
    if (length(gene_indices) > 0) {
      line <- c(pathway_names[j], description_placeholder, gene_names[gene_indices])
      writeLines(paste(line, collapse = "\t"), con = gmt_con)
    }
  }
  
  close(gmt_con)
}

# Process each batch directory
global_network_counter <- 1

for (batch_dir in batch_dirs) {
  # Find all network files in this batch directory (files without extensions)
  all_files <- list.files(batch_dir, full.names = TRUE)
  network_files <- all_files[!grepl("\\.(clabel|rlabel)$", all_files)]
  
  # Extract network numbers and sort
  network_numbers <- as.integer(sub(".*network_", "", basename(network_files)))
  network_files <- network_files[order(network_numbers)]
  
  cat("Processing batch directory:", basename(batch_dir), "with", length(network_files), "networks\n")
  
  for (network_file in network_files) {
    network_base <- basename(network_file)
    
    # Read the network using slam's CLUTO format reader
    tryCatch({
      # Read sparse matrix from BiRewire's CLUTO format
      sparse_mat <- read_stm_CLUTO(network_file)
      
      # Write GMT file with global counter
      gmt_file <- file.path(output_dir, paste0("GeneSet.random", global_network_counter, ".gmt"))
      write_sparse_to_gmt(sparse_mat, gmt_file)
      
      if (global_network_counter %% 100 == 0 || global_network_counter == 1) {
        cat("  Generated GMT file", global_network_counter, "of", num_random_sets, "\n")
      }
      
      global_network_counter <- global_network_counter + 1
      
    }, error = function(e) {
      warning("Error processing ", network_base, " in ", basename(batch_dir), ": ", e$message)
    })
  }
}

# Verify we generated the expected number of GMT files
actual_count <- global_network_counter - 1
if (actual_count != num_random_sets) {
  warning("Expected ", num_random_sets, " GMT files but generated ", actual_count)
}

# Clean up temporary BiRewire files
cat("\nCleaning up temporary BiRewire files...\n")
unlink(birewire_temp_dir, recursive = TRUE)

cat("\n========================================\n")
cat("SUCCESS!\n")
cat("Generated", actual_count, "randomized GMT files\n")
cat("Output location:", output_dir, "\n")
cat("========================================\n")