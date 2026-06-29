#!/usr/bin/env Rscript

# Generate randomized GMT files using KeepPathSize method
# Usage: Rscript generate_keeppathsize_gmts.R <input_gmt> <output_dir> <num_random_sets>

library(GSA)
library(data.table)
library(tidyverse)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 || length(args) > 4) {
  stop("Usage: Rscript generate_keeppathsize_gmts.R <input_gmt> <output_dir> <num_random_sets> [num_cores]")
}

input_gmt <- args[1]
output_dir <- args[2]
num_random_sets <- as.integer(args[3])
num_cores <- if (length(args) == 4) as.integer(args[4]) else 1L

cat("========================================\n")
cat("KeepPathSize GMT Generation\n")
cat("========================================\n")
cat("Input GMT file:", input_gmt, "\n")
cat("Output directory:", output_dir, "\n")
cat("Number of random sets:", num_random_sets, "\n")
cat("Number of cores:", num_cores, "\n")
cat("========================================\n\n")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read GMT file and convert to long format
cat("Reading GMT file...\n")
dat <- GSA.read.gmt(input_gmt)
path.list <- dat$genesets
names(path.list) <- dat$geneset.names

genes_long <- rbindlist(lapply(names(path.list), function(name) {
  data.table(value = path.list[[name]], name = name)
}))

cat("Loaded", length(path.list), "pathways with", length(unique(genes_long$value)), "unique genes\n\n")

# Calculate pathway sizes
set.size <- as.data.frame(table(genes_long$name))
colnames(set.size) <- c("pathway", "size")

cat("Pathway size distribution:\n")
cat("  Min:", min(set.size$size), "\n")
cat("  Max:", max(set.size$size), "\n")
cat("  Mean:", round(mean(set.size$size), 2), "\n")
cat("  Median:", median(set.size$size), "\n\n")

# Get all unique genes for random sampling
all_genes <- unique(genes_long$value)
cat("Total unique genes available for sampling:", length(all_genes), "\n\n")

# Function to generate one randomized GMT file
generate_random_gmt <- function(perm_number, set.size, all_genes, output_dir) {
  cat("Generating randomization", perm_number, "of", num_random_sets, "\n")
  
  # Initialize list to store randomized pathways
  randomized_groups <- vector("list", nrow(set.size))
  
  # For each pathway, randomly sample genes maintaining pathway size
  for (i in 1:nrow(set.size)) {
    group_size <- set.size$size[i]
    randomized_groups[[i]] <- sample(all_genes, group_size, replace = FALSE)
  }
  names(randomized_groups) <- set.size$pathway
  
  # Convert to long format
  master <- rbindlist(lapply(names(randomized_groups), function(name) {
    data.table(value = randomized_groups[[name]], name = name)
  }))
  
  # Function to format pathway as GMT string
  get_pathway_string <- function(name, genes) {
    paste(genes, sep = "\t", collapse = "\t") %>%
      paste(name, "PLACEHOLDER", ., sep = "\t", collapse = "\t") %>%
      return()
  }
  
  # Write GMT file
  gmt_file <- file.path(output_dir, paste0("GeneSet.random", perm_number, ".gmt"))
  fileConn <- file(gmt_file, open = "wt")
  
  for (pathway_name in set.size$pathway) {
    pathway_genes <- master[master$name == pathway_name, "value"] %>%
      unlist() %>%
      unique()
    get_pathway_string(pathway_name, pathway_genes) %>% 
      writeLines(., fileConn)
  }
  close(fileConn)
  
  return(gmt_file)
}

# Generate all randomized GMT files
cat("\n========================================\n")
cat("Starting KeepPathSize randomization\n")
cat("Parallelizing across", num_cores, "core(s)\n")
cat("========================================\n\n")

# Each random set is independent and each worker writes its own GMT file, so
# this loop is embarrassingly parallel. mclapply forks one process per core.
generated_files <- unlist(mclapply(
  1:num_random_sets,
  function(perm) generate_random_gmt(perm, set.size, all_genes, output_dir),
  mc.cores = num_cores,
  mc.preschedule = TRUE
))

# Surface any worker that errored (mclapply returns try-error objects rather
# than aborting the whole run).
failed <- which(!file.exists(generated_files) | is.na(generated_files))
if (length(failed) > 0) {
  warning(length(failed), " randomization(s) failed to produce a GMT file")
}

cat("\n========================================\n")
cat("Verification\n")
cat("========================================\n")

# Verify generated files
actual_count <- length(list.files(output_dir, pattern = "*.gmt"))
cat("Expected", num_random_sets, "GMT files\n")
cat("Generated", actual_count, "GMT files\n")

if (actual_count != num_random_sets) {
  warning("Expected ", num_random_sets, " GMT files but generated ", actual_count)
}

cat("\n========================================\n")
cat("SUCCESS!\n")
cat("Generated", actual_count, "randomized GMT files\n")
cat("Output location:", output_dir, "\n")
cat("========================================\n")