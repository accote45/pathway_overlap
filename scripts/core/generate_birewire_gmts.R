#!/usr/bin/env Rscript

# Generate randomized GMT files using BiRewire
# Usage: Rscript generate_birewire_gmts.R <input_gmt> <output_dir> <num_random_sets>

library(BiRewire)
library(data.table)
library(GSA)

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

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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

# Run BiRewire sampling
cat("Running BiRewire sampling with K =", num_random_sets, "\n")
cat("This may take several hours for large gene sets.\n")
cat("Parameters: max.iter='n', accuracy=0.00005, step=5000\n\n")

step <- 5000
random_networks <- birewire.sampler.bipartite(
  bin_matrix,
  max.iter = "n",
  accuracy = 0.00005,
  verbose = TRUE,
  K = num_random_sets
)

cat("\n\nBiRewire sampling complete!\n")
cat("Converting randomized networks to GMT format...\n\n")

# Convert each randomized network to GMT format
for (i in 1:num_random_sets) {
  if (i %% 100 == 0 || i == 1) {
    cat("Processing random set", i, "of", num_random_sets, "\n")
  }
  
  # Get the randomized matrix
  random_mat <- random_networks[[i]]
  
  # Prepare pathway names and gene names
  pathway_names <- colnames(random_mat)
  gene_names <- rownames(random_mat)
  description_placeholder <- "PLACEHOLDER"
  
  # Write GMT file
  gmt_file <- file.path(output_dir, paste0("GeneSet.random", i, ".gmt"))
  gmt_con <- file(gmt_file, "w")
  
  # Loop over each pathway (column)
  for (j in 1:ncol(random_mat)) {
    # Get the genes (rows) that are non-zero
    gene_indices <- which(random_mat[, j] != 0)
    if (length(gene_indices) > 0) {
      # Prepare the line in GMT format
      line <- c(pathway_names[j], description_placeholder, gene_names[gene_indices])
      writeLines(paste(line, collapse = "\t"), con = gmt_con)
    }
  }
  
  close(gmt_con)
}

cat("\n========================================\n")
cat("SUCCESS!\n")
cat("Generated", num_random_sets, "randomized GMT files\n")
cat("Output location:", output_dir, "\n")
cat("========================================\n")