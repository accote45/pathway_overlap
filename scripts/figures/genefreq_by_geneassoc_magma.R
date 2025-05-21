# Load required libraries
library(tidyverse)
library(ggplot2)
library(data.table)
library(GSA)
library(ggpubr)

# Function to create plot for a specific trait
create_gene_pathway_plot <- function(trait_name) {
  # Read the gene set data (this stays constant)
  dat <- GSA.read.gmt('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/pathway_db/msigdb/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt')
  path.list <- dat$genesets
  names(path.list) <- dat$geneset.names
  
  # Convert the path list to a data table
  genes_long <- rbindlist(lapply(names(path.list), function(name) {
    data.table(value = path.list[[name]], name = name)
  }))
  
  # Calculate gene frequency
  gene_freq <- as.data.frame(table(genes_long$value))
  gene_freq$GENE <- gene_freq$Var1
  
  # Read trait-specific gene-level association
  trait_file_path <- file.path('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_real', 
                              trait_name, 
                              paste0(trait_name, "_real_set.gsa.genes.out"))
  
  # Check if file exists
  if (!file.exists(trait_file_path)) {
    warning(paste("File not found:", trait_file_path))
    return(NULL)
  }
  
  # Read gene-level association data
  gene <- read.table(trait_file_path, header=TRUE)
  
  # Merge gene frequency with gene-level data
  master <- merge(gene_freq, gene, by="GENE")
  
  # Create output directory if it doesn't exist
  output_dir <- file.path('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/pathway_plots')
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Create plot
  output_file <- file.path(output_dir, paste0(trait_name, "_pathway_gene_plot.pdf"))
  pdf(output_file)
  p <- ggplot(master, aes(x=ZSTAT, y=Freq)) + 
    geom_point(alpha=0.5) + 
    theme_classic() + 
    ylab("Number of pathways per gene") + 
    xlab(paste0(trait_name, " MAGMA gene-level z statistic")) +
    ggtitle(paste0("Gene pathway association for ", trait_name)) +
    # Add correlation statistics
    stat_cor(
      method = "pearson",
      label.x.npc = 0.85,  # Position at right side
      label.y.npc = 0.3,    # Position at top
      size = 5,             # Text size
      cor.coef.name = "r",  # Label for coefficient
      p.accuracy = 0.001,   # P-value decimal places
      r.accuracy = 0.001,
      label.sep = "\n"      # Correlation decimal places
    )
  print(p)
  dev.off()
  
  cat("Created plot for", trait_name, "at", output_file, "\n")
  return(p)
}

# List all trait folders
trait_dir <- '/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_real/'
trait_folders <- list.dirs(trait_dir, full.names = FALSE, recursive = FALSE)

# Filter out any empty or invalid folder names
trait_folders <- trait_folders[trait_folders != ""]

# Create plots for all traits
cat("Found", length(trait_folders), "traits to process\n")

# Create a list to store all plots
all_plots <- list()

# Process each trait
for (trait in trait_folders) {
  cat("Processing", trait, "...\n")
  plot_result <- create_gene_pathway_plot(trait)
  
  if (!is.null(plot_result)) {
    all_plots[[trait]] <- plot_result
  }
}



# plot pathsizexFPR corr by genefreqxgeneassoc corr

library(ggplot2)
library(tidyverse)
library(ggrepel)

dat <- read.table('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/explore/pathsizecorr_genefreqcorr_compare.txt', header=TRUE)

cor.test(temp$pathsize_fpr_corr,temp$genefreq_zstat_corr)

pdf('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/explore/pathsize_fpr_corr_genefreq_zstat_corr.pdf')
ggplot(temp, aes(x=pathsize_fpr_corr, y=genefreq_zstat_corr,label=trait)) +
  geom_point() +
  theme_classic() +
  geom_text_repel(size=3) +
  xlab('Pathway size x FPR correlation') +
  ylab('Gene frequency x z-statistic correlation')
dev.off()