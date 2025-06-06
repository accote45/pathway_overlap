# Plot mean(BETA) for MAGMA gene set enrichment for random pathways of different sizes (keeppath + birewire randomization)

library(tidyverse)
library(data.table)
library(ggplot2)

paths <- readRDS('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/randomized_gene_sets/random_birewire/all_pathways.rds')

# read in all MAGMA results
files <- list.files('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_random/birewire/msigdbgenes/bmi', pattern = 'gsa.out', full.names = TRUE)
magma_results_birewire <- lapply(files, fread)

files <- list.files('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_random/keeppathsize/msigdbgenes/bmi', pattern = 'gsa.out', full.names = TRUE)
magma_results_keeppath <- lapply(files, fread)

magma_combined <- do.call(rbind, magma_results_birewire)

avg_beta <- magma_combined %>%
  group_by(NGENES) %>%
  summarise(mean_beta = mean(BETA, na.rm = TRUE)) %>%
  ungroup()


# random subsample of pathways to plot
set.seed(123)
random_paths <- sample(nrow(magma_combined), 100000)
temp <- magma_combined[random_paths, ]

# Create mean BETA values for specific pathway sizes (100, 200, 300, etc.)
target_sizes <- seq(100, max(magma_combined$NGENES, na.rm = TRUE), by = 100)
mean_points <- magma_combined %>%
  filter(NGENES %in% target_sizes) %>%
  group_by(NGENES) %>%
  summarise(mean_beta = mean(BETA, na.rm = TRUE), .groups = 'drop')

pdf('test.pdf')
ggplot(temp, aes(x = NGENES, y = BETA)) +
  geom_point(alpha = 0.5) +
  geom_point(data = mean_points, 
             aes(x = NGENES, y = mean_beta), 
             color = "red", 
             shape = 18,  # Diamond shape
             size = 4) +
  labs(x = 'Number of Genes in Pathway',
       y = 'MAGMA gene set enrichment - BETA',
       subtitle = 'Red diamonds show mean BETA for pathways of size X') +
  theme_minimal()
dev.off()

pdf('test.pdf')
ggplot(temp, aes(x = NGENES, y = -log10(P))) +
  geom_point(alpha = 0.5) +
  labs(x = 'Number of Genes in Pathway',
       y = 'MAGMA gene set enrichment - -log10(P)')+
  theme_minimal()
dev.off()






















