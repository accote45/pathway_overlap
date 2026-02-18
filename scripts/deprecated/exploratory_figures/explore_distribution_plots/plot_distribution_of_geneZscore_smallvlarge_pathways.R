## plot distribution of Z scores for small vs large random pathways

library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(GSA)

# load data
scores <- read.table('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_random/keeppathsize/msigdbgenes/cad/cad_set_random129.gsa.genes.out', header = TRUE)

# read in pathway assignments
dat <- GSA.read.gmt('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/randomized_gene_sets/random_keeppathsize/GeneSet.random129.gmt')

path.list <- dat$genesets
names(path.list) <- dat$geneset.names

# Create gene-pathway mapping
gene_pathway <- rbindlist(lapply(names(path.list), function(pathway_name) {
    data.table(GENE = path.list[[pathway_name]], PATHWAY = pathway_name)
}))

# Merge gene scores with pathway assignments
gene_scores_with_pathways <- merge(
    scores, 
    gene_pathway, 
    by = "GENE", 
    all.x = TRUE
)

# Calculate pathway-level statistics
pathway_stats <- gene_scores_with_pathways %>%
    filter(!is.na(PATHWAY)) %>%  # Remove genes not in any pathway
    group_by(PATHWAY) %>%
    summarise(
        n_genes = n(),
        mean_z = mean(ZSTAT, na.rm = TRUE),
        var_z = var(ZSTAT, na.rm = TRUE),
        sd_z = sd(ZSTAT, na.rm = TRUE),
        median_z = median(ZSTAT, na.rm = TRUE),
        min_z = min(ZSTAT, na.rm = TRUE),
        max_z = max(ZSTAT, na.rm = TRUE),
        .groups = 'drop'
    )

# Add pathway size categories
pathway_stats <- pathway_stats %>%
    mutate(
        pathway_size = case_when(
            n_genes < 50 ~ 'small',
            n_genes >= 50 & n_genes < 200 ~ 'medium', 
            n_genes >= 200 ~ 'large'
        ),
        pathway_size = factor(pathway_size, levels = c('small', 'medium', 'large'))
    )

# Filter for small and large pathways only (if desired)
pathway_stats_filtered <- pathway_stats %>%
    filter(pathway_size %in% c('small', 'large'))

# View the results
head(pathway_stats_filtered)

# Summary by pathway size
summary_by_size <- pathway_stats_filtered %>%
    group_by(pathway_size) %>%
    summarise(
        n_pathways = n(),
        mean_of_means = mean(mean_z, na.rm = TRUE),
        mean_of_vars = mean(var_z, na.rm = TRUE),
        sd_of_means = sd(mean_z, na.rm = TRUE),
        sd_of_vars = sd(var_z, na.rm = TRUE),
        .groups = 'drop'
    )


# Plot distributions
p1 <- ggplot(pathway_stats_filtered, aes(x = mean_z, fill = pathway_size)) +
    geom_density(alpha = 0.7) +
    labs(
        title = "Distribution of Pathway Mean Z-scores",
        x = "Mean Z-score per Pathway",
        y = "Density",
        fill = "Pathway Size"
    ) +
    theme_minimal()

p2 <- ggplot(pathway_stats_filtered, aes(x = var_z, fill = pathway_size)) +
    geom_density(alpha = 0.7) +
    labs(
        title = "Distribution of Pathway Z-score Variance", 
        x = "Variance of Z-scores per Pathway",
        y = "Density",
        fill = "Pathway Size"
    ) +
    theme_minimal()

p3 <- ggplot(pathway_stats_filtered, aes(x = pathway_size, y = mean_z, fill = pathway_size)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(
        title = "Pathway Mean Z-scores by Size",
        x = "Pathway Size",
        y = "Mean Z-score per Pathway"
    ) +
    theme_minimal() +
    theme(legend.position = "none")

p4 <- ggplot(pathway_stats_filtered, aes(x = pathway_size, y = var_z, fill = pathway_size)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(
        title = "Pathway Z-score Variance by Size",
        x = "Pathway Size", 
        y = "Variance of Z-scores per Pathway"
    ) +
    theme_minimal() +
    theme(legend.position = "none")

# Display plots
pdf('test.pdf')
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

# Statistical tests
wilcox_means <- wilcox.test(mean_z ~ pathway_size, data = pathway_stats_filtered)
wilcox_vars <- wilcox.test(var_z ~ pathway_size, data = pathway_stats_filtered)

cat("\nWilcoxon test for difference in mean Z-scores:\n")
cat("p-value:", format(wilcox_means$p.value, scientific = TRUE), "\n")

cat("\nWilcoxon test for difference in Z-score variances:\n") 
cat("p-value:", format(wilcox_vars$p.value, scientific = TRUE), "\n")










