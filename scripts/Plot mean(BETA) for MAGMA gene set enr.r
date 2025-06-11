# Plot mean(BETA) for MAGMA gene set enrichment for random pathways of different sizes (keeppath + birewire randomization)

library(tidyverse)
library(data.table)
library(ggplot2)

paths <- readRDS('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/randomized_gene_sets/random_birewire/all_pathways.rds')

# Function to read MAGMA results for all traits
read_all_traits <- function(base_path, method_name) {
  # Get all trait directories
  trait_dirs <- list.dirs(base_path, recursive = FALSE, full.names = FALSE)
  
  all_results <- list()
  
  for (trait in trait_dirs) {
    files <- list.files(base_path, pattern = 'gsa.out', full.names = TRUE)
    
    if (length(files) > 0) {
      trait_results <- lapply(files, fread)
      trait_combined <- do.call(rbind, trait_results) %>%
        mutate(method = method_name,
               trait = trait)
      all_results[[trait]] <- trait_combined
    }
  }
  
  return(do.call(rbind, all_results))
}

# Read in all MAGMA results for both methods
magma_results_birewire <- read_all_traits(
  '/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_random/birewire/msigdbgenes/'
)

magma_results_keeppath <- read_all_traits(
  '/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_random/keeppathsize/msigdbgenes/'
)

# Combine both methods
magma_all_methods <- rbind(magma_results_birewire, magma_results_keeppath)

# Check what traits we have
cat("Available traits:\n")
print(unique(magma_all_methods$trait))

# Create mean BETA values for specific pathway sizes (100, 200, 300, etc.)
target_sizes <- seq(100, 1000, by = 100)
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



# Create forest plot data for both methods
target_sizes <- seq(100, 1000, by = 100)

forest_dat_grouped <- magma_all_methods %>%
  filter(NGENES %in% target_sizes) %>%
  group_by(NGENES, method) %>%
  summarise(
    mean_beta = mean(BETA, na.rm = TRUE),
    sd_beta = sd(BETA, na.rm = TRUE),
    se_beta = sd(BETA, na.rm = TRUE) / sqrt(n()),  # Standard error
    n_pathways = n(),
    .groups = 'drop'
  )

# Create grouped forest plot
pdf('forest_plot_grouped_methods.pdf', width = 12, height = 8)

# Side-by-side with error bars (SD)
p1 <- ggplot(forest_dat_grouped, aes(x = factor(NGENES), y = mean_beta, color = method)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = mean_beta - sd_beta, ymax = mean_beta + sd_beta), 
                position = position_dodge(width = 0.5), width = 0.3) +
  labs(x = 'Number of Genes in Pathway',
       y = 'Mean BETA (± SD)',
       title = 'Forest Plot: Mean BETA by Pathway Size and Method',
       subtitle = 'Error bars show ± 1 SD',
       color = 'Randomization Method') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("birewire" = "#E69F00", "keeppathsize" = "#0072B2"))

# Faceted by method
p2 <- ggplot(forest_dat_grouped, aes(x = factor(NGENES), y = mean_beta)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_beta - sd_beta, ymax = mean_beta + sd_beta), width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_wrap(~method, ncol = 1) +
  labs(x = 'Number of Genes in Pathway',
       y = 'Mean BETA (± SD)',
       title = 'Forest Plot: Mean BETA by Pathway Size',
       subtitle = 'Comparison of randomization methods') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)
print(p2)
dev.off()


# Explore relationship between effect size and significance
pdf('effect_size_analysis.pdf', width = 12, height = 8)

# Plot 4: Beta vs P-value colored by pathway size
p4 <- ggplot(sample_n(magma_combined_birewire, 50000), 
             aes(x = BETA, y = -log10(P), color = NGENES)) +
  geom_point(alpha = 0.6) +
  scale_color_viridis_c(name = "Pathway\nSize") +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  labs(x = 'Effect Size (BETA',
       y = '-log10(P-value)',
       title = 'Effect Size vs Statistical Significance',
       subtitle = 'Colored by pathway size') +
  theme_minimal()

# Plot 5: Same effect sizes, different significance by pathway size
effect_bins <- magma_combined_birewire %>%
  mutate(
    beta_bin = round(BETA, 2),  # Bin similar effect sizes
    size_cat = case_when(
      NGENES < 50 ~ "Small (<50)",
      NGENES >= 50 & NGENES < 200 ~ "Medium (50-199)",
      NGENES >= 200 ~ "Large (≥200)"
    )
  ) %>%
  group_by(beta_bin, size_cat) %>%
  summarise(
    mean_p = mean(P, na.rm = TRUE),
    n_obs = n(),
    .groups = 'drop'
  )
p5 <- ggplot(effect_bins, aes(x = beta_bin, y = -log10(mean_p), color = size_cat)) +
  geom_point(aes(size = n_obs), alpha = 0.7) +
 # geom_smooth(method = "loess", se = FALSE) +
  labs(x = 'Effect Size (BETA)',
       y = 'Mean -log10(P-value)',
       title = 'Same Effect Sizes, Different P-values by Pathway Size',
       subtitle = 'Point size = number of observations',
       color = "Pathway Size") +
  theme_minimal()

print(p4)
print(p5)
dev.off()





























# Add this to your existing script after loading the data

# Calculate statistical power metrics for each pathway size
power_analysis <- magma_combined %>%
  group_by(NGENES) %>%
  summarise(
    n_pathways = n(),
    mean_beta = mean(BETA, na.rm = TRUE),
    se_beta = sd(BETA, na.rm = TRUE) / sqrt(n()),  # Standard error of mean
    mean_p = mean(P, na.rm = TRUE),
    prop_significant = mean(P < 0.05, na.rm = TRUE),  # FPR proxy
    prop_highly_sig = mean(P < 0.001, na.rm = TRUE),
    median_p = median(P, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(n_pathways >= 1000)  # Only include sizes with sufficient pathways

# Create power-related plots
pdf('power_analysis_plots.pdf', width = 12, height = 8)

# Plot 1: FPR vs Pathway Size
p1 <- ggplot(power_analysis, aes(x = NGENES, y = prop_significant)) +
  geom_point(aes(size = n_pathways), alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE, color = "red") +
  labs(x = 'Number of Genes in Pathway',
       y = 'Proportion of Pathways with P < 0.05',
       title = 'False Positive Rate vs Pathway Size',
       subtitle = 'Point size = number of pathways per size') +
  theme_minimal()

# Plot 2: Standard Error vs Pathway Size (Power relationship)
p2 <- ggplot(power_analysis, aes(x = NGENES, y = se_beta)) +
  geom_point(aes(size = n_pathways), alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  labs(x = 'Number of Genes in Pathway',
       y = 'Standard Error of Beta',
       title = 'Statistical Precision vs Pathway Size',
       subtitle = 'Lower SE = Higher Power') +
  theme_minimal()

# Plot 3: Distribution of P-values by pathway size categories
size_categories <- magma_combined %>%
  mutate(
    size_cat = case_when(
      NGENES < 50 ~ "Small (<50)",
      NGENES >= 50 & NGENES < 100 ~ "Medium (50-99)",
      NGENES >= 100 & NGENES < 200 ~ "Large (100-199)",
      NGENES >= 200 ~ "Very Large (≥200)"
    ),
    size_cat = factor(size_cat, levels = c("Small (<50)", "Medium (50-99)", 
                                           "Large (100-199)", "Very Large (≥200)"))
  )

p3 <- ggplot(size_categories, aes(x = P, fill = size_cat)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  facet_wrap(~size_cat, scales = "free_y") +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
  labs(x = 'P-value',
       y = 'Count',
       title = 'Distribution of P-values by Pathway Size',
       subtitle = 'Red line = α = 0.05') +
  theme_minimal() +
  theme(legend.position = "none")

print(p1)
print(p2)
print(p3)
dev.off()







# Explore relationship between effect size and significance
pdf('effect_size_analysis.pdf', width = 12, height = 8)

# Plot 4: Beta vs P-value colored by pathway size
p4 <- ggplot(sample_n(magma_combined, 50000), 
             aes(x = abs(BETA), y = -log10(P), color = NGENES)) +
  geom_point(alpha = 0.6) +
  scale_color_viridis_c(name = "Pathway\nSize") +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  labs(x = '|Effect Size (BETA)|',
       y = '-log10(P-value)',
       title = 'Effect Size vs Statistical Significance',
       subtitle = 'Colored by pathway size') +
  theme_minimal()

# Plot 5: Same effect sizes, different significance by pathway size
effect_bins <- magma_combined %>%
  mutate(
    beta_bin = round(BETA, 2),  # Bin similar effect sizes
    size_cat = case_when(
      NGENES < 50 ~ "Small (<50)",
      NGENES >= 50 & NGENES < 200 ~ "Medium (50-199)",
      NGENES >= 200 ~ "Large (≥200)"
    )
  ) %>%
  filter(abs(beta_bin) < 1) %>%  # Focus on reasonable effect sizes
  group_by(beta_bin, size_cat) %>%
  summarise(
    mean_p = mean(P, na.rm = TRUE),
    n_obs = n(),
    .groups = 'drop'
  ) %>%
  filter(n_obs >= 5)  # Only bins with sufficient observations

p5 <- ggplot(effect_bins, aes(x = beta_bin, y = -log10(mean_p), color = size_cat)) +
  geom_point(aes(size = n_obs), alpha = 0.7) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(x = 'Effect Size (BETA)',
       y = 'Mean -log10(P-value)',
       title = 'Same Effect Sizes, Different P-values by Pathway Size',
       subtitle = 'Point size = number of observations',
       color = "Pathway Size") +
  theme_minimal()

print(p4)
print(p5)
dev.off()




# Better power analysis - look at the precision of individual pathway tests
pathway_precision <- magma_combined %>%
  mutate(
    t_stat = BETA / SE,    # Approximate t-statistic
    size_cat = case_when(
      NGENES < 50 ~ "Small (<50)",
      NGENES >= 50 & NGENES < 200 ~ "Medium (50-199)", 
      NGENES >= 200 ~ "Large (≥200)"
    )
  )

# Plot the relationship between pathway size and test precision
pdf('pathway_test_precision.pdf', width = 12, height = 8)

p_precision <- ggplot(sample_n(pathway_precision, 50000), 
                     aes(x = NGENES, y = abs(t_stat), color = size_cat)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(x = 'Number of Genes in Pathway',
       y = '|T-statistic| (Effect Size / SE)',
       title = 'Statistical Power: Larger Pathways Have Higher T-statistics',
       subtitle = 'Same effect size becomes more significant with more genes',
       color = "Pathway Size") +
  theme_minimal()

print(p_precision)
dev.off()

# Statistical test
cat("Correlation between pathway size and |t-statistic|:\n")
cor_test <- cor.test(pathway_precision$NGENES, abs(pathway_precision$t_stat))
cat("Correlation:", round(cor_test$estimate, 4), "\n")
cat("P-value:", format(cor_test$p.value, scientific = TRUE), "\n")

# For one-tailed tests, the relationship is:
# t = BETA / SE  (not abs(BETA) / SE)
# P = 1 - pnorm(t)  # One-tailed probability

pathway_precision_corrected <- magma_combined_birewire %>%
  mutate(
    # For one-tailed test
    t_stat = BETA / SE,  # Can be negative
    expected_p = 1 - pnorm(t_stat),  # One-tailed p-value
    size_cat = case_when(
      NGENES < 50 ~ "Small (<50)",
      NGENES >= 50 & NGENES < 200 ~ "Medium (50-199)", 
      NGENES >= 200 ~ "Large (≥200)"
    )
  )

# Plot corrected relationship
p_precision_corrected <- ggplot(sample_n(pathway_precision_corrected, 50000), 
                               aes(x = NGENES, y = t_stat, color = size_cat)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(x = 'Number of Genes in Pathway',
       y = 'T-statistic (BETA / SE)',
       title = 'T-statistics vs Pathway Size (One-tailed test)',
       subtitle = 'Positive t-stats more likely to be significant',
       color = "Pathway Size") +
  theme_minimal()

print(p_precision_corrected)
dev.off()

# Plot mean(BETA) for MAGMA gene set enrichment for random pathways of different sizes (keeppath + birewire randomization)

library(tidyverse)
library(data.table)
library(ggplot2)

paths <- readRDS('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/randomized_gene_sets/random_birewire/all_pathways.rds')

# Function to read MAGMA results for all traits
read_all_traits <- function(base_path, method_name) {
  # Get all trait directories
  trait_dirs <- list.dirs(base_path, recursive = FALSE, full.names = FALSE)
  
  all_results <- list()
  
  for (trait in trait_dirs) {
    trait_path <- file.path(base_path, trait)
    files <- list.files(trait_path, pattern = 'gsa.out', full.names = TRUE)
    
    if (length(files) > 0) {
      trait_results <- lapply(files, fread)
      trait_combined <- do.call(rbind, trait_results) %>%
        mutate(method = method_name,
               trait = trait)
      all_results[[trait]] <- trait_combined
    }
  }
  
  return(do.call(rbind, all_results))
}

# Read in all MAGMA results for both methods
magma_results_birewire <- read_all_traits(
  '/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_random/birewire/msigdbgenes/',
  'birewire'
)

magma_results_keeppath <- read_all_traits(
  '/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_random/keeppathsize/msigdbgenes/',
  'keeppathsize'
)

# Combine both methods
magma_all_methods <- rbind(magma_results_birewire, magma_results_keeppath)

# Check what traits we have
cat("Available traits:\n")
print(unique(magma_all_methods$trait))

# Create forest plot data for both methods and all traits
target_sizes <- seq(100, 1000, by = 100)

forest_dat_grouped <- magma_all_methods %>%
  filter(NGENES %in% target_sizes) %>%
  group_by(NGENES, method, trait) %>%
  summarise(
    mean_beta = mean(BETA, na.rm = TRUE),
    sd_beta = sd(BETA, na.rm = TRUE),
    se_beta = sd(BETA, na.rm = TRUE) / sqrt(n()),  # Standard error
    n_pathways = n(),
    .groups = 'drop'
  )

# Create grouped forest plot showing all traits
pdf('forest_plot_all_traits_grouped.pdf', width = 16, height = 12)

# Forest plot faceted by trait, grouped by method
p1 <- ggplot(forest_dat_grouped, aes(x = factor(NGENES), y = mean_beta, color = method)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = mean_beta - sd_beta, ymax = mean_beta + sd_beta), 
                position = position_dodge(width = 0.5), width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_wrap(~ trait, scales = "free_y", ncol = 4,
             labeller = labeller(trait = function(x) gsub("_", "\n", x))) +
  labs(x = 'Number of Genes in Pathway',
       y = 'Mean BETA (± SD)',
       title = 'Forest Plot: Mean BETA by Pathway Size, Method, and Trait',
       subtitle = 'Error bars show ± 1 SD',
       color = 'Randomization\nMethod') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 10),
        legend.position = "top") +
  scale_color_manual(values = c("birewire" = "#0072B2", "keeppathsize" = "#E69F00"))

print(p1)

# Alternative: Side-by-side comparison for each trait
p2 <- ggplot(forest_dat_grouped, aes(x = method, y = mean_beta, fill = method)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean_beta - sd_beta, ymax = mean_beta + sd_beta), 
                position = position_dodge(width = 0.8), width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_grid(trait ~ factor(NGENES), scales = "free_y",
             labeller = labeller(trait = function(x) gsub("_", "\n", x))) +
  labs(x = 'Randomization Method',
       y = 'Mean BETA (± SD)',
       title = 'Forest Plot: Method Comparison Across Traits and Pathway Sizes',
       fill = 'Randomization\nMethod') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 8),
        legend.position = "top") +
  scale_fill_manual(values = c("birewire" = "#0072B2", "keeppathsize" = "#E69F00"))

print(p2)

# Summary plot: Mean across all pathway sizes for each trait
summary_dat <- forest_dat_grouped %>%
  group_by(trait, method) %>%
  summarise(
    overall_mean_beta = mean(mean_beta, na.rm = TRUE),
    overall_sd = sqrt(mean(sd_beta^2, na.rm = TRUE)),  # Pooled SD
    .groups = 'drop'
  )

p3 <- ggplot(summary_dat, aes(x = trait, y = overall_mean_beta, fill = method)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = overall_mean_beta - overall_sd, 
                    ymax = overall_mean_beta + overall_sd), 
                position = position_dodge(width = 0.8), width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(x = 'GWAS Trait',
       y = 'Overall Mean BETA (± SD)',
       title = 'Summary: Mean BETA Across All Pathway Sizes by Trait',
       fill = 'Randomization\nMethod') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "top") +
  scale_fill_manual(values = c("birewire" = "#0072B2", "keeppathsize" = "#E69F00")) +
  scale_x_discrete(labels = function(x) gsub("_", "\n", x))

print(p3)

dev.off()

# Print summary statistics
cat("\nSummary by trait and method:\n")
summary_stats <- magma_all_methods %>%
  group_by(trait, method) %>%
  summarise(
    n_tests = n(),
    mean_beta = mean(BETA, na.rm = TRUE),
    mean_p = mean(P, na.rm = TRUE),
    prop_significant = mean(P < 0.05, na.rm = TRUE),
    .groups = 'drop'
  )

print(summary_stats)







