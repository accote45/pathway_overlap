## plot effect size of MAGMA enrichment for random pathways of size 10-20 vs size 190-200

library(tidyverse)
library(data.table)

# load data

read_files <- function(path) {
  setwd(path)
  files <- list.files(pattern="*gsa.out")
  ldf <- lapply(files, function(f) tryCatch(read.table(f, header=T), error=function(e) NULL))
  ldf <- ldf[!sapply(ldf, is.null)]  # Remove any NULL entries
  return(ldf)
}

# Load data from both folders
birewire_data <- read_files("/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_random/birewire/msigdbgenes/mdd")
keeppathsize_data <- read_files("sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_random/keeppathsize/msigdbgenes/mdd")

# Process the data to extract effect sizes by pathway size
process_magma_results <- function(magma_list, randomization_type) {
  
  # Combine all results into one data frame
  all_results <- rbindlist(magma_list, idcol = "random_set")
  
  # Filter for the size ranges we're interested in
  size_filtered <- all_results[
    (NGENES >= 10 & NGENES <= 20) |    # Small pathways
    (NGENES >= 190 & NGENES <= 200)   # Large pathways
  ]
  
  # Create size categories and randomization type
  size_filtered[, size_category := ifelse(NGENES <= 20, "Small (10-20 genes)", "Large (190-200 genes)")]
  size_filtered[, randomization_method := randomization_type]
  
  return(size_filtered)
}

# Process both datasets
birewire_processed <- process_magma_results(birewire_data, "Keep gene frequency and pathway size")
keeppathsize_processed <- process_magma_results(keeppathsize_data, "Keep pathway size")

# Combine both datasets
magma_combined <- rbind(birewire_processed, keeppathsize_processed)

# Summary statistics for both methods
summary_stats_combined <- magma_combined[, .(
  n_pathways = .N,
  mean_beta = mean(BETA, na.rm = TRUE),
  median_beta = median(BETA, na.rm = TRUE),
  sd_beta = sd(BETA, na.rm = TRUE),
  mean_pvalue = mean(P, na.rm = TRUE),
  prop_significant = mean(P < 0.05, na.rm = TRUE)
), by = .(size_category, randomization_method)]

print("Summary statistics by pathway size and randomization method:")
print(summary_stats_combined)

# Define the color palette
color_palette <- c("Keep gene frequency and pathway size" = "#FF8C00",  # Dark orange
                   "Keep pathway size" = "#228B22")                      # Forest green

# Create combined effect size comparison plot with means
p1_combined <- ggplot(magma_combined, aes(x = size_category, y = BETA, fill = randomization_method)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3, position = position_dodge(width = 0.8)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "darkred",
               position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "MAGMA Effect Sizes: Small vs Large Random Pathways",
    subtitle = "Comparison of different randomization methods",
    x = "Pathway Size Category",
    y = "Beta (Effect Size)",
    fill = "Randomization Method",
    caption = "Black dashed line = no effect, Dark red diamonds = means"
  ) +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Create combined density plot
p2_combined <- ggplot(magma_combined, aes(x = BETA, fill = randomization_method)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~size_category, scales = "free_y") +
  labs(
    title = "Distribution of Effect Sizes by Randomization Method",
    x = "Beta (Effect Size)",
    y = "Density",
    fill = "Randomization Method"
  ) +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Create combined p-value comparison with means
p3_combined <- ggplot(magma_combined, aes(x = size_category, y = -log10(P), fill = randomization_method)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3, position = position_dodge(width = 0.8)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "darkred",
               position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    title = "MAGMA Significance: Small vs Large Random Pathways",
    subtitle = "Comparison of different randomization methods",
    x = "Pathway Size Category",
    y = "-log10(P-value)",
    fill = "Randomization Method",
    caption = "Black dashed line = p < 0.05 threshold, Dark red diamonds = means"
  ) +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Create FPR comparison plot
p4_fpr <- ggplot(prop_analysis, aes(x = size_category, y = proportion_significant, 
                                   fill = randomization_method)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.8) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  labs(
    title = "False Positive Rate by Randomization Method",
    x = "Pathway Size Category",
    y = "Proportion Significant (FPR)",
    fill = "Randomization Method",
    caption = "Black dashed line = expected 5% FPR"
  ) +
  scale_fill_manual(values = color_palette) +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Statistical tests for differences
# Test effect size differences between methods for large pathways
large_pathway_test <- wilcox.test(
  BETA ~ randomization_method, 
  data = magma_combined[size_category == "Large (190-200 genes)"],
  alternative = "two.sided"
)

# Test effect size differences between methods for small pathways
small_pathway_test <- wilcox.test(
  BETA ~ randomization_method, 
  data = magma_combined[size_category == "Small (10-20 genes)"],
  alternative = "two.sided"
)

# Test effect size differences for small vs large pathways (same method)
large_vs_small_test <- wilcox.test(
  BETA ~ size_category, 
  data = magma_combined[randomization_method == "Keep pathway size"],
  alternative = "two.sided"
)

cat("\nWilcoxon test for randomization method differences:\n")
cat("Large pathways p-value:", format(large_pathway_test$p.value, scientific = TRUE), "\n")
cat("Small pathways p-value:", format(small_pathway_test$p.value, scientific = TRUE), "\n")

# Get proportion of significant pathways by method and size
prop_analysis <- magma_combined %>% 
  group_by(size_category, randomization_method) %>% 
  summarize(
    total_pathways = n(),
    significant_count = sum(P < 0.05, na.rm = TRUE),
    proportion_significant = mean(P < 0.05, na.rm = TRUE),
    .groups = 'drop'
  )

print("\nProportion of significant pathways by method and size:")
print(prop_analysis)


# Save all plots
pdf('combined_magma_analysis.pdf', width = 12, height = 8)
print(p1_combined)
print(p2_combined) 
print(p3_combined)
print(p4_fpr)
dev.off()

# Additional analysis: Effect size differences
effect_size_comparison <- magma_combined[, .(
  mean_effect = round(mean(BETA, na.rm = TRUE), 6),
  median_effect = round(median(BETA, na.rm = TRUE), 6),
  q25 = round(quantile(BETA, 0.25, na.rm = TRUE), 6),
  q75 = round(quantile(BETA, 0.75, na.rm = TRUE), 6)
), by = .(size_category, randomization_method)]

print("\nEffect Size Comparison:")
print(effect_size_comparison)

# Show the key finding: FPR differences between methods
fpr_summary <- prop_analysis %>%
  select(size_category, randomization_method, proportion_significant) %>%
  pivot_wider(names_from = randomization_method, values_from = proportion_significant) %>%
  mutate(
    fpr_difference = `Keep gene frequency and pathway size` - `Keep pathway size`,
    fold_change = `Keep gene frequency and pathway size` / `Keep pathway size`
  )

print("\nFPR Comparison Summary:")
print(fpr_summary)