library(jsonlite)
library(tidyverse)
library(data.table)
library(patchwork)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]
tool_base <- args[2]  # e.g., magma or prset
birewire_results <- args[3]
keeppathsize_results <- args[4]

cat("Processing OpenTargets comparison for", trait, "using", tool_base, "\n")
cat("Birewire file:", birewire_results, "\n")
cat("Keeppathsize file:", keeppathsize_results, "\n")

# Define trait mapping
trait_mapping <- list(
  "t2d" = "MONDO_0005148",
  "cad" = "EFO_0001645", 
  "bmi" = "EFO_0004340",
  "ad" = "MONDO_0004975",
  "mdd" = "MONDO_0002050",
  "scz" = "MONDO_0005090",
  "ibd" = "EFO_0000555",
  "breast" = "MONDO_0007254",
  "HDL_cholesterol" = "EFO_0004612",
  "Lymphocyte_count" = "EFO_0004587",
  "Platelet_crit" = "EFO_0007985",
  "Alkaline_phosphatase" = "EFO_0004533"
)

# Read birewire and keeppathsize results
birewire_data <- fread(birewire_results)
keeppathsize_data <- fread(keeppathsize_results)

cat("Read", nrow(birewire_data), "pathways from birewire results\n")
cat("Read", nrow(keeppathsize_data), "pathways from keeppathsize results\n")

# Determine column names based on tool type
if (tool_base == "magma") {
  pathway_col <- "FULL_NAME"
  pval_col <- "empirical_pval"
} else if (tool_base == "prset") {
  pathway_col <- "Set"
  pval_col <- "empirical_pval"
} else {
  pathway_col <- "pathway"
  pval_col <- "empirical_pval"
}

# Find common pathways
common_pathways <- intersect(
  birewire_data[[pathway_col]],
  keeppathsize_data[[pathway_col]]
)

cat("Found", length(common_pathways), "common pathways between randomization methods\n")

# Create comparison table
results <- birewire_data %>%
  filter(!!sym(pathway_col) %in% common_pathways) %>%
  select(pathway = !!sym(pathway_col), birewire_pval = !!sym(pval_col)) %>%
  left_join(
    keeppathsize_data %>%
    filter(!!sym(pathway_col) %in% common_pathways) %>%
    select(pathway = !!sym(pathway_col), keeppathsize_pval = !!sym(pval_col)),
    by = "pathway"
  ) %>%
  mutate(
    trait = trait,
    tool = tool_base,
    log10_ratio = log10(birewire_pval / keeppathsize_pval),
    difference = birewire_pval - keeppathsize_pval,
    significant_diff = abs(log10_ratio) > 1,  # >10x difference
    significant_any = birewire_pval < 0.05 | keeppathsize_pval < 0.05,
    discordant = (birewire_pval < 0.05 & keeppathsize_pval >= 0.05) | 
                (birewire_pval >= 0.05 & keeppathsize_pval < 0.05)
  ) %>%
  arrange(pmin(birewire_pval, keeppathsize_pval))

# Write results table
output_file <- paste0(trait, "_", tool_base, "_opentargets_comparison.tsv")
fwrite(results, output_file, sep = "\t")

cat("Results written to", output_file, "\n")

# Create visualization
p1 <- ggplot(results, aes(x = birewire_pval, y = keeppathsize_pval)) +
  geom_point(aes(color = significant_diff), alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = paste("Comparison of Randomization Methods for", trait),
       x = "BireWire Empirical P-value",
       y = "KeepPathSize Empirical P-value") +
  theme_bw() +
  annotate("rect", xmin = 0, xmax = 0.05, ymin = 0, ymax = 0.05, 
           fill = "lightgreen", alpha = 0.2) +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "red") + 
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "red")

# Save visualization
pdf_file <- paste0(trait, "_", tool_base, "_opentargets_comparison.pdf")
ggsave(pdf_file, p1, width = 8, height = 6)

cat("Visualization saved to", pdf_file, "\n")

# Print summary statistics
cat("Summary statistics:\n")
cat("Pathways with >10x difference between methods:", sum(results$significant_diff), "\n")
cat("Pathways significant in any method:", sum(results$significant_any), "\n")
cat("Pathways with discordant significance:", sum(results$discordant), "\n")
cat("Correlation between methods:", cor(results$birewire_pval, results$keeppathsize_pval, method="spearman"), "\n")

