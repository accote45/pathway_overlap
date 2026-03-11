library(tidyverse)
library(readxl)
library(ggplot2)

# Define tools and their corresponding Excel files
tools <- c("magma", "prset", "gsamixer", "pascal")
base_path <- "/sc/arion/projects/psychgen/cotea02_prset/pathway_overlap"

# Define excluded traits per tool
excluded_traits <- list(
  prset = c("SCZ", "AD", "IBD")  # Excludes SCZ, AD, IBD for PRSet
)

# Read data from all tools
all_tool_data <- map(tools, function(tool) {
  excel_path <- file.path(base_path, paste0("combined_empirical_pvalues_", tool, ".xlsx"))
  
  if (!file.exists(excel_path)) {
    cat(sprintf("Warning: File not found - %s\n", excel_path))
    return(NULL)
  }
  
  sheets <- excel_sheets(excel_path) %>%
    set_names() %>%
    map(~read_excel(excel_path, sheet = .x))
  
  # Filter out excluded traits for this tool
  if (tool %in% names(excluded_traits)) {
    exclude_list <- excluded_traits[[tool]]
    sheets <- sheets[!toupper(names(sheets)) %in% exclude_list]
    cat(sprintf("  Excluded traits for %s: %s\n", tool, paste(exclude_list, collapse = ", ")))
  }
  
  return(list(tool = tool, data = sheets))
}) %>%
  set_names(tools) %>%
  compact()  # Remove NULL entries for missing files

# Calculate correlations for all tools
summary_long <- imap_dfr(all_tool_data, function(tool_info, tool_name) {
  cat(sprintf("Processing %s...\n", tool_name))
  
  imap_dfr(tool_info$data, function(dat, sheet_name) {
    # Check for required columns
    if (!all(c("ngenes", "original_rank", "gsr_rank") %in% names(dat))) {
      cat(sprintf("  Warning: Missing columns in %s - %s\n", tool_name, sheet_name))
      return(NULL)
    }
    
    # Helper function to calculate correlation
    get_corr <- function(data, n, rank_col) {
      # Ensure we don't exceed available rows
      n_available <- min(n, nrow(data))
      if (n_available < 10) {
        return(tibble(rho = NA, p_value = NA, n_used = n_available))
      }
      
      test <- cor.test(data[1:n_available,]$ngenes, 
                      data[[rank_col]][1:n_available], 
                      method = 'spearman')
      tibble(rho = test$estimate, p_value = test$p.value, n_used = n_available)
    }
    
    bind_rows(
      get_corr(dat, 100, "original_rank") %>% mutate(n_pathways = 100, rank_type = "original"),
      get_corr(dat, 100, "gsr_rank") %>% mutate(n_pathways = 100, rank_type = "gsr"),
      get_corr(dat, 500, "original_rank") %>% mutate(n_pathways = 500, rank_type = "original"),
      get_corr(dat, 500, "gsr_rank") %>% mutate(n_pathways = 500, rank_type = "gsr")
    ) %>%
      mutate(trait = sheet_name, tool = tool_name, .before = 1)
  })
})

# Save combined results
write.table(summary_long, "pathsize_rank_correlation_results_all_tools.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Print summary
cat("\n=== Data Summary ===\n")
summary_long %>%
  filter(!is.na(rho)) %>%
  count(tool, rank_type, n_pathways) %>%
  print()

# Add significance labels
summary_long_with_sig <- summary_long %>%
  mutate(
    significant = p_value < 0.05,
    sig_label = case_when(
      is.na(p_value) ~ "NA",
      p_value < 0.001 ~ "p < 0.001",
      p_value < 0.01 ~ "p < 0.01",
      p_value < 0.05 ~ "p < 0.05",
      TRUE ~ "n.s."
    ),
    # Clean tool names for plotting
    tool_label = factor(tool, 
                       levels = c("magma", "prset", "gsamixer", "pascal"),
                       labels = c("MAGMA", "PRSet", "GSA-MiXeR", "PASCAL")),
    # Set rho to NA for non-significant results to make them gray
    rho_colored = ifelse(significant, rho, NA)
  )

# Create comprehensive heatmap comparing all tools
p_heatmap <- ggplot(summary_long_with_sig %>% filter(!is.na(rho)), 
                    aes(x = as.factor(n_pathways), y = trait, fill = rho_colored)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.5) +
  scale_fill_gradient2(low = "#0571b0", mid = "white", high = "#ca0020",
                       midpoint = 0, name = "Spearman's ρ\n(p < 0.05)",
                       limits = c(-1, 1),
                       na.value = "gray90") +
  facet_grid(tool_label ~ rank_type, 
             labeller = labeller(rank_type = c("original" = "Original Ranking", 
                                              "gsr" = "GSR Ranking"))) +
  labs(
    title = "Pathway Size-Rank Correlations Across All Enrichment Tools",
    subtitle = "Only significant correlations (p < 0.05) are colored; non-significant cells are gray",
    x = "Number of Top Pathways",
    y = "Trait"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "right",
    panel.spacing = unit(0.5, "lines")
  )

ggsave("pathsize_correlation_heatmap_all_tools.png", p_heatmap, 
       width = 12, height = 14, dpi = 300)

# Create scatterplot comparing GSR vs Original correlations
comparison_data <- summary_long_with_sig %>%
  filter(!is.na(rho)) %>%
  select(trait, tool, tool_label, n_pathways, rank_type, rho) %>%
  pivot_wider(names_from = rank_type, values_from = rho) %>%
  mutate(
    rho_change = gsr - original
  )

# Create enhanced scatterplot with interpretation regions
p_scatter <- ggplot(comparison_data, aes(x = original, y = gsr)) +
  # Add shaded region for bias reduction
  annotate("rect", xmin = -1, xmax = 0, ymin = -1, ymax = 0, 
           fill = "lightblue", alpha = 0.1) +
  annotate("text", x = -0.7, y = -0.95, 
           label = "Negative bias\n(larger pathways ranked higher)", 
           size = 3, color = "gray40", hjust = 0) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_point(aes(color = tool_label, shape = as.factor(n_pathways)), 
             size = 3, alpha = 0.7) +
  scale_color_brewer(palette = "Set1", name = "Tool") +
  scale_shape_manual(values = c("100" = 16, "500" = 17), 
                     name = "Top N Pathways") +
  labs(
    title = "GSR vs Original Pathway Size-Rank Correlations",
    subtitle = "Negative ρ = pathway size bias (larger pathways ranked higher); Points above diagonal = GSR reduces bias",
    x = "Original Ranking Correlation (ρ)",
    y = "GSR Ranking Correlation (ρ)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  coord_fixed(ratio = 1, xlim = c(-1, 1), ylim = c(-1, 1))

ggsave("pathsize_correlation_scatter_gsr_vs_original.png", p_scatter, 
       width = 10, height = 8, dpi = 300)

cat("\n=== Plots saved ===\n")
cat("1. pathsize_correlation_heatmap_all_tools.png/pdf\n")
cat("2. pathsize_correlation_scatter_gsr_vs_original.png/pdf\n")
cat("3. pathsize_rank_correlation_results_all_tools.txt\n")



