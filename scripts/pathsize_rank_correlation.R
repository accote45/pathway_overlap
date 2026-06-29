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
    if (!all(c("ngenes", "pathway_name", "original_rank", "gsr_rank") %in% names(dat))) {
      cat(sprintf("  Warning: Missing columns in %s - %s\n", tool_name, sheet_name))
      return(NULL)
    }
    
    # For gsamixer: filter to union of top-N pathways from other tools
    # For all others: take top N rows directly
    get_subset <- function(data, n) {
      if (tool_name == "gsamixer") {
        union_paths <- character(0)
        for (ot in c("magma", "prset", "pascal")) {
          if (!ot %in% names(all_tool_data)) next
          other_dat <- all_tool_data[[ot]]$data[[sheet_name]]
          if (is.null(other_dat)) next
          union_paths <- union(union_paths, head(other_dat, n)$pathway_name)
        }
        data[data$pathway_name %in% union_paths, ]
      } else {
        data[1:min(n, nrow(data)), ]
      }
    }
    
    get_corr <- function(data, n, rank_col) {
      subset_data <- get_subset(data, n)
      if (nrow(subset_data) < 10) {
        return(tibble(rho = NA, p_value = NA, n_used = nrow(subset_data)))
      }
      test <- cor.test(subset_data$ngenes, subset_data[[rank_col]], method = 'spearman')
      tibble(rho = test$estimate, p_value = test$p.value, n_used = nrow(subset_data))
    }
    
    bind_rows(
      get_corr(dat, 100, "original_rank") %>% mutate(n_pathways = 100, rank_type = "original"),
      get_corr(dat, 100, "gsr_rank")      %>% mutate(n_pathways = 100, rank_type = "gsr"),
      get_corr(dat, 500, "original_rank") %>% mutate(n_pathways = 500, rank_type = "original"),
      get_corr(dat, 500, "gsr_rank")      %>% mutate(n_pathways = 500, rank_type = "gsr")
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
                   levels = c("magma", "pascal", "prset", "gsamixer"),
                   labels = c("MAGMA", "PASCAL", "PRSet", "GSA-MiXeR")),
    # Set rho to NA for non-significant results to make them gray
    rho_colored = ifelse(significant, rho, NA)
  )

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
    subtitle = "Negative \u03c1 = pathway size bias (larger pathways ranked higher); Points above diagonal = GSR reduces bias",
    x = "Original Ranking Correlation (\u03c1)",
    y = "GSR Ranking Correlation (\u03c1)"
  ) +
  coord_fixed(ratio = 1, xlim = c(-1, 1), ylim = c(-1, 1)) +
  theme_bw(base_size = 11) +
  theme(
    plot.title        = element_text(face = "bold", color = "black", size = 13),
    plot.subtitle     = element_text(size = 9, color = "gray40"),
    plot.margin       = margin(8, 10, 8, 8),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.text         = element_text(color = "black"),
    axis.title        = element_text(color = "black", size=11),
    axis.ticks        = element_line(color = "black"),
    axis.line         = element_line(color = "black"),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    legend.title      = element_text(color = "black"),
    legend.text       = element_text(color = "black"),
    legend.position   = "right"
  )

ggsave("pathsize_correlation_scatter_gsr_vs_original.png", p_scatter,
       width = 10, height = 8, dpi = 300)

# --- Dumbbell Plot: connect Original → GSR, faceted by tool × n_pathways ---

dumbbell_points <- summary_long_with_sig %>%
  filter(!is.na(rho)) %>%
  mutate(rank_label = factor(rank_type, levels = c("original", "gsr"),
                             labels = c("Original", "GSR")))

dumbbell_segments <- summary_long_with_sig %>%
  filter(!is.na(rho)) %>%
  select(trait, tool, tool_label, n_pathways, rank_type, rho, significant) %>%
  pivot_wider(names_from = rank_type, values_from = c(rho, significant)) %>%
  rename(original = rho_original, gsr = rho_gsr) %>%
  filter(!is.na(original) & !is.na(gsr)) %>%
  mutate(significant = significant_original | significant_gsr)

p_dumbbell <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_segment(data = dumbbell_segments,
               aes(x = original, xend = gsr, y = trait, yend = trait,
                   alpha = significant),
               color = "gray60", linewidth = 0.6) +
  geom_point(data = dumbbell_points,
             aes(x = rho, y = trait, color = rank_label, alpha = significant),
             size = 3) +
  scale_color_manual(
    values = c("Original" = "#377EB8", "GSR" = "#E41A1C"),
    name = "Rank Type"
  ) +
  scale_alpha_manual(
    values = c("TRUE" = 1, "FALSE" = 0.25),
    name = "p<0.05",
    labels = c("TRUE" = "Yes", "FALSE" = "No")
  ) +
  facet_grid(tool_label ~ n_pathways,
             scales = "free_y",
             labeller = labeller(n_pathways = c("100" = "Top 100", "500" = "Top 500"))) +
  scale_y_discrete(limits = rev) +
  labs(
    title = "Pathway Size-Rank Correlations by Trait and Tool",
    subtitle = "Lines connect Original \u2192 GSR; Blue = Original, Red = GSR; Faded = n.s. (p \u2265 0.05)",
    x = "Spearman's \u03c1",
    y = "Trait"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title        = element_text(face = "bold", color = "black", size = 13),
    plot.subtitle     = element_text(size = 9, color = "black"),
    plot.margin       = margin(8, 10, 8, 8),
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold", color = "black", size = 10),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.text         = element_text(color = "black", size = 10),
    axis.title        = element_text(color = "black", size = 11),
    axis.ticks        = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = unit(2, "pt"),
    axis.line         = element_line(color = "black"),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    legend.title      = element_text(color = "black"),
    legend.text       = element_text(color = "black"),
    legend.position   = "right",
    panel.spacing     = unit(0.8, "lines")
  )

ggsave("pathsize_correlation_dotplot_by_tool.png", p_dumbbell,
       width = 12, height = 14, dpi = 300)

cat("\n=== Plots saved ===\n")
cat("2. pathsize_correlation_scatter_gsr_vs_original.png/pdf\n")
cat("3. pathsize_rank_correlation_results_all_tools.txt\n")
cat("4. pathsize_correlation_dotplot_by_tool.png\n")



