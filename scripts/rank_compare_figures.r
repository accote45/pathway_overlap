args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Args: tool_base birewire_results_pattern [top_ns_csv]")
}
tool_base <- args[[1]]
birewire_results_pattern <- args[[2]]  # e.g., "%s/%s_%s_birewire_empirical_pvalues.txt"
top_ns <- if (length(args) >= 3) as.integer(strsplit(args[[3]], ",")[[1]]) else c(50,100,250,500)

library(data.table)
library(dplyr)
library(ggplot2)

# Define all available traits
traits <- c("ad", "Alkaline_phosphatase", "breast", "bmi", "cad", "HDL_cholesterol", 
           "ibd", "Lipoprotein_A", "mdd", "Platelet_crit", "scz", "t2d")

# Function to process one trait
process_trait <- function(trait, tool_base, birewire_results_pattern) {
  # Construct the file path
  birewire_results_file <- sprintf(birewire_results_pattern, trait, trait, tool_base)
  
  # Check if file exists
  if (!file.exists(birewire_results_file)) {
    message("Skipping trait ", trait, ": file not found - ", birewire_results_file)
    return(NULL)
  }
  
  # Read BireWire results
  BW <- fread(birewire_results_file)
  pcol <- intersect(names(BW), c("pathway_name","name","pathway","Pathway"))
  if (length(pcol) == 0) {
    message("Skipping trait ", trait, ": Could not find pathway name column")
    return(NULL)
  }
  setnames(BW, pcol[1], "pathway")
  
  to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
  for (nm in intersect(c("empirical_pval","std_effect_size","p_value","beta_value"), names(BW))) {
    BW[[nm]] <- to_num(BW[[nm]])
  }
  
  # Create ranks
  BW_ord <- BW[order(empirical_pval, -std_effect_size)]
  BW_ord[, rank_birewire := seq_len(.N)]
  BW_ranks <- BW_ord[, .(pathway, rank_birewire)]
  
  PB_ord <- BW[order(p_value, -beta_value)]
  PB_ord[, rank_pvaluebeta := seq_len(.N)]
  PB_ranks <- PB_ord[, .(pathway, rank_pvaluebeta)]
  
  DR <- merge(BW_ranks, PB_ranks, by = "pathway")
  
  # Add delta_rank and trait/tool metadata
  DR[, delta_rank := rank_birewire - rank_pvaluebeta]
  DR[, trait := trait]
  DR[, tool_base := tool_base]
  
  # Remove base model
  DR <- DR[DR$pathway != "Base"]
  
  # Get pathway size
  size_candidates <- intersect(names(BW), c("ngenes","num_genes","n_genes","size","pathway_size","genes_in_set","n_genes_in_set","set_size","n"))
  if (length(size_candidates) > 0) {
    size_col <- size_candidates[1]
    sizes_dt <- unique(BW[, .(pathway, pathway_size = as.integer(get(size_col)))])
    DR <- merge(DR, sizes_dt, by = "pathway", all.x = TRUE)
  } else {
    DR[, pathway_size := NA_integer_]
  }
  
  return(DR)
}

# Process all traits
all_results <- list()
for (trait in traits) {
  result <- process_trait(trait, tool_base, "%s/%s_%s_birewire_empirical_pvalues.txt")
  if (!is.null(result)) {
    all_results[[trait]] <- result
  }
}

# Combine all results
if (length(all_results) == 0) {
  stop("No valid trait data found")
}

DR_all <- rbindlist(all_results, fill = TRUE)

# Create the grid plot with one subplot per trait
p_grid <- ggplot(DR_all, aes(x = rank_pvaluebeta, y = rank_birewire, color = pathway_size)) +
  geom_point(alpha = 0.75, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 1.2) +
  scale_color_viridis_c(option = "D", na.value = "grey80") +
  labs(x = "Rank (PvalueBeta baseline)", y = "Rank (BireWire)",
       color = "Pathway size\n(# genes)",
       title = sprintf("%s — Rank vs Rank colored by pathway size", tool_base)) +
  facet_wrap(~ trait, scales = "free", ncol = 4) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(size = 9))

# Save the grid plot
out_grid <- sprintf("%s_all_traits_rank_vs_rank_grid.pdf", tool_base)
ggsave(out_grid, p_grid, width = 16, height = 12, device = cairo_pdf)
message("Wrote: ", out_grid)

# Also calculate and print correlation summary for all traits
summ_all <- DR_all[, {
  if (.N >= 3) {
    ct <- suppressWarnings(cor.test(rank_pvaluebeta, rank_birewire, method = "kendall", exact = FALSE))
    list(kendall_tau = unname(ct$estimate), p_value = unname(ct$p.value), n = .N)
  } else {
    list(kendall_tau = NA_real_, p_value = NA_real_, n = .N)
  }
}, by = trait]

print("Kendall correlation summary by trait:")
options(scipen = 999)  # Disable scientific notation
print(summ_all, digits = 15)
options(scipen = 0)   # Reset to default