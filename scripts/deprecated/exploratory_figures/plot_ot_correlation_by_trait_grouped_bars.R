suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(stringr)
  library(scales)
})

# ---------- Configuration (edit as needed) ----------
INPUT_ROOT <- "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/opentargets_correlation"
TOOL_BASE  <- "magma"  # "magma" or "prset"
METHODS    <- c("BireWire_empP", "RawP")  # Which methods to plot
METRIC     <- "mean"    # "mean" or "density"
CORR_TYPES <- c("spearman", "kendall")    # Which correlation types to render
SUBSETS    <- c("All Pathways",
                "Top 50 Pathways", "Top 100 Pathways", "Top 250 Pathways", "Top 500 Pathways")
ABS_VALUES <- TRUE       # plot absolute correlations (|rho| / |tau|)
OUT_DIR    <- file.path(INPUT_ROOT, "ot_corr_grouped_by_trait")

# Method labels for legend
method_labels_map <- c(
  "BireWire_empP" = "Empirical p-value (account for gene frequency and pathway size)",
  "RawP"          = "Raw p-value"
)

# Optional trait label mapping
trait_labels_map <- c(
  "ibd"    = "IBD",
  "cad"    = "CAD",
  "breast" = "Breast cancer",
  "scz"    = "SCZ",
  "t2d"    = "T2D",
  "ad"     = "AD",
  "mdd"    = "MDD"
)
# ----------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(paste0("[plot-ot] ", sprintf(...), "\n"))

# Discover all summary CSVs
files <- list.files(INPUT_ROOT, pattern = "_rank_correlation_summary\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(files) == 0) {
  stop("No *_rank_correlation_summary.csv files found under: ", INPUT_ROOT)
}

# Parse trait and tool_base from filename: <trait>_<tool>_rank_correlation_summary.csv
parse_meta <- function(fp) {
  bn <- basename(fp)
  m <- str_match(bn, "^(.+?)_([A-Za-z0-9]+)_rank_correlation_summary\\.csv$")
  tibble(file = fp, trait = m[,2], tool_base = m[,3])
}
meta <- map_dfr(files, parse_meta) %>% filter(!is.na(trait), !is.na(tool_base))

# Keep only selected tool_base
meta <- meta %>% filter(tool_base == TOOL_BASE)
if (nrow(meta) == 0) stop("No files for tool_base = '", TOOL_BASE, "'.")

# Load all files
read_one <- function(row) {
  df <- suppressMessages(fread(row$file)) %>% as_tibble()
  # Ensure required columns exist; handle older/newer scripts
  needed <- c("method", "subset", "n_pathways")
  if (!all(needed %in% names(df))) {
    stop("Missing required columns in file: ", row$file)
  }
  df %>%
    mutate(trait = row$trait, tool_base = row$tool_base)
}
dat <- meta %>% pmap_dfr(read_one)

# Filter desired subsets and methods
dat <- dat %>%
  filter(subset %in% SUBSETS, method %in% METHODS)

if (nrow(dat) == 0) stop("No rows left after filtering by METHODS and SUBSETS.")

# Decide which columns contain correlation and p-value
get_cols <- function(metric = "mean", corr_type = "spearman") {
  stopifnot(metric %in% c("mean","density"))
  stopifnot(corr_type %in% c("spearman","kendall"))
  if (metric == "mean" && corr_type == "spearman") {
    list(corr = "rank_mean_score_correlation", pval = "rank_mean_score_pvalue",
         corr_label = "Spearman |rho| with OpenTargets mean score")
  } else if (metric == "density" && corr_type == "spearman") {
    list(corr = "rank_evidence_density_correlation", pval = "rank_evidence_density_pvalue",
         corr_label = "Spearman |rho| with OpenTargets evidence density")
  } else if (metric == "mean" && corr_type == "kendall") {
    list(corr = "kendall_mean_score_tau", pval = "kendall_mean_score_pvalue",
         corr_label = "Kendall |tau| with OpenTargets mean score")
  } else { # density + kendall
    list(corr = "kendall_evidence_density_tau", pval = "kendall_evidence_density_pvalue",
         corr_label = "Kendall |tau| with OpenTargets evidence density")
  }
}

# Stars for significance
p_to_stars <- function(p) {
  case_when(
    is.na(p)          ~ "",
    p < 1e-4          ~ "***",
    p < 1e-3          ~ "**",
    p < 5e-2          ~ "*",
    TRUE              ~ ""
  )
}

# Safe recode for methods (avoid length mismatch errors)
apply_method_labels <- function(x) {
  y <- recode(x, !!!method_labels_map)
  ifelse(is.na(y) | y == "", as.character(x), y)
}

# Trait label mapping with fallback
apply_trait_labels <- function(x) {
  y <- recode(x, !!!trait_labels_map)
  ifelse(is.na(y) | y == "", str_to_title(gsub("_", " ", x)), y)
}

# Make a nice slug
slugify <- function(s) {
  s %>%
    tolower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
}

# Plot one subset and correlation type
plot_one <- function(df, subset_name, corr_type) {
  cols <- get_cols(metric = METRIC, corr_type = corr_type)
  if (!all(c(cols$corr, cols$pval) %in% names(df))) {
    msg("Skipping: missing columns for %s / %s in data", subset_name, corr_type)
    return(invisible(NULL))
  }

  d <- df %>%
    filter(subset == subset_name) %>%
    transmute(
      trait,
      trait_label = apply_trait_labels(trait),
      method,
      method_label = apply_method_labels(method),
      n_pathways,
      corr = .data[[cols$corr]],
      pval = .data[[cols$pval]]
    ) %>%
    filter(!is.na(corr), !is.na(pval)) %>%
    mutate(
      plot_val = if (ABS_VALUES) abs(corr) else corr,
      stars = p_to_stars(pval)
    )

  if (nrow(d) == 0) {
    msg("No rows to plot for subset '%s' and corr_type '%s'", subset_name, corr_type)
    return(invisible(NULL))
  }

  # Order traits by provided mapping first, otherwise alphabetical
  trait_order <- unique(d$trait_label)
  # Keep mapping order for known labels
  known_order <- trait_labels_map[trait_labels_map %in% trait_order]
  other_order <- setdiff(sort(trait_order), unname(known_order))
  final_order <- c(unname(known_order), other_order)

  d <- d %>%
    mutate(
      trait_label = factor(trait_label, levels = final_order),
      method_label = factor(method_label, levels = unique(apply_method_labels(METHODS)))
    )

  y_lims <- if (ABS_VALUES) c(0, 1) else c(-1, 1)
  y_lab  <- cols$corr_label

  p <- ggplot(d, aes(x = trait_label, y = plot_val, fill = method_label)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, color = NA) +
    geom_text(aes(label = stars, group = method_label),
              position = position_dodge(width = 0.8),
              vjust = -0.2, size = 4, fontface = "bold", color = "black") +
    scale_fill_brewer(palette = "Set2", name = NULL) +
    scale_y_continuous(limits = y_lims, expand = expansion(mult = c(0, 0.08))) +
    labs(
      x = NULL,
      y = y_lab,
      title = sprintf("OpenTargets correlation by trait (%s, %s)", str_to_title(corr_type), subset_name),
      subtitle = sprintf("Tool: %s | Metric: %s | Values: %s",
                         TOOL_BASE, METRIC, ifelse(ABS_VALUES, "absolute", "raw"))
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(face = "bold")
    )

  out_file <- file.path(
    OUT_DIR,
    sprintf("ot_corr_grouped_bars_%s_%s_%s_%s.pdf",
            TOOL_BASE, METRIC, corr_type, slugify(subset_name))
  )
  ggsave(out_file, p, width = 12, height = 4.5, limitsize = FALSE)
  msg("Saved: %s", out_file)
}

# Drive plots
for (subset_name in SUBSETS) {
  for (ct in CORR_TYPES) {
    plot_one(dat, subset_name, ct)
  }
}

msg("Done.")