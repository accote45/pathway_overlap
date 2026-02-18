library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

options(dplyr.summarise.inform = FALSE)

# ------------------------------------------------------------------------------
# Config
# ------------------------------------------------------------------------------
INPUT_ROOT <- "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/opentargets_correlation"
TOOL_BASE  <- "magma"

# Include whatever subsets you want to evaluate
SUBSETS <- c("Top 100 Pathways","Top 250 Pathways","Top 500 Pathways")

# Methods to compare (must match your CSV 'method' values)
METHODS3 <- c("BireWire_EmpPvalStdBeta","KeepPathSize_EmpPvalStdBeta","PvalueBeta")

# Abbreviations for tables/prints
ABBR <- c(
  "BireWire_EmpPvalStdBeta"     = "EmpStd (GF+PS)",
  "KeepPathSize_EmpPvalStdBeta" = "EmpStd (PS)",
  "PvalueBeta"                  = "RawP + ES"
)

n_traits = 7

# Significance thresholds
alpha <- 0.05/(n_traits*length(METHODS3)*length(SUBSETS))  # Bonferroni
alpha_lenient <- 0.05    # p<0.05


# ------------------------------------------------------------------------------
# Read all correlation summaries
# ------------------------------------------------------------------------------
files <- list.files(INPUT_ROOT, pattern = "_rank_correlation_summary\\.csv$", recursive = TRUE, full.names = TRUE)
stopifnot(length(files) > 0)

parse_one <- function(f) {
  bn <- basename(f)
  mm <- regmatches(bn, regexec("^(.+?)_([A-Za-z0-9]+)_rank_correlation_summary\\.csv$", bn))[[1]]
  if (length(mm) != 3) return(NULL)
  d <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  d[, trait := mm[2]]
  d[, tool_base := mm[3]]
  d
}
dat <- rbindlist(Filter(Negate(is.null), lapply(files, parse_one)), fill = TRUE)

# Keep only desired tool, subsets, and methods
dat <- dat %>%
  filter(tool_base == TOOL_BASE, subset %in% SUBSETS, method %in% METHODS3)

# ------------------------------------------------------------------------------
# Core comparison: BireWire more negative than comparators
# ------------------------------------------------------------------------------
summarise_three <- function(df, corr_col, pval_col, alpha = Inf, methods3 = METHODS3) {
  need_cols <- c(corr_col, pval_col, "trait","subset","method")
  if (!all(need_cols %in% names(df))) {
    warning(sprintf("Missing columns for corr=%s pval=%s; skipping.", corr_col, pval_col))
    return(NULL)
  }

  df3 <- dplyr::filter(df, method %in% methods3)

  # Optional restriction by significance at the trait x subset level
  if (is.finite(alpha)) {
    df3 <- df3 %>%
      group_by(trait, subset) %>%
      filter(any(.data[[pval_col]] < alpha, na.rm = TRUE)) %>%
      ungroup()
  }
  if (nrow(df3) == 0) return(NULL)

  # Wide with the three methods
  wide <- df3 %>%
    select(trait, subset, method, value = .data[[corr_col]]) %>%
    pivot_wider(names_from = method, values_from = value) %>%
    filter(if_all(all_of(methods3), ~ !is.na(.)))

  if (nrow(wide) == 0) return(NULL)

  # Signed correlations (more negative is better)
  bw  <- wide[[methods3[1]]]
  kps <- wide[[methods3[2]]]
  pvb <- wide[[methods3[3]]]

  # Deltas: competitor − BireWire; positive means BireWire is more negative (better)
  delta_bw_kps  <- kps - bw
  delta_bw_pvb  <- pvb - bw
  best_other    <- pmin(kps, pvb)
  delta_bw_both <- best_other - bw

  rows <- tibble::tibble(
    trait  = wide$trait,
    subset = wide$subset,
    delta_bw_kps  = delta_bw_kps,
    delta_bw_pvb  = delta_bw_pvb,
    delta_bw_both = delta_bw_both,
    win_over_both = delta_bw_both > 0
  )

  # Per-subset summary
  summary <- rows %>%
    group_by(subset) %>%
    summarise(
      n_rows_considered             = n(),
      wins_bw_over_keepsize         = sum(delta_bw_kps  > 0),
      wins_bw_over_pvaluebeta       = sum(delta_bw_pvb  > 0),
      wins_bw_over_both_competitors = sum(delta_bw_both > 0),
      ties_strict                   = sum(delta_bw_both == 0),
      median_delta_bw_vs_keepsize   = median(delta_bw_kps),
      median_delta_bw_vs_pvaluebeta = median(delta_bw_pvb),
      median_delta_bw_vs_both       = median(delta_bw_both),
      .groups = "drop"
    )

  # Labels
  corr_sym <- if (grepl("kendall", corr_col, TRUE)) "τ" else "ρ"
  bw_short  <- ABBR[methods3[1]]
  kps_short <- ABBR[methods3[2]]
  pvb_short <- ABBR[methods3[3]]

  summary <- summary %>%
    rename(`N (trait–subset pairs)` = n_rows_considered) %>%
    rename(
      !!paste0(bw_short, " > ", kps_short, " (count)") := wins_bw_over_keepsize,
      !!paste0(bw_short, " > ", pvb_short, " (count)") := wins_bw_over_pvaluebeta,
      !!paste0(bw_short, " > both (count)")           := wins_bw_over_both_competitors,
      `Ties (= best competitor)`                      = ties_strict
    ) %>%
    rename(
      !!paste0("Median Δ", corr_sym, " (", kps_short, " − ", bw_short, ")") := median_delta_bw_vs_keepsize,
      !!paste0("Median Δ", corr_sym, " (", pvb_short, " − ", bw_short, ")") := median_delta_bw_vs_pvaluebeta,
      !!paste0("Median Δ", corr_sym, " (Best − ", bw_short, ")")            := median_delta_bw_vs_both
    )

  list(summary = summary, rows = rows)
}

# Small helper to print trait-level scoreboards
print_trait_scoreboard <- function(rows_tbl, label) {
  if (is.null(rows_tbl) || nrow(rows_tbl) == 0) {
    cat("\n=== Trait scoreboard (", label, ") ===\nNo rows.\n", sep = ""); return(invisible(NULL))
  }
  cat("\n=== Trait scoreboard (", label, ") ===\n", sep = "")
  rows_tbl %>%
    group_by(trait) %>%
    summarise(
      wins_over_both_across_subsets = sum(win_over_both),
      total_subsets_evaluated       = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(wins_over_both_across_subsets)) %>%
    print()
}

# ------------------------------------------------------------------------------
# Run evaluations
# ------------------------------------------------------------------------------
# Spearman mean score
spearman_strict <- summarise_three(dat, "rank_mean_score_correlation", "rank_mean_score_pvalue", alpha = alpha, methods3 = METHODS3)
spearman_p05    <- summarise_three(dat, "rank_mean_score_correlation", "rank_mean_score_pvalue", alpha = alpha_lenient, methods3 = METHODS3)

spearman_summary_strict <- if (is.null(spearman_strict)) NULL else spearman_strict$summary
spearman_rows_strict    <- if (is.null(spearman_strict)) NULL else spearman_strict$rows
spearman_summary_p05    <- if (is.null(spearman_p05))    NULL else spearman_p05$summary
spearman_rows_p05       <- if (is.null(spearman_p05))    NULL else spearman_p05$rows

# ------------------------------------------------------------------------------
# Legend for outputs
# ------------------------------------------------------------------------------
cat(
  "\nLegend: ",
  ABBR[METHODS3[1]], " = BireWire (empirical p + standardized effect; gene frequency + pathway size). ",
  ABBR[METHODS3[2]], " = KeepPathSize (empirical p + standardized effect; pathway size). ",
  ABBR[METHODS3[3]], " = RawP + ES. ",
  "Δ = comparator − ", ABBR[METHODS3[1]], " (positive means ", ABBR[METHODS3[1]], " more negative/better).\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# Print results
# ------------------------------------------------------------------------------
cat("\n=== Spearman (mean score, p<0.05 filter) ===\n")
if (is.null(spearman_summary_p05)) message("No Spearman p<0.05 summary.") else print(spearman_summary_p05)

cat("\n=== Spearman (mean score, Bonferroni) ===\n")
if (is.null(spearman_strict)) message("No Spearman all-rows summary.") else print(spearman_strict)

# ------------------------------------------------------------------------------
# Trait-level scoreboards (Spearman mean score only)
# ------------------------------------------------------------------------------
print_trait_scoreboard(spearman_rows_strict, paste0("Spearman mean score (Bonferroni alpha = ", signif(alpha, 3), ")"))
print_trait_scoreboard(spearman_rows_p05, "Spearman mean score (p<0.05)")

# ------------------------------------------------------------------------------
# Write CSV outputs (Spearman mean score summaries)
# ------------------------------------------------------------------------------
OUTPUT_DIR <- file.path(INPUT_ROOT, "summaries")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

write_if_any <- function(tbl, path) {
  if (!is.null(tbl) && is.data.frame(tbl) && nrow(tbl) > 0) {
    data.table::fwrite(tbl, path)
    message("Wrote: ", path)
  } else {
    message("Skipped writing ", basename(path), " (no data).")
  }
}

write_if_any(spearman_summary_strict, file.path(OUTPUT_DIR, "spearman_mean_score_summary_significant_only.csv"))
write_if_any(spearman_summary_p05,    file.path(OUTPUT_DIR, "spearman_mean_score_summary_p05.csv"))

# ------------------------------------------------------------------------------
# Table 2 builder (reusable)
# ------------------------------------------------------------------------------
build_table2 <- function(rows_tbl) {
  if (is.null(rows_tbl) || nrow(rows_tbl) == 0) return(NULL)
  tab <- rows_tbl %>%
    group_by(subset) %>%
    summarise(
      N = n(),
      BW_better_than_RawP = sum(delta_bw_pvb > 0),
      RawP_better_than_BW = sum(delta_bw_pvb < 0),
      BW_better_than_KPS  = sum(delta_bw_kps > 0),
      BW_better_than_both = sum(win_over_both),
      prop_BW_gt_RawP     = BW_better_than_RawP / N,
      med_delta_RawP      = median(delta_bw_pvb),
      med_delta_KPS       = median(delta_bw_kps),
      .groups = "drop"
    ) %>%
    select(
      subset, N,
      `EmpStd (GF+PS) > RawP + ES (count)` = BW_better_than_RawP,
      `RawP + ES > EmpStd (GF+PS) (count)` = RawP_better_than_BW,
      `EmpStd (GF+PS) > EmpStd (PS) (count)` = BW_better_than_KPS,
      `EmpStd (GF+PS) > both (count)` = BW_better_than_both,
      `Prop BW > RawP` = prop_BW_gt_RawP,
      `Median Δρ (RawP − BW)` = med_delta_RawP,
      `Median Δρ (KPS − BW)` = med_delta_KPS
    )
  tab
}

# Build tables
tbl2_p05  <- build_table2(spearman_rows_p05)
tbl2_bonf <- build_table2(spearman_rows_strict)

# Write tables
write_if_any(tbl2_p05,  file.path(OUTPUT_DIR, "spearman_mean_score_table2_p05.csv"))
write_if_any(tbl2_bonf, file.path(OUTPUT_DIR, "spearman_mean_score_table2_bonf.csv"))

# ------------------------------------------------------------------------------
# Trait-level tables for manuscript (optional objects)
# ------------------------------------------------------------------------------
tbl_trait_p05 <- if (!is.null(spearman_rows_p05)) {
  spearman_rows_p05 %>%
    group_by(trait) %>%
    summarise(
      `EmpStd (GF+PS) > both (count)` = sum(win_over_both),
      `Subsets evaluated` = n(),
      `Proportion best` = `EmpStd (GF+PS) > both (count)` / `Subsets evaluated`,
      .groups = "drop"
    ) %>% arrange(desc(`EmpStd (GF+PS) > both (count)`))
} else NULL

tbl_trait_bonf <- if (!is.null(spearman_rows_strict)) {
  spearman_rows_strict %>%
    group_by(trait) %>%
    summarise(
      `EmpStd (GF+PS) > both (count)` = sum(win_over_both),
      `Subsets evaluated` = n(),
      `Proportion best` = `EmpStd (GF+PS) > both (count)` / `Subsets evaluated`,
      .groups = "drop"
    ) %>% arrange(desc(`EmpStd (GF+PS) > both (count)`))
} else NULL



# ------------------------------------------------------------------------------
# Figure D: horizontal bar chart with trait label next to the bar (p<0.05 only)
# Δρ = (RawP + ES) − (EmpStd (GF+PS)); right favors BireWire
# ------------------------------------------------------------------------------
subsets_to_show <- c("Top 100 Pathways","Top 250 Pathways","Top 500 Pathways")
subset_pal <- c(
  "Top 100 Pathways" = "#D73027",  # red
  "Top 250 Pathways" = "#FDC827",  # yellow
  "Top 500 Pathways" = "#4575B4"   # blue
)

if (!is.null(spearman_rows_p05) && nrow(spearman_rows_p05) > 0) {
  figD_dat <- spearman_rows_p05 %>%
    filter(subset %in% subsets_to_show) %>%
    mutate(
      subset = factor(subset, levels = subsets_to_show),
      # Re-label traits (now includes AD and T2D; safer CAD match)
      trait_label = dplyr::case_when(
        grepl("IBD|inflamm", trait, ignore.case = TRUE) ~ "IBD",
        grepl("\\bCAD\\b|coronary", trait, ignore.case = TRUE) ~ "CAD",
        grepl("\\bAD\\b|alzheimer", trait, ignore.case = TRUE) ~ "AD",
        grepl("\\bT2D\\b|type\\s*2.*diab", trait, ignore.case = TRUE) ~ "T2D",
        grepl("SCZ|schizo", trait, ignore.case = TRUE) ~ "SCZ",
        grepl("MDD|depress", trait, ignore.case = TRUE) ~ "MDD",
        grepl("breast", trait, ignore.case = TRUE) ~ "Breast cancer",
        TRUE ~ trait
      ),
      delta = delta_bw_pvb
    ) %>%
    arrange(subset, desc(abs(delta)))  # optional: largest effects at top

  if (nrow(figD_dat) > 0) {
    xmax <- max(abs(figD_dat$delta), na.rm = TRUE)
    xpad <- if (is.finite(xmax)) xmax * 0.12 else 0.1
    pdf('horiz_barchart_nomsig_results_comparison.pdf')
    p4 <- ggplot(figD_dat, aes(x = delta, y = trait_label, fill = subset)) +
      geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
      geom_col(width = 0.6, alpha = 0.95) +
      # Put the trait label next to the bar end (right for positive, left for negative)
      geom_text(
        aes(
          x = delta + ifelse(delta >= 0,  xpad * 0.05, -xpad * 0.05),
          label = trait_label
        ),
        hjust = ifelse(figD_dat$delta >= 0, 0, 1),
        size = 3.1,
        color = "black"
      ) +
      facet_grid(subset ~ ., scales = "free_y", space = "free_y") +
      scale_fill_manual(values = subset_pal, guide = "none") +
      coord_cartesian(xlim = c(-xmax - xpad, xmax + xpad)) +
      labs(
        x = "Δρ (RawP − BireWire); right favors BireWire",
        y = NULL,
        title = "Per-trait differences by subset (Spearman mean score, p<0.05)"
      ) +
      theme_classic(base_size = 12) +
      theme(
        strip.text.y = element_text(angle = 0),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_blank(),   # hide axis labels since we draw them next to bars
        axis.ticks.y = element_blank()
      )
    print(p4)
    dev.off()
  } else {
    message("No rows for Figure D after excluding 'All Pathways'.")
  }
} else {
  message("No data for Figure D (p<0.05).")
}

# ------------------------------------------------------------------------------
# Figures: grouped bar charts of Spearman ρ by method, grouped by trait
# ------------------------------------------------------------------------------

# Short trait labels
label_trait <- function(x) {
  dplyr::case_when(
    grepl("IBD|inflamm", x, ignore.case = TRUE) ~ "IBD",
    grepl("\\bCAD\\b|coronary", x, ignore.case = TRUE) ~ "CAD",
    grepl("\\bAD\\b|alzheimer", x, ignore.case = TRUE) ~ "AD",
    grepl("\\bT2D\\b|type\\s*2.*diab", x, ignore.case = TRUE) ~ "T2D",
    grepl("SCZ|schizo", x, ignore.case = TRUE) ~ "SCZ",
    grepl("MDD|depress", x, ignore.case = TRUE) ~ "MDD",
    grepl("breast", x, ignore.case = TRUE) ~ "Breast cancer",
    TRUE ~ x
  )
}

# Method display order and colors
method_levels <- unname(ABBR[METHODS3])
method_pal <- c(
  "EmpStd (GF+PS)" = "#1B9E77",
  "EmpStd (PS)"    = "#D95F02",
  "RawP + ES"      = "#524c9d"
)

make_methods_barplot <- function(subset_name, alpha_filter, title_suffix, filter_by_any_method = FALSE) {
  df <- dat %>%
    filter(subset == subset_name, method %in% METHODS3)

  if (filter_by_any_method) {
    df <- df %>%
      group_by(trait) %>%
      filter(any(rank_mean_score_pvalue < alpha_filter, na.rm = TRUE)) %>%
      ungroup()
  }

  df <- df %>%
    transmute(
      subset,
      trait_label  = label_trait(trait),
      method_label = ABBR[method],
      rho   = rank_mean_score_correlation,
      pval  = rank_mean_score_pvalue,
      # asterisks drawn on bars
      sig_label = dplyr::case_when(
        pval < alpha         ~ "***",     # Bonferroni
        pval < 0.001         ~ "**",
        pval < alpha_lenient ~ "*",
        TRUE                 ~ ""
      ),
      # legend categories (no “ns” so legend is clean)
      sig_cat_legend = dplyr::case_when(
        pval < alpha         ~ "*** Bonferroni",
        pval < 0.001         ~ "** p<0.001",
        pval < alpha_lenient ~ "* p<0.05",
        TRUE                 ~ NA_character_
      )
    ) %>%
    filter(!is.na(rho))

  if (nrow(df) == 0) return(NULL)

  # Order traits: preferred first (now includes AD and T2D)
  preferred <- c("IBD","CAD","Breast cancer","SCZ","T2D","AD","MDD")
  trait_levels <- c(preferred[preferred %in% df$trait_label],
                    setdiff(sort(unique(df$trait_label)), preferred))
  df$trait_label  <- factor(df$trait_label, levels = trait_levels)
  df$method_label <- factor(df$method_label, levels = method_levels)

  ypad <- 0.02
  dodge <- position_dodge(width = 0.72)

  p <- ggplot(df, aes(x = trait_label, y = rho, fill = method_label)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey70") +
    geom_col(position = dodge, width = 0.68, alpha = 0.95) +
    # draw asterisks on bars
    geom_text(
      aes(
        label = sig_label,
        y = ifelse(rho >= 0, rho + ypad, rho - ypad)
      ),
      position = dodge,
      size = 3.1,
      vjust = ifelse(df$rho >= 0, 0, 1)
    ) +
    # invisible points to produce a clean legend for significance
    geom_point(
      aes(shape = sig_cat_legend),
      position = dodge, size = 0, alpha = 0, show.legend = TRUE
    ) +
    scale_shape_manual(
      name   = "Significance",
      breaks = c("*** Bonferroni", "** p<0.001", "* p<0.05"),
      values = c(
        "*** Bonferroni" = 8,
        "** p<0.001"     = 4,
        "* p<0.05"       = 3
      ),
      labels = c(
        paste0("*** Bonferroni (p<", signif(alpha, 3), ")"),
        "** p<0.001",
        "* p<0.05"
      ),
      drop = TRUE,
      guide = guide_legend(
        override.aes = list(shape = NA, size = 0, alpha = 0)  # hide symbols in legend
      )
    ) +
    scale_fill_manual(values = method_pal, name = "Method") +
    labs(
      title = paste0("Spearman rho by method, grouped by trait (", subset_name, ")"),
      x = NULL,
      y = "Spearman rho (rank vs OT evidence)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid.minor = element_blank()
    )

  p
}

# Generate and save plots: one per subset (all traits; significance annotated with * ** ***)
for (ss in SUBSETS) {
  p_all <- make_methods_barplot(
    subset_name = ss,
    alpha_filter = alpha_lenient,
    title_suffix = "all traits; asterisks denote significance",
    filter_by_any_method = FALSE
  )
  if (!is.null(p_all)) {
    if (interactive()) print(p_all)  # avoid opening default PDF device in Rscript
    ggsave(
      filename = file.path(OUTPUT_DIR, paste0("spearman_methods_bars_", gsub("[^A-Za-z0-9]+","_", ss), ".pdf")),
      plot = p_all, width = 9, height = 4.2, device = cairo_pdf
    )
  } else {
    message("No data for ", ss, ".")
  }
}

# ------------------------------------------------------------------------------
# Dumbbell plots: two panels (facets) — BW vs RawP and BW vs KPS
#   Left endpoint = comparator; Right endpoint = EmpStd (GF+PS)
#   (OT analysis does not include BMI; no BMI label fixes needed)
# ------------------------------------------------------------------------------
make_methods_dumbbellplot <- function(subset_name, alpha_filter, title_suffix, filter_by_any_method = FALSE) {
  df <- dat %>%
    filter(subset == subset_name, method %in% METHODS3)

  if (filter_by_any_method) {
    df <- df %>%
      group_by(trait) %>%
      filter(any(rank_mean_score_pvalue < alpha_filter, na.rm = TRUE)) %>%
      ungroup()
  }

  # Prepare long -> wide for the three methods
  df0 <- df %>%
    transmute(
      trait_label  = label_trait(trait),
      method_id    = method,
      rho          = rank_mean_score_correlation,
      pval         = rank_mean_score_pvalue
    )
  if (nrow(df0) == 0) return(NULL)

  dfw <- df0 %>%
    tidyr::pivot_wider(
      names_from = method_id,
      values_from = c(rho, pval)
    )

  # Require all three methods present per trait
  need_cols <- c(
    "rho_BireWire_EmpPvalStdBeta",
    "rho_KeepPathSize_EmpPvalStdBeta",
    "rho_PvalueBeta",
    "pval_BireWire_EmpPvalStdBeta"
  )
  dfw <- dfw %>% dplyr::filter(if_all(all_of(need_cols), ~ !is.na(.)))
  if (nrow(dfw) == 0) return(NULL)

  # Trait ordering consistent with bar plots
  preferred <- c("IBD","CAD","Breast cancer","SCZ","T2D","AD","MDD")
  trait_levels <- c(preferred[preferred %in% dfw$trait_label],
                    setdiff(sort(unique(dfw$trait_label)), preferred))
  dfw$trait_label <- factor(dfw$trait_label, levels = trait_levels)

  # Build dumbbell data for two comparisons
  dumb_dat <- bind_rows(
    dfw %>%
      transmute(
        trait_label,
        comparison = "BW vs RawP",
        x    = .data[["rho_PvalueBeta"]],
        xend = .data[["rho_BireWire_EmpPvalStdBeta"]],
        pval_bw = .data[["pval_BireWire_EmpPvalStdBeta"]]
      ),
    dfw %>%
      transmute(
        trait_label,
        comparison = "BW vs KPS",
        x    = .data[["rho_KeepPathSize_EmpPvalStdBeta"]],
        xend = .data[["rho_BireWire_EmpPvalStdBeta"]],
        pval_bw = .data[["pval_BireWire_EmpPvalStdBeta"]]
      )
  ) %>%
    mutate(
      sig_label = dplyr::case_when(
        pval_bw < alpha         ~ "***",          # Bonferroni
        pval_bw < 0.001         ~ "**",
        pval_bw < alpha_lenient ~ "*",
        TRUE                    ~ ""
      )
    )

  # Labels from ABBR for figure text
  bw_lab  <- ABBR["BireWire_EmpPvalStdBeta"]     # "EmpStd (GF+PS)"
  kps_lab <- ABBR["KeepPathSize_EmpPvalStdBeta"] # "EmpStd (PS)"
  pvb_lab <- ABBR["PvalueBeta"]                  # "RawP + ES"

  # Keep data levels stable, but show ABBR in legend and facet strips
  comp_levels <- c("BW vs RawP","BW vs KPS")
  dumb_dat$comparison <- factor(dumb_dat$comparison, levels = comp_levels)

  # Colors: left endpoint = comparator; right endpoint = BW
  comp_pal <- setNames(
    as.character(c(method_pal[pvb_lab], method_pal[kps_lab])),
    comp_levels
  )
  bw_color <- unname(method_pal[bw_lab])

  # Facet strip labels use ABBR text
  strip_labs <- c(
    "BW vs RawP" = paste0(bw_lab, " vs ", pvb_lab),
    "BW vs KPS"  = paste0(bw_lab, " vs ", kps_lab)
  )

  # Two-panel plot (facets), one dumbbell per trait in each panel
  xpad <- 0.02
  p <- ggplot(dumb_dat, aes(y = trait_label)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey70") +
    geom_segment(aes(x = x, xend = xend, yend = trait_label), color = "grey60", linewidth = 0.9) +
    geom_point(aes(x = x, color = comparison), size = 2.6) +
    geom_point(aes(x = xend), color = bw_color, size = 2.6, show.legend = FALSE) +
    geom_text(
      aes(
        x = xend + ifelse(xend >= 0, xpad, -xpad),
        label = sig_label,
        hjust = ifelse(xend >= 0, 0, 1)
      ),
      size = 3.0
    ) +
    facet_wrap(~ comparison, nrow = 1, labeller = as_labeller(strip_labs)) +
    scale_color_manual(
      values = comp_pal,
      breaks = comp_levels,
      labels = c(pvb_lab, kps_lab),
      name   = "Comparator (left point)"
    ) +
    labs(
      title = paste0("Spearman \u03C1 dumbbells: ", bw_lab, " vs ", pvb_lab, " and ", kps_lab, " (", subset_name, ")"),
      subtitle = paste0("Left = comparator; Right = ", bw_lab, ". Asterisks show ", bw_lab, " significance."),
      x = "Spearman \u03C1 (rank vs OT evidence)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  p
}

# Generate and save dumbbell plots: one per subset
for (ss in SUBSETS) {
  p_db <- make_methods_dumbbellplot(
    subset_name = ss,
    alpha_filter = alpha_lenient,
    title_suffix = "all traits; asterisks denote BW significance",
    filter_by_any_method = FALSE
  )
  if (!is.null(p_db)) {
    if (interactive()) print(p_db)
    ggsave(
      filename = file.path(OUTPUT_DIR, paste0("spearman_methods_dumbbells_", gsub("[^A-Za-z0-9]+","_", ss), ".pdf")),
      plot = p_db, width = 10, height = 4.2, device = cairo_pdf
    )
  } else {
    message("No dumbbell data for ", ss, ".")
  }
}





