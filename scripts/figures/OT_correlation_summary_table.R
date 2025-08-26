library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

options(dplyr.summarise.inform = FALSE)

# ------------------------------------------------------------------------------
# Config
# ------------------------------------------------------------------------------
INPUT_ROOT <- "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/opentargets_correlation"
TOOL_BASE  <- "magma"

# Include whatever subsets you want to evaluate; add "Top 50 Pathways" if present
SUBSETS <- c("Top 50 Pathways","Top 100 Pathways","Top 250 Pathways","Top 500 Pathways","All Pathways")

# Methods to compare (must match your CSV 'method' values)
METHODS3 <- c("BireWire_EmpPvalStdBeta","KeepPathSize_EmpPvalStdBeta","PvalueBeta")

# Abbreviations for tables/prints
ABBR <- c(
  "BireWire_EmpPvalStdBeta"     = "EmpStd (GF+PS)",
  "KeepPathSize_EmpPvalStdBeta" = "EmpStd (PS)",
  "PvalueBeta"                  = "RawP + ES"
)

# Significance threshold; “strict” summaries keep (trait, subset) with at least one method p < alpha
alpha <- 0.000925925925

# ------------------------------------------------------------------------------
# Read all correlation summaries
# ------------------------------------------------------------------------------
files <- list.files(INPUT_ROOT, pattern = "_rank_correlation_summary\\.csv$", recursive = TRUE, full.names = TRUE)
stopifnot(length(files) > 0)

rows <- vector("list", length(files))
k <- 0
for (f in files) {
  bn <- basename(f)
  m  <- regexec("^(.+?)_([A-Za-z0-9]+)_rank_correlation_summary\\.csv$", bn)
  mm <- regmatches(bn, m)[[1]]
  if (length(mm) != 3) next
  d <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(d)) next
  d[, trait := mm[2]]
  d[, tool_base := mm[3]]
  k <- k + 1
  rows[[k]] <- d
}
if (k == 0) stop("No parsable CSVs found")
dat <- rbindlist(rows[seq_len(k)], fill = TRUE) |> as_tibble()

# Keep only desired tool, subsets, and methods
dat <- dat %>%
  filter(tool_base == TOOL_BASE, subset %in% SUBSETS, method %in% METHODS3)

# ------------------------------------------------------------------------------
# Helper: summarise “BireWire better (more negative) than both competitors”
# Returns a list: $summary (per-subset table with publication-friendly labels)
# and $rows (row-level deltas).
# If alpha is finite, keep only (trait, subset) where at least one method has p < alpha.
# If alpha is Inf or NULL, no p-value filtering (use all rows).
# Deltas are defined so that POSITIVE = BireWire more negative (better).
# ------------------------------------------------------------------------------
summarise_three <- function(df, corr_col, pval_col, alpha = Inf, methods3 = METHODS3) {
  need_cols <- c(corr_col, pval_col, "trait","subset","method")
  if (!all(need_cols %in% names(df))) {
    warning(sprintf("Missing columns for corr=%s pval=%s; skipping.", corr_col, pval_col))
    return(NULL)
  }

  # Keep only the three methods
  df3 <- dplyr::filter(df, method %in% methods3)

  # Optional restriction by significance at the trait x subset level
  if (!is.null(alpha) && is.finite(alpha)) {
    df3 <- df3 |>
      dplyr::group_by(trait, subset) |>
      dplyr::filter(any(.data[[pval_col]] < alpha, na.rm = TRUE)) |>
      dplyr::ungroup()
  }
  if (nrow(df3) == 0) {
    message(sprintf("No rows to consider for %s (after filtering).", corr_col))
    return(NULL)
  }

  # Build wide safely
  df3_small <- dplyr::select(df3, trait, subset, method, value = .data[[corr_col]])
  wide <- tidyr::pivot_wider(df3_small, names_from = method, values_from = value)

  # Ensure all three methods are present
  if (!all(methods3 %in% names(wide))) {
    missing <- setdiff(methods3, names(wide))
    message("Missing method columns: ", paste(missing, collapse = ", "), " — cannot compare.")
    return(NULL)
  }

  # Keep rows where all three methods have non-NA values
  wide <- dplyr::filter(
    wide,
    !is.na(.data[[methods3[1]]]) &
    !is.na(.data[[methods3[2]]]) &
    !is.na(.data[[methods3[3]]])
  )
  if (nrow(wide) == 0) {
    message("No comparable rows (all three methods present) for ", corr_col)
    return(NULL)
  }

  # Signed correlations (more negative is better)
  bw   <- wide[[methods3[1]]]
  kps  <- wide[[methods3[2]]]
  pvb  <- wide[[methods3[3]]]

  # Deltas: competitor − BireWire; positive means BireWire is more negative (better)
  delta_bw_kps  <- kps - bw
  delta_bw_pvb  <- pvb - bw
  best_other    <- pmin(kps, pvb)     # most negative competitor
  delta_bw_both <- best_other - bw    # >0 only if BW is more negative than both

  rows <- tibble::tibble(
    trait  = wide$trait,
    subset = wide$subset,
    delta_bw_kps  = delta_bw_kps,
    delta_bw_pvb  = delta_bw_pvb,
    delta_bw_both = delta_bw_both,
    win_over_both = delta_bw_both > 0
  )

  # Base summary
  summary <- rows |>
    dplyr::group_by(subset) |>
    dplyr::summarise(
      n_rows_considered               = dplyr::n(),
      wins_bw_over_keepsize           = sum(delta_bw_kps  > 0),
      wins_bw_over_pvaluebeta         = sum(delta_bw_pvb  > 0),
      wins_bw_over_both_competitors   = sum(delta_bw_both > 0),
      ties_strict                     = sum(delta_bw_both == 0),
      median_delta_bw_vs_keepsize     = median(delta_bw_kps),
      median_delta_bw_vs_pvaluebeta   = median(delta_bw_pvb),
      median_delta_bw_vs_both         = median(delta_bw_both),
      .groups = "drop"
    )

  # Publication-friendly column labels using abbreviations
  corr_sym <- if (grepl("kendall", corr_col, ignore.case = TRUE)) "τ" else "ρ"
  bw_short  <- ABBR[methods3[1]]
  kps_short <- ABBR[methods3[2]]
  pvb_short <- ABBR[methods3[3]]

  summary <- summary |>
    dplyr::rename(
      `N (trait–subset pairs)` = n_rows_considered
    ) |>
    dplyr::rename(
      !!paste0(bw_short, " > ", kps_short, " (count)") := wins_bw_over_keepsize,
      !!paste0(bw_short, " > ", pvb_short, " (count)") := wins_bw_over_pvaluebeta,
      !!paste0(bw_short, " > both (count)")           := wins_bw_over_both_competitors,
      `Ties (= best competitor)`                      = ties_strict
    ) |>
    dplyr::rename(
      !!paste0("Median Δ", corr_sym, " (", kps_short, " − ", bw_short, ")") := median_delta_bw_vs_keepsize,
      !!paste0("Median Δ", corr_sym, " (", pvb_short, " − ", bw_short, ")") := median_delta_bw_vs_pvaluebeta,
      !!paste0("Median Δ", corr_sym, " (Best − ", bw_short, ")")            := median_delta_bw_vs_both
    )

  list(summary = summary, rows = rows)
}

# Small helper to print trait-level scoreboards
print_trait_scoreboard <- function(rows_tbl, label) {
  if (is.null(rows_tbl)) {
    cat("\n=== Trait scoreboard (", label, ") ===\nNo rows.\n", sep = "")
    return(invisible(NULL))
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
# Produce summaries (MEAN SCORE)
# ------------------------------------------------------------------------------
# Spearman (mean score): significant-only
spearman_strict <- summarise_three(
  df = dat,
  corr_col = "rank_mean_score_correlation",
  pval_col = "rank_mean_score_pvalue",
  alpha = alpha,
  methods3 = METHODS3
)
spearman_summary_strict <- if (is.null(spearman_strict)) NULL else spearman_strict$summary
spearman_rows_strict    <- if (is.null(spearman_strict)) NULL else spearman_strict$rows

# Kendall (tau, mean score): significant-only
kendall_strict <- summarise_three(
  df = dat,
  corr_col = "kendall_mean_score_tau",
  pval_col = "kendall_mean_score_pvalue",
  alpha = alpha,
  methods3 = METHODS3
)
kendall_summary_strict <- if (is.null(kendall_strict)) NULL else kendall_strict$summary
kendall_rows_strict    <- if (is.null(kendall_strict)) NULL else kendall_strict$rows

# Spearman (mean score): ALL rows (no p-value filter)
spearman_all <- summarise_three(
  df = dat,
  corr_col = "rank_mean_score_correlation",
  pval_col = "rank_mean_score_pvalue",
  alpha = Inf,
  methods3 = METHODS3
)
spearman_summary_all <- if (is.null(spearman_all)) NULL else spearman_all$summary
spearman_rows_all    <- if (is.null(spearman_all)) NULL else spearman_all$rows

# Kendall (tau, mean score): ALL rows (no p-value filter)
kendall_all <- summarise_three(
  df = dat,
  corr_col = "kendall_mean_score_tau",
  pval_col = "kendall_mean_score_pvalue",
  alpha = Inf,
  methods3 = METHODS3
)
kendall_summary_all <- if (is.null(kendall_all)) NULL else kendall_all$summary
kendall_rows_all    <- if (is.null(kendall_all)) NULL else kendall_all$rows

# ------------------------------------------------------------------------------
# Produce summaries (EVIDENCE DENSITY)
# ------------------------------------------------------------------------------
# Spearman (evidence density): significant-only
spearman_density_strict <- summarise_three(
  df = dat,
  corr_col = "rank_evidence_density_correlation",
  pval_col = "rank_evidence_density_pvalue",
  alpha = alpha,
  methods3 = METHODS3
)
spearman_density_summary_strict <- if (is.null(spearman_density_strict)) NULL else spearman_density_strict$summary
spearman_density_rows_strict    <- if (is.null(spearman_density_strict)) NULL else spearman_density_strict$rows

# Kendall (evidence density): significant-only
kendall_density_strict <- summarise_three(
  df = dat,
  corr_col = "kendall_evidence_density_tau",
  pval_col = "kendall_evidence_density_pvalue",
  alpha = alpha,
  methods3 = METHODS3
)
kendall_density_summary_strict <- if (is.null(kendall_density_strict)) NULL else kendall_density_strict$summary
kendall_density_rows_strict    <- if (is.null(kendall_density_strict)) NULL else kendall_density_strict$rows

# Spearman (evidence density): ALL rows
spearman_density_all <- summarise_three(
  df = dat,
  corr_col = "rank_evidence_density_correlation",
  pval_col = "rank_evidence_density_pvalue",
  alpha = Inf,
  methods3 = METHODS3
)
spearman_density_summary_all <- if (is.null(spearman_density_all)) NULL else spearman_density_all$summary
spearman_density_rows_all    <- if (is.null(spearman_density_all)) NULL else spearman_density_all$rows

# Kendall (evidence density): ALL rows
kendall_density_all <- summarise_three(
  df = dat,
  corr_col = "kendall_evidence_density_tau",
  pval_col = "kendall_evidence_density_pvalue",
  alpha = Inf,
  methods3 = METHODS3
)
kendall_density_summary_all <- if (is.null(kendall_density_all)) NULL else kendall_density_all$summary
kendall_density_rows_all    <- if (is.null(kendall_density_all)) NULL else kendall_density_all$rows

# ------------------------------------------------------------------------------
# Legend for outputs
# ------------------------------------------------------------------------------
cat(
  "\nLegend: ",
  ABBR[METHODS3[1]], " = BireWire (empirical p + standardized effect; accounts for gene frequency + pathway size); ",
  ABBR[METHODS3[2]], " = KeepPathSize (empirical p + standardized effect; accounts for pathway size); ",
  ABBR[METHODS3[3]], " = RawP + ES. ",
  "Δ = comparator − ", ABBR[METHODS3[1]], "; positive Δ means ", ABBR[METHODS3[1]], " is more negative (better).\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# Print results
# ------------------------------------------------------------------------------
cat("\n=== Spearman (mean score, restricted to significant) ===\n")
if (is.null(spearman_summary_strict)) message("No Spearman strict summary.") else print(spearman_summary_strict)

cat("\n=== Kendall (mean score, restricted to significant) ===\n")
if (is.null(kendall_summary_strict)) message("No Kendall strict summary.") else print(kendall_summary_strict)

cat("\n=== Spearman (mean score, ALL correlations) ===\n")
if (is.null(spearman_summary_all)) message("No Spearman all-rows summary.") else print(spearman_summary_all)

cat("\n=== Kendall (mean score, ALL correlations) ===\n")
if (is.null(kendall_summary_all)) message("No Kendall all-rows summary.") else print(kendall_summary_all)

cat("\n=== Spearman (evidence density, restricted to significant) ===\n")
if (is.null(spearman_density_summary_strict)) message("No Spearman density strict summary.") else print(spearman_density_summary_strict)

cat("\n=== Kendall (evidence density, restricted to significant) ===\n")
if (is.null(kendall_density_summary_strict)) message("No Kendall density strict summary.") else print(kendall_density_summary_strict)

cat("\n=== Spearman (evidence density, ALL correlations) ===\n")
if (is.null(spearman_density_summary_all)) message("No Spearman density all-rows summary.") else print(spearman_density_summary_all)

cat("\n=== Kendall (evidence density, ALL correlations) ===\n")
if (is.null(kendall_density_summary_all)) message("No Kendall density all-rows summary.") else print(kendall_density_summary_all)

# ------------------------------------------------------------------------------
# Trait-level scoreboards
# ------------------------------------------------------------------------------
print_trait_scoreboard(spearman_rows_strict,         "Spearman mean score (significant-only)")
print_trait_scoreboard(spearman_rows_all,            "Spearman mean score (all rows)")
print_trait_scoreboard(kendall_rows_strict,          "Kendall mean score (significant-only)")
print_trait_scoreboard(kendall_rows_all,             "Kendall mean score (all rows)")

print_trait_scoreboard(spearman_density_rows_strict, "Spearman evidence density (significant-only)")
print_trait_scoreboard(spearman_density_rows_all,    "Spearman evidence density (all rows)")
print_trait_scoreboard(kendall_density_rows_strict,  "Kendall evidence density (significant-only)")
print_trait_scoreboard(kendall_density_rows_all,     "Kendall evidence density (all rows)")



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

# Spearman mean score: summaries
write_if_any(
  spearman_summary_all,
  file.path(OUTPUT_DIR, "spearman_mean_score_summary_all.csv")
)
write_if_any(
  spearman_summary_strict,
  file.path(OUTPUT_DIR, "spearman_mean_score_summary_significant_only.csv")
)

# Optional: also write row-level deltas
# write_if_any(spearman_rows_all,
#              file.path(OUTPUT_DIR, "spearman_mean_score_rows_all.csv"))
# write_if_any(spearman_rows_strict,
#              file.path(OUTPUT_DIR, "spearman_mean_score_rows_significant_only.csv"))