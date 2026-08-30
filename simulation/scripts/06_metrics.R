#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 06 - Per-replicate metrics (instructions Sec. 6)
#
# Joins the ground-truth table to the three methods' pathway statistics and
# computes, for each method:
#   1. contaminated-null FPR      (the money metric)
#   2. clean-null FPR             (calibration check)
#   3. power by effect tier       (does GSR over-correct?)
#   4. ranking accuracy vs planted truth: AUROC, AUPRC, Spearman, top-K precision
#   5. mechanism: slope of -log10(p) on the number of signal-carrying hubs
#
# Methods:
#   Original = raw MAGMA competitive p from the real simulated database
#   PS       = keeppathsize empirical p / standardised effect size
#   GSR      = birewire     empirical p / standardised effect size
#
# Ranking statistics are reported twice: on -log10(p) (what significance calls
# use, but tied at 1/(n_null+1) granularity for the empirical methods) and on
# the standardised effect size (continuous, so it is the primary ranking stat
# for the empirical methods). SES is undefined for Original, which uses p only.
# ---------------------------------------------------------------------------

parse_args_defaults <- list(
  ground_truth   = "",
  real_gsa       = "",     # MAGMA .gsa.out on the real simulated database
  emp_birewire   = "",     # calc_empirical.r output, birewire
  emp_keeppath   = "",     # calc_empirical.r output, keeppathsize
  condition      = "cond",
  replicate      = 1,
  alpha          = 0.05,
  top_k          = 30,
  out_prefix     = "metrics"
)
source(file.path(Sys.getenv("SIM_SCRIPTS", unset = "."), "sim_common.R"))
args <- parse_args(parse_args_defaults)
log_header("06 metrics | ", args$condition, " rep ", args$replicate)

# ---- rank-based helpers (no extra package dependencies) ------------------
auroc <- function(score, pos) {
  # Mann-Whitney U with mid-ranks, so ties count as half-credit
  r <- rank(score, ties.method = "average")
  n1 <- sum(pos); n0 <- sum(!pos)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

auprc <- function(score, pos) {
  # Average precision. Tied scores are resolved by giving every member of a tie
  # group the precision at the END of the group, which is the standard
  # tie-safe (non-optimistic) convention.
  n1 <- sum(pos)
  if (n1 == 0) return(NA_real_)
  o <- order(score, decreasing = TRUE)
  s <- score[o]; p <- pos[o]
  grp <- cumsum(c(TRUE, diff(s) != 0))
  tp_end <- cumsum(p)[cumsum(rle(grp)$lengths)][grp]
  n_end  <- (seq_along(p))[cumsum(rle(grp)$lengths)][grp]
  sum((tp_end / n_end)[p]) / n1
}

topk_precision <- function(score, pos, k) {
  o <- order(score, decreasing = TRUE)
  mean(pos[o][seq_len(min(k, length(pos)))])
}

# ---- load ----------------------------------------------------------------
gt <- fread(args$ground_truth)

read_real <- function(f) {
  raw <- readLines(f, warn = FALSE)
  hdr <- which(!startsWith(raw, "#") & nzchar(trimws(raw)))[1]
  d <- fread(text = paste(raw[hdr:length(raw)], collapse = "\n"), header = TRUE)
  if (!"FULL_NAME" %in% names(d)) d[, FULL_NAME := VARIABLE]
  d[, .(pathway_id = FULL_NAME, ngenes = NGENES, p = P, beta = BETA)]
}

stats <- list()
real <- read_real(args$real_gsa)
stats[["Original"]] <- real[, .(pathway_id, p, ses = NA_real_, ngenes)]

for (m in c(GSR = "emp_birewire", PS = "emp_keeppath")) {
  f <- args[[m]]
  if (!nzchar(f) || !file.exists(f)) next
  e <- fread(f)
  nm <- names(which(c(GSR = "emp_birewire", PS = "emp_keeppath") == m))
  stats[[nm]] <- e[, .(pathway_id = pathway_name, p = empirical_pval,
                       ses = std_effect_size, ngenes)]
}

pathway_stats <- rbindlist(lapply(names(stats), function(nm)
  merge(gt, stats[[nm]], by = "pathway_id")[, method := nm]), fill = TRUE)
pathway_stats[, `:=`(condition = args$condition, replicate = as.integer(args$replicate))]
pathway_stats[, significant := p < args$alpha]
pathway_stats[, neglog10p := -log10(pmax(p, .Machine$double.xmin))]

missing <- pathway_stats[, .(n = .N), by = method][n != nrow(gt)]
if (nrow(missing) > 0) {
  warning("Method(s) with != ", nrow(gt), " pathways: ",
          paste(missing$method, missing$n, collapse = "; "))
}

# ---- metrics -------------------------------------------------------------
metrics <- pathway_stats[, {
  cn <- class == "contaminated_null"; cl <- class == "clean_null"
  pos <- is_truly_enriched
  # Both ranking statistics are ALWAYS emitted, as NA where undefined (Original
  # has no SES), so every method contributes the same columns to the by-group.
  rank_stats <- list(neglog10p = neglog10p, ses = ses)

  base <- list(
    fpr_contaminated_null = if (any(cn)) mean(significant[cn]) else NA_real_,
    fpr_clean_null        = if (any(cl)) mean(significant[cl]) else NA_real_,
    fpr_all_null          = mean(significant[!pos]),
    power_low             = mean(significant[class == "enriched_low"]),
    power_med             = mean(significant[class == "enriched_med"]),
    power_high            = mean(significant[class == "enriched_high"]),
    power_overall         = mean(significant[pos]),
    n_contaminated_null   = sum(cn),
    n_clean_null          = sum(cl)
  )

  rank_metrics <- unlist(lapply(names(rank_stats), function(sn) {
    s <- rank_stats[[sn]]
    ok <- is.finite(s)
    vals <- if (sum(ok) < 3 || length(unique(s[ok])) < 2) {
      list(NA_real_, NA_real_, NA_real_, NA_real_)
    } else {
      list(auroc(s[ok], pos[ok]),
           auprc(s[ok], pos[ok]),
           suppressWarnings(cor(s[ok], truth_rank[ok], method = "spearman")),
           topk_precision(s[ok], pos[ok], as.integer(args$top_k)))
    }
    setNames(vals, paste0(c("auroc_", "auprc_", "spearman_truth_", "topk_precision_"), sn))
  }), recursive = FALSE)

  # Mechanism (Sec. 6.5): within NULL pathways only, so the only source of
  # signal is the hub genes. Slope > 0 means hub contamination drives the score.
  nl <- !pos
  slope <- slope_all <- NA_real_
  if (sum(nl) > 2 && length(unique(n_signal_hubs[nl])) > 1) {
    slope <- unname(coef(lm(neglog10p[nl] ~ n_signal_hubs[nl]))[2])
  }
  if (length(unique(n_signal_hubs)) > 1) {
    slope_all <- unname(coef(lm(neglog10p ~ n_signal_hubs))[2])
  }
  c(base, rank_metrics,
    list(hub_slope_nulls = slope, hub_slope_all = slope_all))
}, by = .(condition, replicate, method)]

long <- melt(metrics, id.vars = c("condition", "replicate", "method"),
             variable.name = "metric", value.name = "value")

fwrite(pathway_stats, paste0(args$out_prefix, "_pathway_stats.tsv"), sep = "\t")
fwrite(long, paste0(args$out_prefix, "_metrics.tsv"), sep = "\t")

cat("\n")
print(dcast(long[metric %in% c("fpr_contaminated_null", "fpr_clean_null",
                               "power_low", "power_med", "power_high",
                               "auprc_ses", "hub_slope_nulls")],
            metric ~ method, value.var = "value"))
cat("\n06 metrics: done\n")
