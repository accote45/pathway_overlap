#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 07 - Aggregate per-replicate metrics across the whole sweep (Sec. 6).
#
# Emits the tidy CSVs that back every figure:
#   summary_metrics.csv        - mean, sd, 95% CI and n per condition/method/metric
#   all_replicate_metrics.csv  - every replicate's value (for boxplots)
#   all_pathway_stats.csv      - per-pathway rows pooled over replicates
#   hub_regression.csv         - pooled -log10(p) vs n_signal_hubs fit per method
#   sanity_checks.txt
# ---------------------------------------------------------------------------
parse_args_defaults <- list(
  metrics_glob  = "*_metrics.tsv",
  pathway_glob  = "*_pathway_stats.tsv",
  in_dir        = ".",
  out_dir       = ".",
  conf_level    = 0.95
)
source(file.path(Sys.getenv("SIM_SCRIPTS", unset = "."), "sim_common.R"))
args <- parse_args(parse_args_defaults)
log_header("07 aggregate")
dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

read_all <- function(pat) {
  f <- list.files(args$in_dir, pattern = glob2rx(pat), full.names = TRUE, recursive = TRUE)
  if (length(f) == 0) stop("No files matching ", pat, " under ", args$in_dir)
  cat("Reading", length(f), "files matching", pat, "\n")
  rbindlist(lapply(f, fread), fill = TRUE)
}

met <- read_all(args$metrics_glob)
pw  <- read_all(args$pathway_glob)

# Order methods and conditions consistently everywhere downstream.
method_levels <- c("Original", "PS", "GSR")
met[, method := factor(method, levels = intersect(method_levels, unique(method)))]
pw[,  method := factor(method, levels = intersect(method_levels, unique(method)))]

n_reps <- met[, uniqueN(replicate), by = condition]
cat("Replicates per condition:\n"); print(n_reps)

z <- qnorm(1 - (1 - args$conf_level) / 2)
summ <- met[!is.na(value), .(
  mean = mean(value),
  sd   = sd(value),
  n    = .N,
  se   = sd(value) / sqrt(.N)
), by = .(condition, method, metric)]
summ[, `:=`(ci_lo = mean - z * se, ci_hi = mean + z * se)]

# Pooled mechanism regression (Sec. 6.5) across all replicates of a condition,
# restricted to null pathways so hub genes are the only signal source.
hubreg <- pw[is_truly_enriched == FALSE, {
  if (length(unique(n_signal_hubs)) < 2 || .N < 5) {
    .(slope = NA_real_, se = NA_real_, p = NA_real_, r2 = NA_real_, n = .N)
  } else {
    fit <- lm(neglog10p ~ n_signal_hubs)
    cf <- summary(fit)$coefficients
    .(slope = cf[2, 1], se = cf[2, 2], p = cf[2, 4],
      r2 = summary(fit)$r.squared, n = .N)
  }
}, by = .(condition, method)]

fwrite(summ,   file.path(args$out_dir, "summary_metrics.csv"))
fwrite(met,    file.path(args$out_dir, "all_replicate_metrics.csv"))
fwrite(pw,     file.path(args$out_dir, "all_pathway_stats.csv"))
fwrite(hubreg, file.path(args$out_dir, "hub_regression.csv"))

# ---- sanity checks (Sec. 8) ----------------------------------------------
chk <- character(0)
add <- function(ok, msg) chk <<- c(chk, sprintf("[%s] %s",
                                     if (is.na(ok)) "INFO" else if (isTRUE(ok)) "PASS" else "CHECK", msg))

clean <- summ[metric == "fpr_clean_null"]
for (i in seq_len(nrow(clean))) {
  add(abs(clean$mean[i] - 0.05) < 0.03,
      sprintf("clean-null FPR %s/%s = %.3f (target ~0.05)",
              clean$condition[i], clean$method[i], clean$mean[i]))
}
# the zero-overlap control condition(s): "ov00", "ov0", "ov00_something"
ov0 <- summ[grepl("^ov0+($|_)", condition)]
if (nrow(ov0) > 0) {
  for (mm in c("fpr_clean_null", "power_high", "auprc_ses")) {
    v <- ov0[metric == mm]
    if (nrow(v) > 1) {
      add(diff(range(v$mean, na.rm = TRUE)) < 0.10,
          sprintf("0%%-overlap control: %s spread across methods = %.3f (should be ~0)",
                  mm, diff(range(v$mean, na.rm = TRUE))))
    }
  }
}
for (i in seq_len(nrow(hubreg))) {
  if (is.na(hubreg$slope[i])) {
    add(NA, sprintf("hub slope %s/%s = not estimable (no variation in n_signal_hubs)",
                    hubreg$condition[i], hubreg$method[i]))
  } else {
    # GSR is expected to be flat; Original/PS are expected to be positive
    expect_flat <- hubreg$method[i] == "GSR"
    ok <- if (expect_flat) hubreg$p[i] > 0.05 || abs(hubreg$slope[i]) < 0.01 else hubreg$slope[i] > 0
    add(ok, sprintf("hub slope %s/%s = %+.4f (se %.4f, p %.3g)%s",
                    hubreg$condition[i], hubreg$method[i], hubreg$slope[i],
                    hubreg$se[i], hubreg$p[i],
                    if (expect_flat) " [expect ~0]" else " [expect > 0]"))
  }
}
writeLines(chk, file.path(args$out_dir, "sanity_checks.txt"))
cat("\n", paste(chk, collapse = "\n"), "\n\n07 aggregate: done\n")
