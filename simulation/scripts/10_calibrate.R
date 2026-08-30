#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 10 - Calibrate the per-tier mean Z-shifts mu (instructions Sec. 3).
#
# Sec. 3 forbids hard-coded magic numbers: mu must be calibrated ONCE, in the
# 0%-overlap condition, to land on interpretable operating points and then be
# reused across every condition.
#
# The pilot sweep runs the 0%-overlap architecture with ALL enriched pathways
# sharing one mu, over a grid of mu values. This script fits Original's power
# curve across that grid and inverts it:
#     mu_high = mu at power 0.80
#     mu_med  = mu at power 0.50
#     mu_low  = mu at power 0.20
#
# Note the clean-null FPR target of ~0.05 is NOT tunable through mu -- it is a
# property of the test, not the signal. It is therefore reported as a pass/fail
# calibration CHECK here rather than being solved for.
# ---------------------------------------------------------------------------
parse_args_defaults <- list(
  pilot_metrics = "",        # concatenated 06 outputs from the pilot sweep
  pilot_grid    = "",        # tsv: condition, mu  (the grid point each ran at)
  method        = "Original",
  target_low    = 0.20,
  target_med    = 0.50,
  target_high   = 0.80,
  out_prefix    = "calibration"
)
source(file.path(Sys.getenv("SIM_SCRIPTS", unset = "."), "sim_common.R"))
args <- parse_args(parse_args_defaults)
log_header("10 calibrate")

met <- if (dir.exists(args$pilot_metrics)) {
  rbindlist(lapply(list.files(args$pilot_metrics, pattern = "_metrics\\.tsv$",
                              full.names = TRUE, recursive = TRUE), fread), fill = TRUE)
} else fread(args$pilot_metrics)
grid <- fread(args$pilot_grid)          # condition -> mu
met <- merge(met, grid, by = "condition")

pw <- met[method == args$method & metric == "power_overall" & !is.na(value)]
if (nrow(pw) == 0) stop("No power_overall rows for method ", args$method)
curve <- pw[, .(power = mean(value), sd = sd(value), n = .N), by = mu][order(mu)]
cat("\nPilot power curve (", args$method, "):\n", sep = ""); print(curve)

if (nrow(curve) < 3) stop("Need at least 3 grid points to invert the power curve")
if (any(diff(curve$power) < -0.05)) {
  warning("Power is not monotone in mu across the pilot grid - check the pilot run")
}

# Logistic fit on the log-mu scale: power is a smooth saturating function of
# effect size, and log-mu keeps the fit well behaved over a multiplicative grid.
fit <- glm(power ~ log(mu), family = quasibinomial(), weights = rep(1, nrow(curve)),
           data = curve)
inv <- function(target) {
  b <- coef(fit)
  eta <- log(target / (1 - target))
  mu <- exp((eta - b[[1]]) / b[[2]])
  # keep the solution inside the pilot grid, extrapolation here is not trustworthy
  list(mu = mu, in_range = mu >= min(curve$mu) && mu <= max(curve$mu))
}

sol <- lapply(c(low = args$target_low, med = args$target_med, high = args$target_high), inv)
mus <- vapply(sol, `[[`, numeric(1), "mu")
in_range <- vapply(sol, `[[`, logical(1), "in_range")

fpr <- met[metric == "fpr_clean_null" & !is.na(value), .(fpr = mean(value)), by = .(method, mu)]

res <- data.table(tier = names(mus), target_power = c(args$target_low, args$target_med, args$target_high),
                  mu = as.numeric(mus), within_pilot_grid = in_range)
fwrite(res, paste0(args$out_prefix, "_mu.tsv"), sep = "\t")
fwrite(curve, paste0(args$out_prefix, "_power_curve.csv"))
fwrite(fpr,   paste0(args$out_prefix, "_clean_null_fpr.csv"))

# Nextflow reads these back as --mu_low/--mu_med/--mu_high
writeLines(sprintf("mu_low=%.6f\nmu_med=%.6f\nmu_high=%.6f",
                   mus[["low"]], mus[["med"]], mus[["high"]]),
           paste0(args$out_prefix, "_mu.env"))

lines <- c(
  sprintf("calibration method            : %s (raw MAGMA p < 0.05), 0%% overlap", args$method),
  sprintf("pilot grid mu                 : %s", paste(sprintf("%.3f", curve$mu), collapse = ", ")),
  sprintf("pilot power                   : %s", paste(sprintf("%.2f", curve$power), collapse = ", ")),
  sprintf("mu_low  (power %.2f)          : %.4f%s", args$target_low, mus[["low"]],
          if (in_range[["low"]]) "" else "   [EXTRAPOLATED - widen the pilot grid]"),
  sprintf("mu_med  (power %.2f)          : %.4f%s", args$target_med, mus[["med"]],
          if (in_range[["med"]]) "" else "   [EXTRAPOLATED - widen the pilot grid]"),
  sprintf("mu_high (power %.2f)          : %.4f%s", args$target_high, mus[["high"]],
          if (in_range[["high"]]) "" else "   [EXTRAPOLATED - widen the pilot grid]"),
  "",
  "clean-null FPR across the pilot grid (target ~0.05, NOT tunable via mu):",
  paste0("  ", fpr$method, " mu=", sprintf("%.3f", fpr$mu), " -> ", sprintf("%.4f", fpr$fpr))
)
writeLines(lines, paste0(args$out_prefix, "_report.txt"))
cat("\n", paste(lines, collapse = "\n"), "\n\n10 calibrate: done\n")
