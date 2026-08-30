#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 08 - Figures 1-5 (instructions Sec. 7), plus the tidy CSV behind each one.
#
# Colour is categorical (method identity), assigned in a FIXED order that never
# changes when a condition is filtered out: Original = blue, PS = orange,
# GSR = aqua. These are slots 1-3 of the validated default categorical palette,
# the documented all-pairs-safe subset (worst-pair CVD dE 9.2, normal-vision
# 24.0 on a light surface), which covers the line / box / scatter forms used
# here. Print figures commit to the light surface only.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

parse_args_defaults <- list(
  summary_csv   = "summary_metrics.csv",
  replicate_csv = "all_replicate_metrics.csv",
  pathway_csv   = "all_pathway_stats.csv",
  hubreg_csv    = "hub_regression.csv",
  out_dir       = "figures",
  overlap_map   = "",   # e.g. "ov00=0,ov05=5,ov10=10" - condition -> overlap %
  nullhub_cond  = "ov05_nullhub",
  main_cond     = "ov05",
  width         = 7.0,
  height        = 4.4
)
source(file.path(Sys.getenv("SIM_SCRIPTS", unset = "."), "sim_common.R"))
args <- parse_args(parse_args_defaults)
log_header("08 figures")
dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

METHOD_LEVELS <- c("Original", "PS", "GSR")
METHOD_COLS <- c(Original = "#2a78d6", PS = "#eb6834", GSR = "#1baf7a")
INK        <- "#22221f"; INK_2 <- "#5b5a52"; GRIDC <- "#e6e5e0"

theme_sim <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      text             = element_text(colour = INK),
      axis.text        = element_text(colour = INK_2),
      axis.title       = element_text(colour = INK),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = GRIDC, linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      strip.text       = element_text(colour = INK, face = "bold", size = base_size - 1),
      legend.position  = "top",
      legend.title     = element_blank(),
      legend.key.size  = unit(0.9, "lines"),
      plot.title       = element_text(face = "bold", size = base_size + 1),
      plot.subtitle    = element_text(colour = INK_2, size = base_size - 1),
      plot.caption     = element_text(colour = INK_2, size = base_size - 2, hjust = 0)
    )
}
scale_method <- function(...) scale_colour_manual(values = METHOD_COLS, drop = FALSE, ...)
fill_method  <- function(...) scale_fill_manual(values = METHOD_COLS, drop = FALSE, ...)

save_fig <- function(p, name, w = args$width, h = args$height) {
  for (ext in c("png", "pdf")) {
    ggsave(file.path(args$out_dir, paste0(name, ".", ext)), p,
           width = w, height = h, dpi = 300, bg = "white")
  }
  cat("  wrote", name, "\n")
}

summ <- fread(args$summary_csv)
reps <- fread(args$replicate_csv)
pw   <- fread(args$pathway_csv)
hubreg <- if (file.exists(args$hubreg_csv)) fread(args$hubreg_csv) else data.table()

# Assign explicitly (not via a loop over a list) so the level order actually
# reaches ggplot -- facet and legend order depend on it.
summ[, method := factor(method, levels = METHOD_LEVELS)]
reps[, method := factor(method, levels = METHOD_LEVELS)]
pw[,   method := factor(method, levels = METHOD_LEVELS)]

# condition -> numeric overlap %, from the CLI map (falls back to parsing "ovNN")
ov <- if (nzchar(args$overlap_map)) {
  kv <- strsplit(strsplit(args$overlap_map, ",")[[1]], "=")
  data.table(condition = vapply(kv, `[`, character(1), 1),
             overlap_pct = as.numeric(vapply(kv, `[`, character(1), 2)))
} else {
  u <- unique(summ$condition)
  data.table(condition = u,
             overlap_pct = as.numeric(sub("^ov0*([0-9]+).*$", "\\1", u)))
}
summ <- merge(summ, ov, by = "condition", all.x = TRUE)
reps <- merge(reps, ov, by = "condition", all.x = TRUE)

is_main <- function(dt) dt[condition != args$nullhub_cond & !is.na(overlap_pct)]

# ---- Figure 1: null FPR vs overlap (the money plot) ---------------------
# Faceted by null class rather than encoding it with linetype: at 0% overlap the
# contaminated-null class is EMPTY by construction (no hub genes exist), and a
# facet shows that absence honestly instead of hiding it under the clean-null
# points. Falls back to a categorical x when only one overlap level is available.
d1 <- is_main(summ)[metric %in% c("fpr_contaminated_null", "fpr_clean_null") & !is.na(mean)]
if (nrow(d1) > 0) {
  d1[, null_class := factor(metric,
        levels = c("fpr_contaminated_null", "fpr_clean_null"),
        labels = c("Contaminated null\n(contains signal hubs)", "Clean null\n(no hub genes)"))]
  n_ov <- uniqueN(d1$overlap_pct)
  p1 <- ggplot(d1, aes(overlap_pct, mean, colour = method, group = method)) +
    geom_hline(yintercept = 0.05, linetype = "22", colour = INK_2, linewidth = 0.4) +
    { if (n_ov > 1) geom_line(linewidth = 0.8) else NULL } +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.6, linewidth = 0.5,
                  position = position_dodge(width = 1.1)) +
    geom_point(size = 2.4, position = position_dodge(width = 1.1)) +
    facet_wrap(~ null_class) +
    scale_method() +
    scale_x_continuous(breaks = sort(unique(d1$overlap_pct)),
                       labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = 0.25)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, NA)) +
    labs(title = "False positive rate among null pathways, by overlap level",
         subtitle = paste("Dashed line is the nominal 0.05 rate.",
                          "At 0% overlap no hub genes exist, so there is no contaminated-null class."),
         x = "Gene-slot overlap", y = "False positive rate",
         caption = "Mean over replicates; bars are 95% CI. Significance at p < 0.05.") +
    theme_sim() +
    theme(plot.subtitle = element_text(colour = INK_2, size = 9, lineheight = 1.15))
  save_fig(p1, "fig1_null_fpr_vs_overlap", w = 7.6, h = 4.4)
  fwrite(d1, file.path(args$out_dir, "fig1_data.csv"))
}

# ---- Figure 2: power by effect tier --------------------------------------
d2 <- is_main(reps)[metric %in% c("power_low", "power_med", "power_high") & !is.na(value)]
if (nrow(d2) > 0) {
  d2[, tier := factor(sub("power_", "", metric), levels = c("low", "med", "high"),
                      labels = c("Low", "Medium", "High"))]
  p2 <- ggplot(d2, aes(tier, value, fill = method, colour = method)) +
    geom_boxplot(alpha = 0.25, outlier.size = 0.5, linewidth = 0.45,
                 position = position_dodge(width = 0.78), width = 0.68) +
    scale_method(); p2 <- p2 + fill_method() +
    facet_wrap(~ paste0(overlap_pct, "% overlap")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(title = "Power by planted effect tier",
         subtitle = "Truly-enriched pathways called significant, by planted effect size",
         x = "Planted effect tier", y = "Power",
         caption = "Box = IQR over replicates; whiskers 1.5 x IQR.") +
    theme_sim()
  save_fig(p2, "fig2_power_by_tier")
  fwrite(d2, file.path(args$out_dir, "fig2_data.csv"))
}

# ---- Figure 3: ranking accuracy vs overlap -------------------------------
d3 <- is_main(summ)[metric %in% c("auprc_ses", "auprc_neglog10p", "spearman_truth_ses") &
                      !is.na(mean)]
if (nrow(d3) > 0) {
  d3[, panel := factor(metric,
      levels = c("auprc_ses", "auprc_neglog10p", "spearman_truth_ses"),
      labels = c("AUPRC (effect size)", "AUPRC (-log10 p)", "Spearman vs planted truth"))]
  p3 <- ggplot(d3, aes(overlap_pct, mean, colour = method, group = method)) +
    geom_line(linewidth = 0.8) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.25, linewidth = 0.5) +
    geom_point(size = 2.4) +
    facet_wrap(~ panel, scales = "free_y") +
    scale_method() +
    scale_x_continuous(breaks = sort(unique(d3$overlap_pct)),
                       labels = function(x) paste0(x, "%")) +
    labs(title = "Ranking accuracy against the planted truth",
         subtitle = "Discriminating the 30 truly-enriched pathways from all null pathways",
         x = "Gene-slot overlap", y = "Accuracy",
         caption = "Original has no null distribution, so it is ranked on -log10 p only.") +
    theme_sim()
  save_fig(p3, "fig3_ranking_accuracy_vs_overlap", w = 9.0)
  fwrite(d3, file.path(args$out_dir, "fig3_data.csv"))
}

# ---- Figure 4: mechanism -- -log10(p) vs number of signal hubs -----------
d4 <- pw[is_truly_enriched == FALSE & condition == args$main_cond]
if (nrow(d4) > 0) {
  lab <- hubreg[condition == args$main_cond]
  lab[, txt := sprintf("slope = %+.3f\n(p = %.2g)", slope, p)]
  p4 <- ggplot(d4, aes(n_signal_hubs, neglog10p, colour = method)) +
    geom_jitter(width = 0.18, height = 0, alpha = 0.18, size = 0.7, stroke = 0) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.9) +
    facet_wrap(~ method) +
    scale_method(guide = "none") +
    labs(title = "Pathway significance against hub contamination",
         subtitle = "Null pathways only, so hub genes are the sole source of signal",
         x = "Number of signal-carrying hub genes in the pathway",
         y = expression(-log[10](p))) +
    theme_sim()
  if (nrow(lab) > 0) {
    p4 <- p4 + geom_text(data = lab, aes(label = txt), x = Inf, y = Inf,
                         hjust = 1.05, vjust = 1.2, size = 3, colour = INK,
                         inherit.aes = FALSE)
    p4 <- p4 + facet_wrap(~ method)
  }
  save_fig(p4, "fig4_mechanism_hub_slope", w = 8.5)
  fwrite(d4, file.path(args$out_dir, "fig4_data.csv"))
  fwrite(hubreg, file.path(args$out_dir, "fig4_slopes.csv"))
}

# ---- Figure 5: null-hub control ------------------------------------------
d5 <- reps[condition %in% c(args$main_cond, args$nullhub_cond) &
             metric == "fpr_contaminated_null" & !is.na(value)]
if (nrow(d5) > 0) {
  d5[, cond_lab := factor(condition, levels = c(args$main_cond, args$nullhub_cond),
        labels = c("Hubs carry signal", "Null-hub control\n(hubs carry no signal)"))]
  p5 <- ggplot(d5, aes(cond_lab, value, fill = method, colour = method)) +
    geom_hline(yintercept = 0.05, linetype = "22", colour = INK_2, linewidth = 0.4) +
    geom_boxplot(alpha = 0.25, outlier.size = 0.5, linewidth = 0.45,
                 position = position_dodge(width = 0.78), width = 0.62) +
    scale_method(); p5 <- p5 + fill_method() +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = "Null-hub control: contaminated-null FPR with and without hub signal",
         subtitle = "Same 5% architecture; only the hub genes' signal is switched off",
         x = NULL, y = "Contaminated-null FPR",
         caption = "Dashed line is the nominal 0.05 rate.") +
    theme_sim()
  save_fig(p5, "fig5_null_hub_control", w = 6.5)
  fwrite(d5, file.path(args$out_dir, "fig5_data.csv"))
}

cat("\n08 figures: done ->", args$out_dir, "\n")
