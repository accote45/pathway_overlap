library(data.table)
library(dplyr)
library(jsonlite)
library(GSA)
library(ggplot2)

# Usage:
# Rscript delta_rank_ot_correlation.R <trait> <tool_base> <birewire_results_file> <gmt_file> <opentargets_json_dir> [top_ns_csv]
#

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Args: trait tool_base birewire_results_file gmt_file opentargets_json_dir [top_ns_csv]")
}
trait <- args[[1]]
tool_base <- args[[2]]
birewire_results_file <- args[[3]]
gmt_file <- args[[4]]
ot_dir <- args[[5]]
top_ns <- if (length(args) >= 6) as.integer(strsplit(args[[6]], ",")[[1]]) else c(50,100,250,500)

trait <- "ad"
tool_base <- "magma"
birewire_results_file <- "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/empirical_pvalues/magma_birewire/ad/ad_magma_birewire_empirical_pvalues.txt"
gmt_file <- "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/pathway_db/msigdb/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt"
ot_dir <- "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/results/drugtarget_test/associationByDatatypeDirect"
top_ns <- c(50,100,250,500)


# Map traits to OpenTargets disease IDs (extend as needed)
trait_mapping <- list(
  "t2d" = "MONDO_0005148",
  "cad" = "EFO_0001645",
  "ad" = "MONDO_0004975",
  "mdd" = "MONDO_0002050",
  "scz" = "MONDO_0005090",
  "ibd" = "EFO_0000555",
  "breast" = "MONDO_0007254"
)
disease_id <- trait_mapping[[trait]]
if (is.null(disease_id)) stop(sprintf("No OpenTargets ID mapping for trait '%s'", trait))

message("Trait: ", trait, "  Tool: ", tool_base)
message("Disease ID: ", disease_id)
message("Reading GMT: ", gmt_file)

# 1) Read GMT to build gene-to-pathway map
dat_gmt <- GSA.read.gmt(gmt_file)
path.list <- dat_gmt$genesets
names(path.list) <- dat_gmt$geneset.names
genes_long <- rbindlist(lapply(names(path.list), function(nm) {
  data.table(targetId = path.list[[nm]], pathway = nm)
}), use.names = TRUE, fill = TRUE)

# 2) Read OpenTargets JSON (associationByDatatypeDirect/*.json), keep disease and selected datatypes
json_files <- list.files(ot_dir, pattern = "\\.json$", full.names = TRUE)
if (length(json_files) == 0) stop("No JSON files found in opentargets_json_dir: ", ot_dir)

message("Reading OpenTargets JSON (this may take a bit)...")
ot_list <- lapply(json_files, function(f) {
  tryCatch({
    dd <- stream_in(file(f), verbose = FALSE)
    as.data.frame(dd)
  }, error = function(e) NULL)
})
ot_df <- bind_rows(Filter(Negate(is.null), ot_list))
if (nrow(ot_df) == 0) stop("Failed to read any OpenTargets JSON rows")

ot_df <- ot_df %>%
  filter(diseaseId == disease_id & datatypeId %in% c("animal_model","known_drug","literature")) %>%
  select(targetId, diseaseId, datatypeId, score, evidenceCount)

if (nrow(ot_df) == 0) stop("No OT rows remain after filtering to disease and datatypes")

# Pick the record per target with max evidenceCount
disease_targets <- ot_df %>%
  group_by(targetId) %>%
  slice_max(evidenceCount, with_ties = FALSE) %>%
  ungroup()

# 3) Aggregate to pathway-level scores
master <- genes_long %>% left_join(disease_targets, by = "targetId") %>% mutate(score = ifelse(is.na(score), 0, score))
pathway_scores <- master %>%
  group_by(pathway) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    max_score = max(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    num_genes = n(),
    num_with_evidence = sum(score > 0, na.rm = TRUE),
    evidence_density = num_with_evidence / pmax(num_genes, 1)
  ) %>%
  ungroup()

# 4) Read empirical results (BireWire file contains p_value, beta_value, empirical_pval, std_effect_size)
BW <- fread(birewire_results_file)
# Harmonize pathway column name
pcol <- intersect(names(BW), c("pathway_name","name","pathway","Pathway"))
if (length(pcol) == 0) stop("Could not find pathway name column in birewire_results_file")
setnames(BW, pcol[1], "pathway")

# Robust numeric coercions
to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
for (nm in intersect(c("empirical_pval","std_effect_size","p_value","beta_value"), names(BW))) {
  BW[[nm]] <- to_num(BW[[nm]])
}

# Sanity checks
need_cols_bw <- c("empirical_pval","std_effect_size","p_value","beta_value")
missing <- setdiff(need_cols_bw, names(BW))
if (length(missing) > 0) {
  stop("Missing required columns in birewire_results_file: ", paste(missing, collapse = ", "))
}

# 5) Create ranks
# BireWire_EmpPvalStdBeta: smaller empirical_pval better; tie-break by larger std_effect_size
BW_ord <- BW[order(empirical_pval, -std_effect_size)]
BW_ord[, rank_birewire := seq_len(.N)]
BW_ranks <- BW_ord[, .(pathway, rank_birewire)]

# PvalueBeta: smaller p_value better; tie-break by larger beta_value
PB_ord <- BW[order(p_value, -beta_value)]
PB_ord[, rank_pvaluebeta := seq_len(.N)]
PB_ranks <- PB_ord[, .(pathway, rank_pvaluebeta)]

# Merge ranks and compute delta
DR <- merge(BW_ranks, PB_ranks, by = "pathway")
DR[, delta_rank := rank_birewire - rank_pvaluebeta]  # negative => BireWire improved (moved up)

# 6) Merge with OT pathway scores
DR_scores <- merge(DR, as.data.table(pathway_scores), by.x = "pathway", by.y = "pathway", all.x = TRUE)

# Keep only pathways with OT scores
DR_scores <- DR_scores[is.finite(mean_score)]
if (nrow(DR_scores) == 0) stop("No overlapping pathways between ranks and OT scores")

# 7) Correlations: All + Top-N by baseline (PvalueBeta) to avoid selection on the target metric
subset_list <- list(`All Pathways` = copy(DR_scores))

# Build top-N subsets based on baseline rank (PvalueBeta)
for (N in sort(unique(top_ns))) {
  subset_list[[sprintf("Top %d Pathways", N)]] <- DR_scores[rank_pvaluebeta <= N]
}

# Helper to compute Spearman (primary) and Kendall (exploratory)
corr_one <- function(df, x, y, method = "spearman") {
  if (nrow(df) < 3) return(list(estimate = NA_real_, p = NA_real_, n = nrow(df)))
  ct <- suppressWarnings(cor.test(df[[x]], df[[y]], method = method, exact = FALSE))
  list(estimate = unname(ct$estimate), p = unname(ct$p.value), n = nrow(df))
}

summ <- rbindlist(lapply(names(subset_list), function(lbl) {
  d <- subset_list[[lbl]]
  sp_mean <- corr_one(d, "delta_rank", "mean_score", "spearman")
  sp_den  <- corr_one(d, "delta_rank", "evidence_density", "spearman")
  kd_mean <- corr_one(d, "delta_rank", "mean_score", "kendall")
  kd_den  <- corr_one(d, "delta_rank", "evidence_density", "kendall")
  data.table(
    trait = trait,
    tool_base = tool_base,
    subset = lbl,
    n_pathways = nrow(d),
    spearman_rho_mean = sp_mean$estimate,
    spearman_p_mean   = sp_mean$p,
    spearman_rho_density = sp_den$estimate,
    spearman_p_density   = sp_den$p,
    kendall_tau_mean   = kd_mean$estimate,
    kendall_p_mean     = kd_mean$p,
    kendall_tau_density = kd_den$estimate,
    kendall_p_density   = kd_den$p
  )
}), fill = TRUE)

# 8) Write outputs
out_prefix <- sprintf("%s_%s", trait, tool_base)
fwrite(DR_scores[order(rank_pvaluebeta)],
       sprintf("%s_delta_rank_with_OT_scores.csv", out_prefix))
fwrite(summ, sprintf("%s_delta_rank_ot_correlation_summary.csv", out_prefix))

# 9) Scatterplots: OT score vs delta rank
plot_dir <- file.path(getwd(), "figs_delta_rank_ot")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

for (lbl in names(subset_list)) {
  d <- subset_list[[lbl]]
  if (nrow(d) < 3) next

  # long format for two OT metrics
  d_long <- data.table::melt(
    d[, .(pathway, delta_rank, mean_score, evidence_density)],
    id.vars = c("pathway","delta_rank"),
    variable.name = "metric",
    value.name = "ot_score"
  )
  # pretty metric labels
  d_long[, metric := factor(metric,
                            levels = c("mean_score","evidence_density"),
                            labels = c("OT mean score","OT evidence density"))]

  # per-metric Spearman stats for annotation
  stats <- d_long[, {
    ct <- suppressWarnings(cor.test(delta_rank, ot_score, method = "spearman"))
    xr <- range(delta_rank, na.rm = TRUE); yr <- range(ot_score, na.rm = TRUE)
    data.table(
      rho = unname(ct$estimate),
      p = unname(ct$p.value),
      n = .N,
      x = xr[1] + 0.05 * diff(xr),
      y = yr[2] - 0.05 * diff(yr)
    )
  }, by = metric]
  stats[, label := sprintf("rho=%.2f, p=%.2g, n=%d", rho, p, n)]

  p <- ggplot(d_long, aes(x = delta_rank, y = ot_score)) +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.2, color = "grey70") +
    geom_point(alpha = 0.6, size = 1.4, color = "#2c7fb8") +
    geom_smooth(method = "loess", se = FALSE, color = "#d95f0e", linewidth = 0.8) +
    geom_text(data = stats, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.2) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    labs(
      title = sprintf("%s — OT score vs Δrank (Δ = BireWire − PvalueBeta)", trait),
      subtitle = lbl,
      x = "Delta rank (negative = improved by BireWire)",
      y = "OpenTargets score"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank()
    )

  out_pdf <- file.path(plot_dir, sprintf("%s_delta_rank_vs_OT_%s.pdf",
                                         out_prefix, gsub(" ", "_", tolower(lbl))))
  ggsave(out_pdf, p, width = 10, height = 4, device = cairo_pdf)
  message("Wrote: ", out_pdf)
}

message("Wrote: ", sprintf("%s_delta_rank_with_OT_scores.csv", out_prefix))
message("Wrote: ", sprintf("%s_delta_rank_ot_correlation_summary.csv", out_prefix))