
library(data.table)
library(dplyr)
library(GSA)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Args: trait tool_base birewire_results_file scores_file [top_ns_csv]")
}
trait <- args[[1]]
tool_base <- args[[2]]
birewire_results_file <- args[[3]]
scores_file <- args[[4]]
top_ns <- if (length(args) >= 5) as.integer(strsplit(args[[5]], ",")[[1]]) else c(50,100,250,500)

message("Trait: ", trait, "  Tool: ", tool_base)

# 1. Load dorothea pathway scores
pathway_scores <- read.csv(scores_file)

# 3) Read empirical results (BireWire file contains p_value, beta_value, empirical_pval, std_effect_size)
BW <- fread(birewire_results_file)
pcol <- intersect(names(BW), c("pathway_name","name","pathway","Pathway"))
if (length(pcol) == 0) stop("Could not find pathway name column in birewire_results_file")
setnames(BW, pcol[1], "pathway")

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
for (nm in intersect(c("empirical_pval","std_effect_size","p_value","beta_value"), names(BW))) {
  BW[[nm]] <- to_num(BW[[nm]])
}
need_cols_bw <- c("empirical_pval","std_effect_size","p_value","beta_value")
missing <- setdiff(need_cols_bw, names(BW))
if (length(missing) > 0) {
  stop("Missing required columns in birewire_results_file: ", paste(missing, collapse = ", "))
}

# 4) Create ranks
BW_ord <- BW[order(empirical_pval, -std_effect_size)]
BW_ord[, rank_birewire := seq_len(.N)]
BW_ranks <- BW_ord[, .(pathway, rank_birewire)]

PB_ord <- BW[order(p_value, -beta_value)]
PB_ord[, rank_pvaluebeta := seq_len(.N)]
PB_ranks <- PB_ord[, .(pathway, rank_pvaluebeta)]

DR <- merge(BW_ranks, PB_ranks, by = "pathway")
DR[, delta_rank := rank_birewire - rank_pvaluebeta]

# 5) Merge with tissue specificity scores
DR_scores <- merge(DR, pathway_scores, by="pathway", all.x = TRUE)

if (nrow(DR_scores) == 0) stop("No overlapping pathways between ranks and tissue scores")


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
  sp_mean <- corr_one(d, "delta_rank", "score", "spearman")
  kd_mean <- corr_one(d, "delta_rank", "score", "kendall")
  data.table(
    trait = trait,
    tool_base = tool_base,
    subset = lbl,
    n_pathways = nrow(d),
    spearman_rho_mean = sp_mean$estimate,
    spearman_p_mean   = sp_mean$p,
    kendall_tau_mean   = kd_mean$estimate,
    kendall_p_mean     = kd_mean$p,
  )
}), fill = TRUE)

# 8) Write outputs
out_prefix <- sprintf("%s_%s", trait, tool_base)
fwrite(DR_scores[order(rank_pvaluebeta)],
       sprintf("%s_delta_rank_with_OT_scores.csv", out_prefix))
fwrite(summ, sprintf("%s_delta_rank_ot_correlation_summary.csv", out_prefix))

message("Wrote: ", sprintf("%s_delta_rank_with_OT_scores.csv", out_prefix))
message("Wrote: ", sprintf("%s_delta_rank_ot_correlation_summary.csv", out_prefix))


