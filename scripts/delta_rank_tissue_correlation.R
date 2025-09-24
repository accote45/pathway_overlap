# Delta-rank vs Tissue Specificity correlation (BireWire − PvalueBeta)
# Usage:
# Rscript delta_rank_tissue_correlation.R <trait> <tool_base> <birewire_results_file> <gmt_file> <tissue_scores_file> [top_ns_csv]


library(data.table)
library(dplyr)
library(GSA)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Args: trait tool_base birewire_results_file gmt_file tissue_scores_file [top_ns_csv]")
}
trait <- args[[1]]
tool_base <- args[[2]]
birewire_results_file <- args[[3]]
gmt_file <- args[[4]]
tissue_scores_file <- args[[5]]
top_ns <- if (length(args) >= 6) as.integer(strsplit(args[[6]], ",")[[1]]) else c(50,100,250,500)


trait <- "cad"
tool_base <- "magma_birewire"
birewire_results_file <- "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/empirical_pvalues/magma_birewire/cad/cad_magma_birewire_empirical_pvalues.txt"
gmt_file <- "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/pathway_db/msigdb/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt"
tissue_scores_file <- "/sc/arion/projects/psychgen/cotea02_prset/judit_revisions/software/1kg_test/GeneExpressionLandscape/data/Exp_Spe_DataTables/specificity"
top_ns <- if (length(args) >= 6) as.integer(strsplit(args[[6]], ",")[[1]]) else c(50,100,250,500)


message("Trait: ", trait, "  Tool: ", tool_base)
message("Reading GMT: ", gmt_file)


# 1. Load pathway data
cat("Loading pathway data from", gmt_file, "...\n")
dat <- GSA.read.gmt(gmt_file)
path.list <- dat$genesets
names(path.list) <- dat$geneset.names

# Convert the path list to a data table
genes_long <- rbindlist(lapply(names(path.list), function(name) {
  data.table(value = path.list[[name]], name = name)
}))

# 2. Load tissue specificity data
cat("Loading tissue specificity data from", tissue_scores_file, "...\n")
ts <- read.csv(tissue_scores_file)
cat("Loaded tissue data with", nrow(ts), "genes and", ncol(ts) - 1, "tissues\n")

# 3. Calculate pathway-level tissue expression scores
cat("Calculating pathway-level tissue expression scores...\n")
masterfin <- merge(genes_long, ts, by.x="value", by.y="Name")

# Get all tissue columns (excluding gene identifiers)
tissue_cols <- setdiff(colnames(masterfin), c("value", "name"))

# Calculate pathway-level tissue expression metrics (mean only, no median)
pathway_tissue_scores <- masterfin %>% 
  group_by(name) %>% 
  summarise(
    # Calculate mean for each tissue
    across(all_of(tissue_cols), 
           list(
             mean = ~mean(.x, na.rm = TRUE)
           ),
           .names = "{.col}_mean"),
    # Count genes per pathway
    num_genes = n()
  ) %>%
  as.data.frame()

cat("Generated tissue expression profiles for", nrow(pathway_tissue_scores), "pathways\n")


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
DR_scores <- merge(DR, pathway_tissue_scores, by.x = "pathway", by.y = "name", all.x = TRUE)


# Keep only pathways with tissue scores
tissue_cols <- grep("_mean$", names(DR_scores), value = TRUE)
DR_scores <- DR_scores[rowSums(is.na(DR_scores[, ..tissue_cols])) < length(tissue_cols)]

if (nrow(DR_scores) == 0) stop("No overlapping pathways between ranks and tissue scores")

# 6) Correlations: All + Top-N by baseline (PvalueBeta)
subset_list <- list(`All Pathways` = copy(DR_scores))
for (N in sort(unique(top_ns))) {
  subset_list[[sprintf("Top %d Pathways", N)]] <- DR_scores[rank_pvaluebeta <= N]
}

# 7) For each tissue, compute Spearman/Kendall correlations
corr_one <- function(df, x, y, method = "spearman") {
  if (nrow(df) < 3) return(list(estimate = NA_real_, p = NA_real_, n = nrow(df)))
  ct <- suppressWarnings(cor.test(df[[x]], df[[y]], method = method, exact = FALSE))
  list(estimate = unname(ct$estimate), p = unname(ct$p.value), n = nrow(df))
}

summ_list <- list()
for (subset_lbl in names(subset_list)) {
  d <- subset_list[[subset_lbl]]
  for (tissue in tissue_cols) {
    sp <- corr_one(d, "delta_rank", tissue, "spearman")
    kd <- corr_one(d, "delta_rank", tissue, "kendall")
    summ_list[[length(summ_list)+1]] <- data.table(
      trait = trait,
      tool_base = tool_base,
      subset = subset_lbl,
      tissue_metric = tissue,
      n_pathways = nrow(d),
      spearman_rho = sp$estimate,
      spearman_p = sp$p,
      kendall_tau = kd$estimate,
      kendall_p = kd$p
    )
  }
}
summ <- rbindlist(summ_list, fill = TRUE)

# 8) Write outputs
out_prefix <- sprintf("%s_%s", trait, tool_base)
fwrite(DR_scores[order(rank_pvaluebeta)],
       sprintf("%s_delta_rank_with_tissue_scores.csv", out_prefix))
fwrite(summ, sprintf("%s_delta_rank_tissue_correlation_summary.csv", out_prefix))

message("Wrote: ", sprintf("%s_delta_rank_with_tissue_scores.csv", out_prefix))
message("Wrote: ", sprintf("%s_delta_rank_tissue_correlation_summary.csv", out_prefix))