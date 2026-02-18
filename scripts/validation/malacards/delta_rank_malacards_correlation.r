# Delta-rank vs MalaCards correlation (BireWire − PvalueBeta)
# Usage:
# Rscript delta_rank_malacards_correlation.R <trait> <tool_base> <birewire_results_file> <gmt_file> <malacards_path_or_file> [top_ns_csv]
#
# trait: e.g., cad
# tool_base: magma|prset
# birewire_results_file: path to file with columns {pathway_name|name|pathway|Pathway}, empirical_pval, std_effect_size, p_value, beta_value
# gmt_file: pathway GMT with ENSEMBL gene IDs
# malacards_path_or_file: a directory (with *ensembl.csv) or a single csv containing ENSEMBL and Score columns
# top_ns_csv: optional e.g. "50,100,250,500"

library(data.table)
library(dplyr)
library(GSA)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Args: trait tool_base birewire_results_file gmt_file malacards_path_or_file [top_ns_csv]")
}
trait <- args[[1]]
tool_base <- args[[2]]
birewire_results_file <- args[[3]]
gmt_file <- args[[4]]
malacards_path <- args[[5]]
top_ns <- if (length(args) >= 6) as.integer(strsplit(args[[6]], ",")[[1]]) else c(50,100,250,500)

message("Trait: ", trait, "  Tool: ", tool_base)
message("Reading GMT: ", gmt_file)

# 1) Read GMT (expects ENSEMBL IDs)
dat_gmt <- GSA.read.gmt(gmt_file)
path.list <- dat_gmt$genesets
names(path.list) <- dat_gmt$geneset.names
genes_long <- rbindlist(lapply(names(path.list), function(nm) {
  data.table(ENSEMBL = toupper(path.list[[nm]]), pathway = nm)
}), use.names = TRUE, fill = TRUE)

# 2) Read MalaCards CSVs (robust ENSEMBL+Score parser)
is_dir <- function(p) isTRUE(file.info(p)$isdir)
numify <- function(x) as.numeric(gsub("[^0-9eE\\+\\-\\.]", "", as.character(x)))

list_malacards_files <- function(base_path, trait) {
  if (!file.exists(base_path)) stop("MalaCards path does not exist: ", base_path)
  if (isTRUE(file.info(base_path)$isdir) == FALSE) {
    # Single file: must contain ENSEMBL + Score columns
    return(normalizePath(base_path))
  }
  files <- list.files(base_path, pattern = "ensembl\\.csv$", full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) stop("No *ensembl.csv files under: ", base_path)

  tci <- tolower(trait)
  bns <- tolower(basename(files))

  # 1) Allow trait followed by letters/digits/._- before 'ensembl.csv' (matches cad1, cad_2, cad-3, etc.)
  strict <- grepl(sprintf("^malacards_%s[a-z0-9._-]*ensembl\\.csv$", tci), bns, perl = TRUE)
  candidates <- files[strict]

  # 2) Token-based match, splitting also at letter–digit boundaries
  if (!length(candidates)) {
    noext <- sub("\\.csv$", "", bns)
    tokens <- strsplit(noext, "[^a-z0-9]+|(?<=[a-z])(?=[0-9])|(?<=[0-9])(?=[a-z])", perl = TRUE)
    has_trait_token <- vapply(tokens, function(v) tci %in% v, logical(1))
    candidates <- files[has_trait_token]
  }

  # 3) Substring fallback
  if (!length(candidates)) {
    candidates <- files[grepl(tci, bns, fixed = TRUE)]
  }

  if (!length(candidates)) {
    stop("No MalaCards files matched trait '", trait, "'. Available: ",
         paste(basename(files), collapse = ", "))
  }
  normalizePath(candidates)
}

read_one_malacards <- function(f) {
  all <- readLines(f, warn = FALSE)
  header_idx <- which(grepl("\\bENSEMBL\\b", all, ignore.case = TRUE) &
                        grepl("\\bScore\\b", all, ignore.case = TRUE))[1]
  if (is.na(header_idx)) stop("Could not find header with 'ENSEMBL' and 'Score' in: ", f)
  dt <- tryCatch(
    suppressWarnings(fread(f, skip = header_idx - 1, header = TRUE)),
    error = function(e) suppressWarnings(as.data.table(read.csv(f, skip = header_idx - 1, header = TRUE, check.names = FALSE)))
  )
  setnames(dt, old = names(dt), new = make.names(names(dt)))
  ensembl_col <- names(dt)[tolower(names(dt)) == "ensembl"]
  if (!length(ensembl_col)) ensembl_col <- names(dt)[grepl("ensembl", names(dt), ignore.case = TRUE)]
  score_col <- names(dt)[tolower(names(dt)) == "score"]
  if (!length(ensembl_col) || !length(score_col)) {
    stop("Missing required ENSEMBL/Score columns in: ", f, " | Columns: ", paste(names(dt), collapse = ", "))
  }
  out <- dt[, .(ENSEMBL = toupper(get(ensembl_col[1])), Score = numify(get(score_col[1])))]
  out <- out[!is.na(ENSEMBL) & nzchar(ENSEMBL)]
  out
}

mc_files <- list_malacards_files(malacards_path, trait)
message("Reading ", length(mc_files), " MalaCards file(s)")
mc_list <- lapply(mc_files, read_one_malacards)
mc_all <- rbindlist(mc_list, use.names = TRUE, fill = TRUE)

# Collapse to one score per gene (max across files)
gene_scores <- mc_all %>%
  filter(!is.na(Score)) %>%
  group_by(ENSEMBL) %>%
  summarise(Score = max(Score, na.rm = TRUE), .groups = "drop")

message("MalaCards rows: ", nrow(mc_all), " | unique ENSEMBL genes: ", nrow(gene_scores))

# Rank-normalize positive scores to (0,1); missing/zero -> 0
gene_scores$Score_original <- gene_scores$Score
positive <- !is.na(gene_scores$Score)
if (length(positive) > 0) {
  r <- rank(gene_scores$Score[positive], ties.method = "average")
  gene_scores$Score <- 0
  gene_scores$Score[positive] <- (r + 1) / (length(positive) + 1)
} else {
  gene_scores$Score <- 0
}

# 3) Aggregate to pathway-level MalaCards scores
master <- genes_long %>%
  left_join(gene_scores, by = "ENSEMBL") %>%
  mutate(Score = ifelse(is.na(Score), 0, Score))
pathway_scores <- master %>%
  group_by(pathway) %>%
  summarise(
    mean_score = mean(Score, na.rm = TRUE),
    max_score = max(Score, na.rm = TRUE),
    median_score = median(Score, na.rm = TRUE),
    num_genes = n(),
    num_with_evidence = sum(Score > 0, na.rm = TRUE),
    evidence_density = num_with_evidence / pmax(num_genes, 1),
    .groups = "drop"
  )

# 4) Read empirical results (BireWire file contains p_value, beta_value, empirical_pval, std_effect_size)
BW <- fread(birewire_results_file)
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

# Merge ranks and compute delta (negative => BireWire improved)
DR <- merge(BW_ranks, PB_ranks, by = "pathway")
DR[, delta_rank := rank_birewire - rank_pvaluebeta]

# 6) Merge with MalaCards pathway scores
DR_scores <- merge(DR, as.data.table(pathway_scores), by = "pathway", all.x = TRUE)

# Keep only pathways with MalaCards scores
DR_scores <- DR_scores[is.finite(mean_score)]
if (nrow(DR_scores) == 0) stop("No overlapping pathways between ranks and MalaCards scores")

# 7) Correlations: All + Top-N by baseline (PvalueBeta)
subset_list <- list(`All Pathways` = copy(DR_scores))
for (N in sort(unique(top_ns))) {
  subset_list[[sprintf("Top %d Pathways", N)]] <- DR_scores[rank_pvaluebeta <= N]
}

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
       sprintf("%s_delta_rank_with_MC_scores.csv", out_prefix))
fwrite(summ, sprintf("%s_delta_rank_malacards_correlation_summary.csv", out_prefix))

# 9) Scatterplots: MalaCards score vs delta rank
plot_dir <- file.path(getwd(), "figs_delta_rank_malacards")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

for (lbl in names(subset_list)) {
  d <- subset_list[[lbl]]
  if (nrow(d) < 3) next

  d_long <- data.table::melt(
    d[, .(pathway, delta_rank, mean_score, evidence_density)],
    id.vars = c("pathway","delta_rank"),
    variable.name = "metric",
    value.name = "mc_score"
  )
  d_long[, metric := factor(metric,
                            levels = c("mean_score","evidence_density"),
                            labels = c("MalaCards mean score","MalaCards evidence density"))]

  stats <- d_long[, {
    ct <- suppressWarnings(cor.test(delta_rank, mc_score, method = "spearman"))
    xr <- range(delta_rank, na.rm = TRUE); yr <- range(mc_score, na.rm = TRUE)
    data.table(
      rho = unname(ct$estimate),
      p = unname(ct$p.value),
      n = .N,
      x = xr[1] + 0.05 * diff(xr),
      y = yr[2] - 0.05 * diff(yr)
    )
  }, by = metric]
  stats[, label := sprintf("rho=%.2f, p=%.2g, n=%d", rho, p, n)]

  p <- ggplot(d_long, aes(x = delta_rank, y = mc_score)) +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.2, color = "grey70") +
    geom_point(alpha = 0.6, size = 1.4, color = "#2c7fb8") +
    geom_smooth(method = "loess", se = FALSE, color = "#d95f0e", linewidth = 0.8) +
    geom_text(data = stats, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.2) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    labs(
      title = sprintf("%s — MalaCards score vs Δrank (Δ = BireWire − PvalueBeta)", trait),
      subtitle = lbl,
      x = "Delta rank (negative = improved by BireWire)",
      y = "MalaCards score"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  out_pdf <- file.path(plot_dir, sprintf("%s_delta_rank_vs_MalaCards_%s.pdf",
                                         out_prefix, gsub(" ", "_", tolower(lbl))))
  ggsave(out_pdf, p, width = 10, height = 4, device = cairo_pdf)
  message("Wrote: ", out_pdf)
}

message("Wrote: ", sprintf("%s_delta_rank_with_MC_scores.csv", out_prefix))
message("Wrote: ", sprintf("%s_delta_rank_malacards_correlation_summary.csv", out_prefix))
