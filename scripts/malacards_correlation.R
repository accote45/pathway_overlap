library(plyr)
  library(tidyverse)
  library(data.table)
  library(GSA)

write_results_csv <- function(data, trait, tool_base, prefix="", suffix="", verbose=TRUE) {
  filename <- paste0(trait, "_", tool_base, prefix, ".csv")
  write.csv(data, filename, row.names=FALSE)
  if (verbose) cat("Wrote", nrow(data), "rows to", filename, "\n")
  return(filename)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript malacards_correlation_stats.R <trait> <tool_base> <malacards_path> <birewire_results> <keeppathsize_results> <gmt_file>")
}

trait <- args[1]
tool_base <- args[2]
malacards_path <- args[3]
birewire_results_file <- args[4]
keeppathsize_results_file <- args[5]
gmt_file <- args[6]


# Top-N cutoffs for subset analyses (keep consistent with OT script)
top_ns <- c(100, 250, 500)

cat("======= Starting MalaCards correlation for", trait, "(", tool_base, ") =======\n")

# 1) Load pathway GMT and build gene->pathway long table
cat("Loading pathway data from", gmt_file, "...\n")
dat <- GSA.read.gmt(gmt_file)
path.list <- dat$genesets
names(path.list) <- dat$geneset.names

genes_long <- rbindlist(lapply(names(path.list), function(name) {
  # Pathway genes now stored under ENSEMBL (must already be ENSEMBL IDs in the GMT)
  data.table(ENSEMBL = toupper(path.list[[name]]), name = name)
}))

# 2) Read MalaCards CSVs (file or directory; combine if multiple)
is_dir <- function(p) isTRUE(file.info(p)$isdir)

list_malacards_files <- function(base_path, trait) {
  if (!file.exists(base_path)) stop("MalaCards path does not exist: ", base_path)
  if (!is_dir(base_path)) {
    # Single file path supplied: enforce *ensembl.csv requirement
    if (!grepl("ensembl\\.csv$", tolower(basename(base_path))))
      stop("Provided file does not match *ensembl.csv: ", base_path)
    return(normalizePath(base_path))
  }
  # Only list files ending with ensembl.csv
  files <- list.files(base_path, pattern = "ensembl\\.csv$", full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) stop("No *ensembl.csv files under: ", base_path)
  tci <- tolower(trait)
  sel <- grepl(sprintf("^malacards_%s.*ensembl\\.csv$", tci), tolower(basename(files))) |
         grepl(tci, tolower(basename(files)))
  out <- files[sel]
  if (length(out) == 0) stop("No *ensembl.csv MalaCards files matched trait '", trait,
                             "'. Available: ", paste(basename(files), collapse = ", "))
  normalizePath(out)
}

numify <- function(x) as.numeric(gsub("[^0-9eE\\+\\-\\.]", "", as.character(x)))

read_one_malacards <- function(f) {
  all <- readLines(f, warn = FALSE)
  # Require header line containing ENSEMBL and Score
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
cat("Reading", length(mc_files), "MalaCards file(s)...\n")
mc_list <- lapply(mc_files, read_one_malacards)
mc_all <- rbindlist(mc_list, use.names = TRUE, fill = TRUE)

# Collapse to one score per gene (max across files)
gene_scores <- mc_all %>%
  filter(!is.na(Score)) %>%
  group_by(ENSEMBL) %>%
  summarise(Score = max(Score, na.rm = TRUE), .groups = "drop")

cat("MalaCards rows:", nrow(mc_all), "| unique ENSEMBL genes:", nrow(gene_scores), "\n")

# Rank-normalize scores: r = ascending rank among positive scores (highest raw -> r = n)
# normalized Score = (r + 1) / (n + 1); genes without score -> 0
gene_scores$Score_original <- gene_scores$Score
positive <- !is.na(gene_scores$Score)
if (length(positive) > 0) {
  r <- rank(gene_scores$Score[positive], ties.method = "average")  # lowest=1, highest=n_pos
  gene_scores$Score <- 0
  gene_scores$Score[positive] <- (r + 1) / (length(positive) + 1)
  cat("Rank-normalized", length(positive), "genes; max normalized =", max(gene_scores$Score), "\n")
} else {
  gene_scores$Score <- 0
}


# 3) Build pathway-level MalaCards scores
masterfin <- genes_long %>%
  left_join(gene_scores, by = "ENSEMBL") %>%
  mutate(Score = ifelse(is.na(Score), 0, Score))

pathway_scores <- masterfin %>% 
  group_by(name) %>% 
  summarise(
    mean_score = mean(Score, na.rm = TRUE),
    num_genes = n(),
    num_with_evidence = sum(Score > 0, na.rm = TRUE),
    evidence_density = num_with_evidence / num_genes,
    max_score = max(Score, na.rm = TRUE),
    median_score = median(Score, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  as.data.frame()

# 4) Load empirical results files (BireWire, KeepPathSize)
birewire_data <- read.table(birewire_results_file, header = TRUE)
keeppath_data <- read.table(keeppathsize_results_file, header = TRUE)

cat("Loaded", nrow(birewire_data), "pathways from BireWire results\n")
cat("Loaded", nrow(keeppath_data), "pathways from KeepPathSize results\n")

# Ensure numeric fields
num_cols <- c("empirical_pval", "p_value", "beta_value", "std_effect_size")
for (nm in num_cols) {
  if (nm %in% names(birewire_data)) birewire_data[[nm]] <- as.numeric(as.character(birewire_data[[nm]]))
  if (nm %in% names(keeppath_data)) keeppath_data[[nm]] <- as.numeric(as.character(keeppath_data[[nm]]))
}

# 5) Rank Correlation Analysis (style kept close to OT script)
cat("\n======= Performing Rank Correlation Analysis (MalaCards) =======\n")
rank_correlation_results <- data.frame()

all_paths_with_scores <- pathway_scores %>% 
  filter(!is.na(mean_score)) %>%
  select(name, mean_score, evidence_density)

ranking_methods <- list(
  list(method_name = "PvalueBeta", 
       data = birewire_data, 
       rank_col = c("p_value", "beta_value"),
       sig_col = "p_value",
       higher_better = c(FALSE, TRUE)),
  list(method_name = "BireWire_EmpPvalStdBeta", 
       data = birewire_data, 
       rank_col = c("empirical_pval", "std_effect_size"),
       sig_col = "empirical_pval",
       higher_better = c(FALSE, TRUE)),
  list(method_name = "KeepPathSize_EmpPvalStdBeta", 
       data = keeppath_data, 
       rank_col = c("empirical_pval", "std_effect_size"),
       sig_col = "empirical_pval",
       higher_better = c(FALSE, TRUE))
)

for (ranking in ranking_methods) {
  method_name <- ranking$method_name
  data <- ranking$data
  rank_cols <- ranking$rank_col
  higher_better <- ranking$higher_better
  
  cat(paste("\nCalculating rank correlation for", method_name, "...\n"))
  
  ranking_data <- data
  
  # Build rank keys (invert if higher is better so smaller is better)
  for (i in seq_along(rank_cols)) {
    col <- rank_cols[i]
    if (!col %in% names(ranking_data)) next
    ranking_data[[paste0("ranking_", i)]] <- if (isTRUE(higher_better[i])) -ranking_data[[col]] else ranking_data[[col]]
  }
  
  # Order by 1 or 2 keys (as in OT script)
  if (length(rank_cols) == 1 && "ranking_1" %in% names(ranking_data)) {
    ranking_data <- ranking_data[order(ranking_data$ranking_1), ]
  } else if (length(rank_cols) >= 2 && all(c("ranking_1", "ranking_2") %in% names(ranking_data))) {
    ranking_data <- ranking_data[order(ranking_data$ranking_1, ranking_data$ranking_2), ]
  } else {
    # If required columns missing, skip
    cat("  Missing columns for", method_name, "- skipping.\n")
    next
  }
  
  ranking_data$pathway_rank <- seq_len(nrow(ranking_data))
  ranked_paths <- ranking_data %>% select(pathway_name, pathway_rank) %>% rename(name = pathway_name)
  
  # Merge with MalaCards pathway scores
  all_merged_ranks <- merge(ranked_paths, all_paths_with_scores, by = "name")
  
  # Subsets: All + Top-N
  subset_list <- list("All Pathways" = all_merged_ranks)
  for (N in sort(unique(top_ns))) {
    topN_ranked_paths <- ranked_paths %>% filter(pathway_rank <= N)
    subset_list[[sprintf("Top %d Pathways", N)]] <- merge(topN_ranked_paths, all_paths_with_scores, by = "name")
  }
  
  add_corr_row <- function(df, subset_label) {
    if (nrow(df) <= 1) return(NULL)
    sp_mean <- suppressWarnings(cor.test(df$pathway_rank, df$mean_score, method = "spearman"))
    sp_den  <- suppressWarnings(cor.test(df$pathway_rank, df$evidence_density, method = "spearman"))
    kd_mean <- suppressWarnings(cor.test(df$pathway_rank, df$mean_score, method = "kendall"))
    kd_den  <- suppressWarnings(cor.test(df$pathway_rank, df$evidence_density, method = "kendall"))
    data.frame(
      trait = trait,
      tool_base = tool_base,
      method = method_name,
      subset = subset_label,
      n_pathways = nrow(df),
      rank_mean_score_correlation = unname(sp_mean$estimate),
      rank_mean_score_pvalue = sp_mean$p.value,
      rank_evidence_density_correlation = unname(sp_den$estimate),
      rank_evidence_density_pvalue = sp_den$p.value,
      kendall_mean_score_tau = unname(kd_mean$estimate),
      kendall_mean_score_pvalue = kd_mean$p.value,
      kendall_evidence_density_tau = unname(kd_den$estimate),
      kendall_evidence_density_pvalue = kd_den$p.value,
      stringsAsFactors = FALSE
    )
  }
  
  for (lbl in names(subset_list)) {
    df_sub <- subset_list[[lbl]]
    if (nrow(df_sub) > 0) {
      cat(sprintf("  Found %d pathways for %s with MalaCards scores\n", nrow(df_sub), lbl))
      row <- add_corr_row(df_sub, lbl)
      if (!is.null(row)) rank_correlation_results <- rbind(rank_correlation_results, row)
    } else {
      cat("  No matching pathways found for", method_name, "- skipping", lbl, "\n")
    }
  }
}

# 6) Write correlation results
if (nrow(rank_correlation_results) > 0) {
  write_results_csv(rank_correlation_results, trait, tool_base, "_malacards_rank_correlation_summary")
  cat("======= Rank Correlation Analysis Complete =======\n")
} else {
  cat("No results to write.\n")
}