#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 01 - Build the simulated gene-by-pathway architecture (instructions Sec. 2)
#
# Emits, for one (condition, replicate):
#   sim_geneset.gmt     - pathway database in the main pipeline's GMT layout
#   ground_truth.tsv    - pathway_id, class, tier, n_hubs, n_signal_hubs, ...
#   gene_signal.tsv     - per-gene causal status and mean Z-shift mu
#   architecture_diag.txt - overlap / degree diagnostics for the log
#
# Design (defaults follow Sec. 2 exactly):
#   100 pathways x 50 genes  = 5000 gene-slots
#   30 truly enriched (10 low / 10 med / 10 high), 70 truly null
#   overlap_frac of slots are taken by hub genes of fixed degree; each hub sits
#   in `hub_slots_in_high` high-tier pathways and the rest in null pathways
#   drawn from a restricted "contaminable" subset, so that clean nulls survive.
# ---------------------------------------------------------------------------

parse_args_defaults <- list(
  gene_pool_file        = "",     # one eligible gene ID per line (from 02_define_ld_blocks)
  out_prefix            = "sim",
  seed                  = 1,

  n_pathways            = 100,
  pathway_size          = 50,
  variable_sizes        = FALSE,  # draw sizes from an MSigDB-like lognormal instead
  size_median           = 37,     # used only when variable_sizes=TRUE
  size_log_sd           = 0.95,   # "        "
  size_min              = 10,
  size_max              = 200,

  n_enriched_per_tier   = 10,     # -> 30 enriched pathways total
  n_contaminable_nulls  = 35,     # null pathways eligible to receive hubs

  overlap_frac          = 0.05,   # fraction of gene-slots held by hub genes
  hub_degree            = 10,     # pathways each hub gene belongs to
  hub_slots_in_high     = 1,      # of those, how many are high-tier enriched
  hub_signal            = TRUE,   # FALSE = the Sec.6 "null-hub" control

  # Per-tier mean Z-shift at the causal SNP. Calibrated by 10_calibrate.R;
  # these defaults are placeholders and MUST be overridden from calibration.
  mu_low                = 0.10,
  mu_med                = 0.20,
  mu_high               = 0.35,
  mu_hub                = -1,     # <0 means "use mu_high" (Sec. 3: mu_hub = mu_high)
  effect_scale          = 1.0,    # global multiplier s (small/medium/large axis)

  frac_causal           = 1.0     # fraction of an enriched pathway's genes made causal
)
source(file.path(Sys.getenv("SIM_SCRIPTS", unset = "."), "sim_common.R"))
args <- parse_args(parse_args_defaults)

set.seed(as.integer(args$seed))
log_header("01 build_architecture | seed=", args$seed)

if (!nzchar(args$gene_pool_file)) stop("gene_pool_file is required")
gene_pool <- unique(readLines(args$gene_pool_file, warn = FALSE))
gene_pool <- gene_pool[nzchar(gene_pool)]
cat("Eligible gene pool:", length(gene_pool), "genes\n")

n_path <- as.integer(args$n_pathways)

# ---- pathway sizes -------------------------------------------------------
if (isTRUE(args$variable_sizes)) {
  sizes <- round(rlnorm(n_path, meanlog = log(args$size_median), sdlog = args$size_log_sd))
  sizes <- pmin(pmax(sizes, args$size_min), args$size_max)
} else {
  sizes <- rep(as.integer(args$pathway_size), n_path)
}
n_slots <- sum(sizes)

# ---- hub budget ----------------------------------------------------------
hub_degree <- as.integer(args$hub_degree)
n_hub_slots <- round(args$overlap_frac * n_slots)
n_hubs <- if (n_hub_slots <= 0) 0L else as.integer(round(n_hub_slots / hub_degree))
n_hub_slots <- n_hubs * hub_degree           # exact, after rounding to whole hubs
n_unique_slots <- n_slots - n_hub_slots

if (n_hubs > 0 && hub_degree > n_path) stop("hub_degree exceeds n_pathways")
if (length(gene_pool) < n_unique_slots + n_hubs) {
  stop("Gene pool too small: need ", n_unique_slots + n_hubs, ", have ", length(gene_pool))
}

# ---- pathway roles (assigned to a random permutation of pathway indices) --
pids <- sprintf("SIMPATH_%03d", seq_len(n_path))
perm <- sample(n_path)
k <- as.integer(args$n_enriched_per_tier)
idx_high <- perm[seq_len(k)]
idx_med  <- perm[k + seq_len(k)]
idx_low  <- perm[2 * k + seq_len(k)]
idx_null <- perm[(3 * k + 1):n_path]
n_contam_elig <- min(as.integer(args$n_contaminable_nulls), length(idx_null))
idx_contaminable <- idx_null[seq_len(n_contam_elig)]

tier <- rep("null", n_path)
tier[idx_high] <- "high"; tier[idx_med] <- "med"; tier[idx_low] <- "low"

# ---- draw genes ----------------------------------------------------------
drawn <- sample(gene_pool, n_unique_slots + n_hubs)
hub_genes <- if (n_hubs > 0) drawn[seq_len(n_hubs)] else character(0)
unique_genes <- drawn[n_hubs + seq_len(n_unique_slots)]

# ---- place hubs ----------------------------------------------------------
# Each hub: `hub_slots_in_high` slots in distinct high-tier pathways, the rest
# in distinct contaminable-null pathways. This is what manufactures the
# "contaminated null" class while leaving the remaining nulls clean.
membership <- vector("list", n_path)   # pathway index -> gene ids
hub_members <- vector("list", n_path)  # pathway index -> hub gene ids
for (i in seq_len(n_path)) { membership[[i]] <- character(0); hub_members[[i]] <- character(0) }

n_in_high <- min(as.integer(args$hub_slots_in_high), length(idx_high))
for (h in seq_len(n_hubs)) {
  g <- hub_genes[h]
  targets <- c(
    sample(idx_high, n_in_high),
    sample(idx_contaminable, hub_degree - n_in_high)
  )
  for (i in targets) {
    membership[[i]] <- c(membership[[i]], g)
    hub_members[[i]] <- c(hub_members[[i]], g)
  }
}

# ---- fill the rest with unique (degree-1) genes --------------------------
# Assign in a random pathway order so no pathway systematically gets leftovers.
cursor <- 0L
for (i in sample(n_path)) {
  need <- sizes[i] - length(membership[[i]])
  if (need < 0) {
    stop("Pathway ", pids[i], " received ", length(membership[[i]]),
         " hubs but has size ", sizes[i], " - lower hub_degree or overlap_frac")
  }
  if (need > 0) {
    membership[[i]] <- c(membership[[i]], unique_genes[cursor + seq_len(need)])
    cursor <- cursor + need
  }
}
stopifnot(cursor == n_unique_slots)
names(membership) <- pids

# ---- signal assignment (Sec. 3) ------------------------------------------
s <- as.numeric(args$effect_scale)
mu_by_tier <- c(low = args$mu_low, med = args$mu_med, high = args$mu_high) * s
mu_hub_val <- (if (args$mu_hub < 0) args$mu_high else args$mu_hub) * s

# Unique genes inherit the tier of their (single) pathway; a fraction
# `frac_causal` of an enriched pathway's own genes are made causal.
gene_tier <- setNames(rep("null", length(unique_genes)), unique_genes)
for (i in seq_len(n_path)) {
  if (tier[i] == "null") next
  own <- setdiff(membership[[i]], hub_genes)
  n_causal <- round(args$frac_causal * length(own))
  if (n_causal > 0) gene_tier[sample(own, n_causal)] <- tier[i]
}

gene_mu <- ifelse(gene_tier == "null", 0, mu_by_tier[gene_tier])
names(gene_mu) <- names(gene_tier)

hub_mu <- setNames(rep(if (isTRUE(args$hub_signal)) mu_hub_val else 0, n_hubs), hub_genes)
hub_tier <- setNames(rep(if (isTRUE(args$hub_signal)) "hub" else "null", n_hubs), hub_genes)

all_mu <- c(gene_mu, hub_mu)
all_tier <- c(gene_tier, hub_tier)

gene_degree <- table(unlist(membership))

gene_signal <- data.table(
  gene_id = names(all_mu),
  is_hub = names(all_mu) %in% hub_genes,
  tier = unname(all_tier[names(all_mu)]),
  mu = unname(all_mu),
  degree = as.integer(gene_degree[names(all_mu)])
)
gene_signal[, causal := mu > 0]

# ---- ground truth table (Sec. 2) -----------------------------------------
n_hubs_in <- vapply(hub_members, length, integer(1))
n_signal_hubs_in <- vapply(hub_members, function(gs) sum(all_mu[gs] > 0), integer(1))
# `contaminated_null` is defined by HUB PRESENCE, not by hub signal, so that the
# null-hub control (hub_signal=FALSE) measures FPR in the same pathway class.
pclass <- ifelse(tier != "null", paste0("enriched_", tier),
                 ifelse(n_hubs_in > 0, "contaminated_null", "clean_null"))

planted_effect <- vapply(seq_len(n_path), function(i) mean(all_mu[membership[[i]]]), numeric(1))
truth_rank <- c(null = 0, low = 1, med = 2, high = 3)[tier]

ground_truth <- data.table(
  pathway_id = pids,
  class = pclass,
  tier = tier,
  size = sizes,
  n_hubs = n_hubs_in,
  n_signal_hubs = n_signal_hubs_in,
  n_causal_genes = vapply(seq_len(n_path), function(i) sum(all_mu[membership[[i]]] > 0), integer(1)),
  planted_effect = planted_effect,
  truth_rank = unname(truth_rank),
  is_truly_enriched = tier != "null"
)

# ---- write ---------------------------------------------------------------
write_gmt(membership, paste0(args$out_prefix, "_geneset.gmt"))
fwrite(ground_truth, paste0(args$out_prefix, "_ground_truth.tsv"), sep = "\t")
fwrite(gene_signal, paste0(args$out_prefix, "_gene_signal.tsv"), sep = "\t")

# ---- diagnostics ---------------------------------------------------------
diag_lines <- c(
  sprintf("seed                     : %s", args$seed),
  sprintf("pathways / slots         : %d / %d", n_path, n_slots),
  sprintf("distinct genes           : %d", length(all_mu)),
  sprintf("hub genes / degree       : %d / %d", n_hubs, hub_degree),
  sprintf("hub slots (overlap frac) : %d (%.4f requested %.4f)",
          n_hub_slots, n_hub_slots / n_slots, args$overlap_frac),
  sprintf("hub degree as %% of paths : %.2f%%", 100 * hub_degree / n_path),
  sprintf("hub_signal               : %s", args$hub_signal),
  sprintf("mu low/med/high/hub      : %.4f / %.4f / %.4f / %.4f",
          mu_by_tier["low"], mu_by_tier["med"], mu_by_tier["high"],
          if (isTRUE(args$hub_signal)) mu_hub_val else 0),
  sprintf("effect_scale s           : %.3f", s),
  sprintf("frac_causal              : %.3f", args$frac_causal),
  sprintf("causal genes / all genes : %d / %d (%.1f%%)",
          sum(all_mu > 0), length(all_mu), 100 * mean(all_mu > 0)),
  "",
  "pathway class counts:",
  paste0("  ", names(table(pclass)), " = ", as.integer(table(pclass)), collapse = "\n"),
  "",
  "signal hubs per null pathway:",
  paste0("  ", paste(capture.output(print(table(n_signal_hubs_in[tier == "null"]))), collapse = "\n  "))
)
writeLines(diag_lines, paste0(args$out_prefix, "_architecture_diag.txt"))
cat(paste(diag_lines, collapse = "\n"), "\n")
cat("\n01 build_architecture: done\n")
