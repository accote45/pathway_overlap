#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 09 - BiRewire randomization diagnostics on the SIMULATED matrix (Sec. 8,
#      and reviewer R2 comment 3).
#
# Shows the constrained randomization space is non-trivial for the simulated
# architecture, not just for real MSigDB:
#   a) birewire.analysis.bipartite - the Jaccard index must fall from 1 and
#      plateau at or before the analytic swap bound N.
#   b) birewire.similarity - two independently randomized databases must be as
#      dissimilar from each other as either is from the original. If they were
#      much more similar to each other, the sampler would be stuck in a corner
#      of the space.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(BiRewire); library(data.table); library(ggplot2)
})
parse_args_defaults <- list(
  gmt        = "",
  out_prefix = "birewire_diag",
  n_pairs    = 20,     # independent randomized pairs for the similarity check
  step       = 10,
  accuracy   = 0.00005,
  seed       = 1
)
source(file.path(Sys.getenv("SIM_SCRIPTS", unset = "."), "sim_common.R"))
args <- parse_args(parse_args_defaults)
set.seed(as.integer(args$seed))
log_header("09 birewire_diagnostics")

pl <- read_gmt(args$gmt)
edges <- gmt_to_edges(pl)
m <- table(edges$gene, edges$pathway)
storage.mode(m) <- "integer"
cat("Simulated incidence matrix:", nrow(m), "genes x", ncol(m), "pathways,",
    sum(m), "edges\n")

# ---- (a) convergence trajectory ------------------------------------------
an <- birewire.analysis.bipartite(m, step = as.integer(args$step),
                                  max.iter = "n", accuracy = args$accuracy,
                                  verbose = FALSE, display = FALSE)
# `$data` is the Jaccard trajectory (a vector, or a matrix over repeats);
# `$N` is BiRewire's analytic bound on the number of successful swaps.
traj <- if (is.matrix(an$data)) colMeans(an$data) else as.numeric(an$data)
iters <- as.integer(names(an$data))
if (length(iters) != length(traj) || anyNA(iters)) iters <- seq_along(traj) * as.integer(args$step)
bound <- if (!is.null(an$N)) as.numeric(an$N) else NA_real_

traj_dt <- data.table(swaps = iters, jaccard = traj)
fwrite(traj_dt, paste0(args$out_prefix, "_convergence.csv"))

plateau <- median(tail(traj, max(3, ceiling(length(traj) / 5))))
reached <- which(traj <= plateau * 1.01)[1]
swaps_to_plateau <- if (is.na(reached)) NA_real_ else iters[reached]

p <- ggplot(traj_dt, aes(swaps, jaccard)) +
  geom_line(colour = "#2a78d6", linewidth = 0.8) +
  { if (is.finite(bound)) geom_vline(xintercept = bound, linetype = "22",
                                     colour = "#5b5a52", linewidth = 0.4) else NULL } +
  { if (is.finite(bound)) annotate("text", x = bound, y = max(traj), hjust = -0.05,
                                   vjust = 1, size = 3, colour = "#5b5a52",
                                   label = "analytic swap bound N") else NULL } +
  labs(title = "BiRewire convergence on the simulated gene-by-pathway matrix",
       subtitle = "Jaccard similarity to the original, as successful swaps accumulate",
       x = "Successful swaps", y = "Jaccard index vs. original") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
for (ext in c("png", "pdf")) {
  ggsave(paste0(args$out_prefix, "_convergence.", ext), p,
         width = 6.5, height = 4, dpi = 300, bg = "white")
}

# ---- (b) similarity between independent randomizations -------------------
rand <- lapply(seq_len(2 * as.integer(args$n_pairs)), function(i)
  birewire.rewire.bipartite(m, max.iter = "n", accuracy = args$accuracy, verbose = FALSE))
idx <- seq_len(as.integer(args$n_pairs))
sim_rand_rand <- vapply(idx, function(i)
  birewire.similarity(rand[[2 * i - 1]], rand[[2 * i]]), numeric(1))
sim_orig_rand <- vapply(idx, function(i)
  birewire.similarity(m, rand[[2 * i - 1]]), numeric(1))

sim_dt <- rbind(
  data.table(comparison = "randomized vs randomized", jaccard = sim_rand_rand),
  data.table(comparison = "original vs randomized",   jaccard = sim_orig_rand)
)
fwrite(sim_dt, paste0(args$out_prefix, "_similarity.csv"))

tt <- t.test(sim_rand_rand, sim_orig_rand)
out <- c(
  sprintf("matrix                        : %d genes x %d pathways, %d edges",
          nrow(m), ncol(m), sum(m)),
  sprintf("analytic swap bound N         : %s", ifelse(is.finite(bound), sprintf("%.0f", bound), "n/a")),
  sprintf("Jaccard start / plateau       : %.4f / %.4f", traj[1], plateau),
  sprintf("swaps to plateau              : %s (bound is %s)",
          ifelse(is.na(swaps_to_plateau), "n/a", format(swaps_to_plateau)),
          ifelse(is.finite(bound), sprintf("%.0f", bound), "n/a")),
  sprintf("plateau reached at/before N   : %s",
          ifelse(is.finite(bound) && !is.na(swaps_to_plateau),
                 as.character(swaps_to_plateau <= bound), "n/a")),
  sprintf("Jaccard rand-vs-rand          : %.4f (sd %.4f)", mean(sim_rand_rand), sd(sim_rand_rand)),
  sprintf("Jaccard orig-vs-rand          : %.4f (sd %.4f)", mean(sim_orig_rand), sd(sim_orig_rand)),
  sprintf("difference                    : %+.4f (Welch p = %.3g)",
          mean(sim_rand_rand) - mean(sim_orig_rand), tt$p.value),
  "",
  "Interpretation: a plateau well below 1 reached at or before N, with",
  "rand-vs-rand similarity comparable to orig-vs-rand, means the degree-",
  "preserving randomization space for this simulated architecture is large",
  "and well mixed - the GSR null is not a near-copy of the original database."
)
writeLines(out, paste0(args$out_prefix, "_report.txt"))
cat(paste(out, collapse = "\n"), "\n\n09 birewire_diagnostics: done\n")
