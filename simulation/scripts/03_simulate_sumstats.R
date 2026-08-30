#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 03 - Simulate GWAS summary statistics, MVN under real LD (instructions Sec. 4)
#
# Per replicate: take this replicate's per-gene signal assignment, pull the real
# LD for the SNPs in those genes' +/-35 kb windows out of the reference panel,
# and draw  Z ~ MVN(m, R)  block by block, where m is zero except at each causal
# gene's most-central SNP, which carries that gene's tier mean shift mu.
#
# Only SNPs in the selected genes' windows are simulated: MAGMA's annotation
# restricts the gene analysis to window SNPs anyway, so this is the full set of
# SNPs the analysis can use, and the gene background becomes exactly the
# simulated pathway universe.
#
# Outputs:
#   <prefix>_sumstats.txt      - SNP CHR BP Z P  (MAGMA reads SNP + P)
#   <prefix>_causal_snps.tsv   - gene_id, causal snp, mu   (Sec. 8 sanity check)
#   <prefix>_sumstats_diag.txt
# ---------------------------------------------------------------------------

parse_args_defaults <- list(
  reference_rds  = "",
  gene_signal    = "",
  bfile          = "",      # PLINK prefix; .bed/.bim/.fam must exist
  out_prefix     = "sim",
  seed           = 1,
  max_block_snps = 2000,    # cap on one MVN draw; larger merged regions are split
  ld_lambda      = 0.01,    # shrinkage toward identity, keeps R positive definite
  independent    = FALSE    # TRUE = Sec.4 "fast fallback": iid N(0,1), NO LD
)
source(file.path(Sys.getenv("SIM_SCRIPTS", unset = "."), "sim_common.R"))
args <- parse_args(parse_args_defaults)

set.seed(as.integer(args$seed))
log_header("03 simulate_sumstats | seed=", args$seed,
           if (isTRUE(args$independent)) " | INDEPENDENT-Z FALLBACK (not for final results)" else "")

ref <- readRDS(args$reference_rds)
bim <- ref$bim
genes_all <- ref$genes
win <- ref$window_bp

sig <- fread(args$gene_signal)
sel <- merge(genes_all, sig[, .(gene_id, mu, tier, causal)], by = "gene_id")
if (nrow(sel) != nrow(sig)) {
  stop(nrow(sig) - nrow(sel), " selected genes are missing from the reference gene table")
}
cat("Selected genes:", nrow(sel), " causal:", sum(sel$causal), "\n")

fam_lines <- readLines(paste0(args$bfile, ".fam"))
n_samples <- length(fam_lines)
bed_path <- paste0(args$bfile, ".bed")
cat("Reference panel samples:", n_samples, "\n")

# ---- merge selected gene windows into LD blocks --------------------------
blocks <- rbindlist(lapply(sort(unique(sel$chr)), function(cc) {
  w <- sel[chr == cc][order(win_start)]
  bs <- w$win_start[1]; be <- w$win_stop[1]
  out <- list()
  for (i in seq_len(nrow(w))[-1]) {
    if (w$win_start[i] <= be) {
      be <- max(be, w$win_stop[i])
    } else {
      out[[length(out) + 1]] <- data.table(chr = cc, start = bs, end = be)
      bs <- w$win_start[i]; be <- w$win_stop[i]
    }
  }
  out[[length(out) + 1]] <- data.table(chr = cc, start = bs, end = be)
  rbindlist(out)
}))

# Attach the contiguous .bim row run for each merged region, then split any
# region whose SNP count exceeds max_block_snps into consecutive chunks.
setkey(bim, chr, bp)
rows_for <- function(cc, s, e) bim[.(cc), on = .(chr)][bp >= s & bp <= e, row]
maxb <- as.integer(args$max_block_snps)
block_list <- list()
for (i in seq_len(nrow(blocks))) {
  rr <- rows_for(blocks$chr[i], blocks$start[i], blocks$end[i])
  if (length(rr) == 0) next
  n_chunk <- ceiling(length(rr) / maxb)
  chunk_id <- rep(seq_len(n_chunk), each = ceiling(length(rr) / n_chunk))[seq_along(rr)]
  for (k in seq_len(n_chunk)) {
    rk <- rr[chunk_id == k]
    if (length(rk) == 0) next
    block_list[[length(block_list) + 1]] <-
      data.table(block = length(block_list) + 1L, chr = blocks$chr[i],
                 first_row = min(rk), n_snps = length(rk))
  }
}
blocks <- rbindlist(block_list)
blocks[, last_row := first_row + n_snps - 1L]
cat("LD blocks:", nrow(blocks), " total SNPs:", sum(blocks$n_snps),
    " max block:", max(blocks$n_snps), "\n")

# ---- assign each gene to the block holding most of its window SNPs -------
bim_blk <- rbindlist(lapply(seq_len(nrow(blocks)), function(i)
  data.table(row = blocks$first_row[i]:blocks$last_row[i], block = blocks$block[i])))
snp_tab <- merge(bim[row %in% bim_blk$row], bim_blk, by = "row")
setkey(snp_tab, chr, bp)

gene_snps <- rbindlist(lapply(seq_len(nrow(sel)), function(i) {
  s <- snp_tab[.(sel$chr[i]), on = .(chr)][bp >= sel$win_start[i] & bp <= sel$win_stop[i]]
  if (nrow(s) == 0) return(NULL)
  s[, .(gene_id = sel$gene_id[i], snp, row, bp, block,
        dist = abs(bp - sel$midpoint[i]))]
}))
gene_block <- gene_snps[, .N, by = .(gene_id, block)][order(gene_id, -N)][, .SD[1], by = gene_id]

# Causal SNP = the SNP nearest the gene midpoint, restricted to that gene's block.
causal <- merge(gene_snps, gene_block[, .(gene_id, block)], by = c("gene_id", "block"))
causal <- causal[order(gene_id, dist)][, .SD[1], by = gene_id]
causal <- merge(causal, sig[, .(gene_id, mu, tier, is_causal = causal)], by = "gene_id")
n_no_snps <- nrow(sel) - length(unique(gene_snps$gene_id))
if (n_no_snps > 0) warning(n_no_snps, " selected gene(s) had no SNPs in their window")

# ---- simulate block by block ---------------------------------------------
setkey(causal, block)
z_out <- vector("list", nrow(blocks))
lambdas <- numeric(nrow(blocks))
n_dropped <- 0L
n_signal_placed <- 0L

for (i in seq_len(nrow(blocks))) {
  fr <- blocks$first_row[i]; np <- blocks$n_snps[i]
  rows_i <- fr:(fr + np - 1L)
  info <- bim[row %in% rows_i][order(row)]

  if (isTRUE(args$independent)) {
    keep <- seq_len(np); L <- NULL; lambdas[i] <- NA_real_
  } else {
    G <- read_bed_run(bed_path, n_samples, fr, np)
    ld <- ld_from_genotypes(G)
    keep <- ld$keep
    n_dropped <- n_dropped + (np - length(keep))
    if (length(keep) == 0) next
    ch <- shrunk_chol(ld$R, lambda = args$ld_lambda)
    L <- ch$L; lambdas[i] <- ch$lambda
  }
  info_k <- info[keep]
  p <- nrow(info_k)

  # mean vector: mu at each causal gene's causal SNP (summed on the rare
  # occasion two genes resolve to the same SNP)
  m <- numeric(p)
  cs <- causal[.(blocks$block[i]), nomatch = NULL][is_causal == TRUE & mu > 0]
  if (nrow(cs) > 0) {
    pos <- match(cs$row, info_k$row)
    ok <- !is.na(pos)
    if (any(ok)) {
      m[pos[ok]] <- m[pos[ok]] + cs$mu[ok]
      n_signal_placed <- n_signal_placed + sum(ok)
    }
    if (any(!ok)) {
      # causal SNP was monomorphic and dropped: fall back to the nearest kept SNP
      for (j in which(!ok)) {
        cand <- which.min(abs(info_k$bp - cs$bp[j]))
        if (length(cand) == 1) {
          m[cand] <- m[cand] + cs$mu[j]
          n_signal_placed <- n_signal_placed + 1L
        }
      }
    }
  }

  z <- if (isTRUE(args$independent)) m + rnorm(p) else as.vector(m + L %*% rnorm(p))
  z_out[[i]] <- data.table(snp = info_k$snp, chr = info_k$chr, bp = info_k$bp, Z = z)
}

sumstats <- rbindlist(z_out)
sumstats[, P := 2 * pnorm(-abs(Z))]
sumstats[P < 1e-300, P := 1e-300]
setorder(sumstats, chr, bp)

fwrite(sumstats[, .(SNP = snp, CHR = chr, BP = bp, Z = round(Z, 6), P = P)],
       paste0(args$out_prefix, "_sumstats.txt"), sep = "\t")
fwrite(causal[, .(gene_id, tier, mu, is_causal, causal_snp = snp, snp_bp = bp, block)],
       paste0(args$out_prefix, "_causal_snps.tsv"), sep = "\t")

diag_lines <- c(
  sprintf("seed                   : %s", args$seed),
  sprintf("mode                   : %s", if (isTRUE(args$independent)) "INDEPENDENT (no LD)" else "MVN under panel LD"),
  sprintf("selected genes         : %d (causal %d)", nrow(sel), sum(sel$causal)),
  sprintf("LD blocks              : %d", nrow(blocks)),
  sprintf("block SNPs min/med/max : %d / %.0f / %d",
          min(blocks$n_snps), median(blocks$n_snps), max(blocks$n_snps)),
  sprintf("SNPs simulated         : %d", nrow(sumstats)),
  sprintf("monomorphic dropped    : %d", n_dropped),
  sprintf("causal SNPs placed     : %d of %d causal genes", n_signal_placed, sum(sel$causal)),
  sprintf("shrinkage lambda (max) : %s",
          if (all(is.na(lambdas))) "n/a" else sprintf("%.4f", max(lambdas, na.rm = TRUE))),
  sprintf("genome-wide sig SNPs   : %d (P < 5e-8)", sum(sumstats$P < 5e-8)),
  sprintf("lambda_GC              : %.4f",
          median(qchisq(sumstats$P, df = 1, lower.tail = FALSE)) / qchisq(0.5, 1))
)
writeLines(diag_lines, paste0(args$out_prefix, "_sumstats_diag.txt"))
cat(paste(diag_lines, collapse = "\n"), "\n\n03 simulate_sumstats: done\n")
