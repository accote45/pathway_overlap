#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 02 - One-time reference preparation (instructions Sec. 1 + Sec. 4)
#
# Reads the MAGMA gene-location file and the g1000_eur PLINK .bim, and emits:
#   snp_loc.txt      - SNP CHR BP for the ONE-TIME `magma --annotate` run
#   reference.rds    - bim table + gene table + per-gene window SNP counts
#   gene_pool.txt    - genes eligible to be drawn into an architecture
#                      (autosomal, >= min_snps SNPs in their +/-35 kb window)
#   reference_diag.txt
#
# Nothing here depends on a replicate, a condition, or the signal model, so it
# is computed once and shared by every replicate in the sweep.
# ---------------------------------------------------------------------------

parse_args_defaults <- list(
  gene_loc   = "",     # MAGMA gene-location file (params.gene_file)
  bim        = "",     # <bfile>.bim of the LD reference panel
  out_prefix = "ref",
  window_kb  = 35,     # +/- window, matches `magma --annotate window=35,35`
  min_snps   = 10,     # genes with fewer window SNPs are not usable
  max_snps   = 5000    # guard against pathological windows
)
source(file.path(Sys.getenv("SIM_SCRIPTS", unset = "."), "sim_common.R"))
args <- parse_args(parse_args_defaults)

log_header("02 prepare_reference")
if (!nzchar(args$gene_loc) || !nzchar(args$bim)) stop("gene_loc and bim are required")

win <- as.numeric(args$window_kb) * 1000

# ---- gene locations ------------------------------------------------------
# MAGMA gene-loc layout: GENE CHR START STOP [STRAND] [NAME]
genes <- fread(args$gene_loc, header = FALSE, select = 1:4,
               col.names = c("gene_id", "chr", "start", "stop"))
genes[, chr := as.character(chr)]
genes <- genes[chr %in% as.character(1:22)]
genes[, chr := as.integer(chr)]
genes[, `:=`(win_start = pmax(1, start - win), win_stop = stop + win,
             midpoint = (start + stop) / 2)]
cat("Autosomal genes in gene-loc:", nrow(genes), "\n")

# ---- reference panel SNPs ------------------------------------------------
# PLINK .bim: CHR SNP CM BP A1 A2. `row` is the 1-based SNP index in the .bed.
bim <- fread(args$bim, header = FALSE, select = c(1, 2, 4),
             col.names = c("chr", "snp", "bp"))
bim[, row := .I]
bim[, chr := suppressWarnings(as.integer(chr))]
bim <- bim[!is.na(chr) & chr %in% 1:22]
setkey(bim, chr, bp)
cat("Autosomal SNPs in .bim:", nrow(bim), "\n")

# ---- SNPs per gene window ------------------------------------------------
# foverlaps on a keyed range join is far cheaper than a per-gene scan.
gwin <- genes[, .(gene_id, chr, ws = win_start, we = win_stop)]
setkey(gwin, chr, ws, we)
snp_iv <- bim[, .(chr, ws = bp, we = bp)]
setkey(snp_iv, chr, ws, we)
hits <- foverlaps(snp_iv, gwin, type = "within", nomatch = NULL)
counts <- hits[, .N, by = gene_id]
genes <- merge(genes, counts, by = "gene_id", all.x = TRUE)
genes[is.na(N), N := 0L]
setnames(genes, "N", "n_window_snps")

eligible <- genes[n_window_snps >= args$min_snps & n_window_snps <= args$max_snps]
cat("Eligible genes:", nrow(eligible), "\n")
if (nrow(eligible) < 6000) {
  warning("Only ", nrow(eligible), " eligible genes - architecture variation across ",
          "replicates will be limited (need >= ~5000 per replicate).")
}

# ---- outputs -------------------------------------------------------------
fwrite(bim[, .(snp, chr, bp)], paste0(args$out_prefix, "_snp_loc.txt"),
       sep = "\t", col.names = FALSE)
writeLines(eligible$gene_id, paste0(args$out_prefix, "_gene_pool.txt"))
saveRDS(list(bim = bim, genes = genes, window_bp = win),
        paste0(args$out_prefix, "_reference.rds"), compress = FALSE)

qs <- quantile(eligible$n_window_snps, c(0, .25, .5, .75, .95, 1))
diag_lines <- c(
  sprintf("gene_loc               : %s", args$gene_loc),
  sprintf("bim                    : %s", args$bim),
  sprintf("window                 : +/- %g kb", args$window_kb),
  sprintf("autosomal genes        : %d", nrow(genes)),
  sprintf("autosomal SNPs         : %d", nrow(bim)),
  sprintf("eligible genes         : %d (min_snps=%g, max_snps=%g)",
          nrow(eligible), args$min_snps, args$max_snps),
  sprintf("window SNPs per gene   : min %.0f q25 %.0f med %.0f q75 %.0f q95 %.0f max %.0f",
          qs[1], qs[2], qs[3], qs[4], qs[5], qs[6]),
  sprintf("total SNPs in eligible windows (with overlap): %d", sum(eligible$n_window_snps))
)
writeLines(diag_lines, paste0(args$out_prefix, "_reference_diag.txt"))
cat(paste(diag_lines, collapse = "\n"), "\n\n02 prepare_reference: done\n")
