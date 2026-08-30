# Shared helpers for the GSR simulation sub-pipeline.
# Sourced by the 0*_ scripts; kept dependency-light (base R + data.table only).

suppressPackageStartupMessages(library(data.table))

# Parse `key=value` command-line arguments into a named list, with defaults.
# Values are coerced to the type of the matching default (numeric/logical/character).
parse_args <- function(defaults = list()) {
  argv <- commandArgs(trailingOnly = TRUE)
  kv <- strsplit(argv[grepl("=", argv, fixed = TRUE)], "=", fixed = TRUE)
  supplied <- setNames(
    vapply(kv, function(x) paste(x[-1], collapse = "="), character(1)),
    vapply(kv, `[`, character(1), 1)
  )

  unknown <- setdiff(names(supplied), names(defaults))
  if (length(unknown) > 0) {
    stop("Unknown argument(s): ", paste(unknown, collapse = ", "),
         "\nKnown: ", paste(names(defaults), collapse = ", "))
  }

  out <- defaults
  for (k in names(supplied)) {
    proto <- defaults[[k]]
    raw <- supplied[[k]]
    out[[k]] <- if (is.logical(proto)) {
      tolower(raw) %in% c("true", "t", "yes", "1")
    } else if (is.numeric(proto)) {
      as.numeric(raw)
    } else {
      raw
    }
  }
  out
}

# Write a pathway list (named list of character vectors) as a GMT file in the
# exact layout the main pipeline uses: name \t PLACEHOLDER \t gene1 \t gene2 ...
write_gmt <- function(path_list, file) {
  con <- file(file, "wt")
  on.exit(close(con))
  for (nm in names(path_list)) {
    writeLines(paste(c(nm, "PLACEHOLDER", path_list[[nm]]), collapse = "\t"), con)
  }
  invisible(file)
}

# Read a GMT into a named list of character vectors (no GSA dependency).
read_gmt <- function(file) {
  lines <- readLines(file, warn = FALSE)
  lines <- lines[nzchar(lines)]
  parts <- strsplit(lines, "\t", fixed = TRUE)
  out <- lapply(parts, function(x) x[-c(1, 2)])
  names(out) <- vapply(parts, `[`, character(1), 1)
  out
}

# Long (gene, pathway) edge table from a pathway list.
gmt_to_edges <- function(path_list) {
  rbindlist(lapply(names(path_list), function(nm) {
    data.table(pathway = nm, gene = path_list[[nm]])
  }))
}

log_header <- function(...) {
  cat("========================================\n")
  cat(..., "\n", sep = "")
  cat("========================================\n")
}

# ---------------------------------------------------------------------------
# Minimal PLINK 1 .bed reader.
#
# Avoids a plink/snpStats dependency on the cluster. The .bed layout is
# SNP-major: 3 magic bytes, then ceil(n/4) bytes per SNP, 2 bits per sample
# read from the low end of each byte. Codes: 00 = hom A1, 01 = missing,
# 10 = het, 11 = hom A2 -> A1 dosages 2, NA, 1, 0.
#
# `snp_rows` must be a CONTIGUOUS increasing run of 1-based .bim row indices,
# which is what genomic blocks give us on a position-sorted .bim; the whole
# run is then pulled with a single seek + read.
# ---------------------------------------------------------------------------
.bed_lut <- local({
  code_to_dosage <- c(2L, NA_integer_, 1L, 0L)   # codes 00,01,10,11
  lut <- matrix(NA_integer_, nrow = 4, ncol = 256)
  for (b in 0:255) {
    for (j in 0:3) lut[j + 1, b + 1] <- code_to_dosage[bitwAnd(bitwShiftR(b, 2 * j), 3L) + 1L]
  }
  lut
})

read_bed_run <- function(bed_path, n_samples, first_row, n_snps) {
  bytes_per_snp <- ceiling(n_samples / 4)
  con <- file(bed_path, "rb")
  on.exit(close(con))
  magic <- readBin(con, "raw", 3)
  if (length(magic) < 3 || magic[1] != as.raw(0x6c) || magic[2] != as.raw(0x1b)) {
    stop("Not a PLINK .bed file: ", bed_path)
  }
  if (magic[3] != as.raw(0x01)) stop("Only SNP-major .bed is supported: ", bed_path)

  seek(con, where = 3 + as.numeric(first_row - 1) * bytes_per_snp, origin = "start")
  raws <- readBin(con, "raw", n = bytes_per_snp * n_snps)
  if (length(raws) < bytes_per_snp * n_snps) stop("Short read from ", bed_path)

  g <- .bed_lut[, as.integer(raws) + 1L]        # 4 x (bytes_per_snp * n_snps)
  dim(g) <- c(4L * bytes_per_snp, n_snps)
  g[seq_len(n_samples), , drop = FALSE]
}

# Column-standardised genotype matrix -> LD correlation matrix.
# Missing calls are mean-imputed; monomorphic columns are reported so the
# caller can drop those SNPs entirely (MAGMA would exclude them anyway).
ld_from_genotypes <- function(G) {
  n <- nrow(G)
  mu <- colMeans(G, na.rm = TRUE)
  na_idx <- which(is.na(G), arr.ind = TRUE)
  if (nrow(na_idx) > 0) G[na_idx] <- mu[na_idx[, 2]]
  Gc <- sweep(G, 2, colMeans(G), "-")
  sds <- sqrt(colSums(Gc^2) / (n - 1))
  keep <- which(sds > 1e-8 & is.finite(sds))
  if (length(keep) == 0) return(list(R = NULL, keep = keep))
  Gs <- sweep(Gc[, keep, drop = FALSE], 2, sds[keep], "/")
  R <- crossprod(Gs) / (n - 1)                  # BLAS-backed
  list(R = R, keep = keep)
}

# Cholesky factor of a shrunk LD matrix. The panel has far fewer samples than
# a large block has SNPs, so R is rank-deficient; shrink toward the identity
# (which also keeps unit variances) and escalate lambda if chol still fails.
shrunk_chol <- function(R, lambda = 0.01, max_lambda = 0.5) {
  p <- nrow(R)
  repeat {
    Rl <- (1 - lambda) * R
    diag(Rl) <- 1
    L <- tryCatch(t(chol(Rl)), error = function(e) NULL)
    if (!is.null(L)) return(list(L = L, lambda = lambda))
    lambda <- lambda * 2
    if (lambda > max_lambda) stop("Cholesky failed up to lambda=", max_lambda, " for p=", p)
  }
}
