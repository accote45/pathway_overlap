#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 05 - Reverse of 04: split a batched MAGMA .gsa.out back into one .gsa.out per
# random database, with the ORIGINAL set names restored, so that the main
# pipeline's scripts/core/calc_empirical.r consumes them completely unchanged.
#
# Also guarantees a FULL_NAME column. MAGMA only emits FULL_NAME when set names
# are long enough to be truncated in VARIABLE; the simulated names are short, so
# without this the unmodified calc_empirical.r (which reads FULL_NAME for MAGMA)
# would not find its pathway column.
# ---------------------------------------------------------------------------
parse_args_defaults <- list(
  gsa_out   = "",
  out_dir   = ".",
  tag       = "",      # emitted as <tag>_set_random<k>.<method>.gsa.out
  method    = "birewire",
  strip_prefix = TRUE  # FALSE = a non-batched (real) file, just normalise columns
)
source(file.path(Sys.getenv("SIM_SCRIPTS", unset = "."), "sim_common.R"))
args <- parse_args(parse_args_defaults)

raw <- readLines(args$gsa_out, warn = FALSE)
hdr_i <- which(!startsWith(raw, "#") & nzchar(trimws(raw)))[1]
if (is.na(hdr_i)) stop("No header row found in ", args$gsa_out)
comments <- raw[seq_len(hdr_i - 1)]
tab <- fread(text = paste(raw[hdr_i:length(raw)], collapse = "\n"), header = TRUE)

if (!"VARIABLE" %in% names(tab)) stop("No VARIABLE column in ", args$gsa_out)
# FULL_NAME is authoritative when present (VARIABLE may be truncated by MAGMA)
if (!"FULL_NAME" %in% names(tab)) tab[, FULL_NAME := VARIABLE]

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

write_one <- function(dt, path) {
  # comment block + header first, then append the rows (fwrite needs a path,
  # not a connection)
  writeLines(c(comments, paste(names(dt), collapse = "\t")), path)
  fwrite(dt, path, sep = "\t", col.names = FALSE, quote = FALSE, append = TRUE)
}

if (!isTRUE(args$strip_prefix)) {
  write_one(tab, file.path(args$out_dir, basename(args$gsa_out)))
  cat("Normalised", nrow(tab), "sets ->", basename(args$gsa_out), "\n")
} else {
  m <- regmatches(tab$FULL_NAME, regexec("^R([0-9]+)__(.*)$", tab$FULL_NAME))
  bad <- vapply(m, length, integer(1)) != 3
  if (any(bad)) {
    stop(sum(bad), " set name(s) lack the R<k>__ batch prefix, e.g. '",
         tab$FULL_NAME[which(bad)[1]], "' - 04/05 are out of sync")
  }
  tab[, rand_iter := as.integer(vapply(m, `[`, character(1), 2))]
  tab[, FULL_NAME := vapply(m, `[`, character(1), 3)]
  tab[, VARIABLE := FULL_NAME]
  for (k in sort(unique(tab$rand_iter))) {
    sub <- tab[rand_iter == k][, rand_iter := NULL]
    write_one(sub, file.path(args$out_dir,
      sprintf("%s_set_random%d.%s.gsa.out", args$tag, k, args$method)))
  }
  cat("Split", nrow(tab), "rows into", length(unique(tab$rand_iter)),
      "per-database .gsa.out files in", args$out_dir, "\n")
}
