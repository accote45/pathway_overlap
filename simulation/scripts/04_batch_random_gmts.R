#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 04 - Pack many randomized GMTs into ONE MAGMA --set-annot file.
#
# Why: the main pipeline runs one MAGMA job per null database. At
# conditions x replicates x 1000 nulls that is ~10^5-10^6 LSF jobs. MAGMA scores
# each set in --set-annot independently (conditional only on the gene covariates,
# not on the other sets), so concatenating K null databases with prefixed set
# names into one file gives bit-identical per-set results in 1/K the jobs.
# `verify_batching` in the pipeline asserts that equivalence on real output.
#
# Set names become  R<k>__<original>  ; 05_split_batch_gsa.R reverses this.
# ---------------------------------------------------------------------------
parse_args_defaults <- list(
  gmt_dir     = "",    # directory of GeneSet.random<k>.gmt files
  first       = 1,     # first random-set index in this batch (inclusive)
  last        = 100,   # last random-set index in this batch (inclusive)
  out_file    = "batch.setannot",
  index_file  = "batch_index.tsv"
)
source(file.path(Sys.getenv("SIM_SCRIPTS", unset = "."), "sim_common.R"))
args <- parse_args(parse_args_defaults)

ks <- seq(as.integer(args$first), as.integer(args$last))
con <- file(args$out_file, "wt"); on.exit(close(con))
n_sets <- 0L; missing <- integer(0)
for (k in ks) {
  f <- file.path(args$gmt_dir, sprintf("GeneSet.random%d.gmt", k))
  if (!file.exists(f)) { missing <- c(missing, k); next }
  lines <- readLines(f, warn = FALSE)
  lines <- lines[nzchar(lines)]
  # prefix only the set-name field, leave PLACEHOLDER + genes untouched
  writeLines(sub("^([^\t]+)", sprintf("R%d__\\1", k), lines), con)
  n_sets <- n_sets + length(lines)
}
if (length(missing) > 0) {
  stop("Missing randomized GMT(s) for index: ", paste(head(missing, 10), collapse = ", "),
       if (length(missing) > 10) sprintf(" ... (%d total)", length(missing)) else "")
}
fwrite(data.table(first = min(ks), last = max(ks), n_random_sets = length(ks), n_sets = n_sets),
       args$index_file, sep = "\t")
cat("Packed", length(ks), "random databases /", n_sets, "gene sets ->", args$out_file, "\n")
