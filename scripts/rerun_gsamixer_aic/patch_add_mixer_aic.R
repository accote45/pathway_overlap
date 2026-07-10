#!/usr/bin/env Rscript
# ============================================================================
# patch_add_mixer_aic.R
#
# Surgically add a `mixer_aic` column to an EXISTING GSA-MiXeR empirical-pvalue
# table, WITHOUT recomputing anything. The AIC value is taken verbatim from the
# raw real GSA-MiXeR result (`loglike_aic` in *_full.go_test_enrich.csv) and
# joined onto the empirical table by pathway id.
#
# Every existing column (std_effect_size, empirical_pval, enrich, ...) is copied
# through untouched. Only `mixer_aic` is appended. This is what unlocks the
# top-500-by-AIC pathway selection in the validation scripts, while leaving the
# empirical / standardized-effect-size corrections exactly as the old run left
# them.
#
# Usage:
#   Rscript patch_add_mixer_aic.R <real_go_test_enrich.csv> <emp_in.txt> <emp_out.txt>
# ============================================================================

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: patch_add_mixer_aic.R <real_go_test_enrich.csv> <emp_in.txt> <emp_out.txt>")
}
real_file <- args[1]
emp_in    <- args[2]
emp_out   <- args[3]

stopifnot(file.exists(real_file), file.exists(emp_in))

# --- Load the raw real GSA-MiXeR result (source of AIC) ---------------------
real <- fread(real_file, header = TRUE)

# pathway id column in the real file (GSA-MiXeR uses "GO")
real_id_candidates <- c("GO", "FULL_NAME", "name", "GENE")
real_id <- real_id_candidates[real_id_candidates %in% names(real)][1]
if (is.na(real_id)) {
  stop("Could not find a pathway-id column in ", real_file,
       " (looked for: ", paste(real_id_candidates, collapse = ", "), ")")
}
if (!"loglike_aic" %in% names(real)) {
  stop("Column 'loglike_aic' not found in ", real_file,
       " -- this file does not carry the native AIC.")
}

# Drop non-pathway rows and de-duplicate on the id (keep first occurrence)
real <- real[!(get(real_id) %in% c("base", "coding_genes", "Base"))]
real <- unique(real, by = real_id)

aic_map <- setNames(as.numeric(real[["loglike_aic"]]), as.character(real[[real_id]]))

# --- Load the existing empirical table (leave every column intact) ----------
emp <- fread(emp_in, header = TRUE)

emp_id_candidates <- c("pathway_name", "FULL_NAME", "GO", "name")
emp_id <- emp_id_candidates[emp_id_candidates %in% names(emp)][1]
if (is.na(emp_id)) {
  stop("Could not find a pathway-id column in ", emp_in,
       " (looked for: ", paste(emp_id_candidates, collapse = ", "), ")")
}

# --- Join AIC on, append as mixer_aic ---------------------------------------
emp[["mixer_aic"]] <- aic_map[as.character(emp[[emp_id]])]

matched  <- sum(!is.na(emp[["mixer_aic"]]))
total    <- nrow(emp)
cat(sprintf("  [patch] %s: matched mixer_aic for %d / %d pathways (id: %s <- %s)\n",
            basename(emp_in), matched, total, emp_id, real_id))
if (matched == 0) {
  stop("No pathways matched between empirical table and real file -- check that ",
       "the pathway ids are on the same naming scheme.")
}
if (matched < total) {
  cat(sprintf("  [patch] WARNING: %d empirical pathways had no AIC (will be dropped ",
              total - matched),
      "by the validation top-500 filter, which requires non-NA mixer_aic).\n", sep = "")
}

fwrite(emp, emp_out, sep = "\t")
cat("  [patch] wrote ", emp_out, "\n", sep = "")
