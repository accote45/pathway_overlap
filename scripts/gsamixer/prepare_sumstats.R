#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(jsonlite)
})

opt <- OptionParser() |>
  add_option("--trait", type="character") |>
  add_option("--input", type="character") |>
  add_option("--output", type="character") |>
  add_option("--traits-config", type="character", default=NULL,
             help="Path to GWAS_input.json. Column mapping must be specified for the given trait.") |>
  parse_args(commandArgs(TRUE))

stopifnot(!is.null(opt$trait), !is.null(opt$input), !is.null(opt$output))

# Verify traits-config is provided
if (is.null(opt$`traits-config`) || !nzchar(opt$`traits-config`)) {
  stop("--traits-config parameter is required. Please provide a JSON file with column mappings.")
}

# Load GWAS
dt <- fread(opt$input, sep = "\t", header = TRUE, nThread = 1, showProgress = FALSE)

# Helpers
to_key <- function(x) gsub("[^a-z0-9]", "", tolower(x))
cn <- to_key(names(dt))
name_by_key <- setNames(names(dt), cn)

# Load JSON mapping (required)
json_map <- list()
cfg <- jsonlite::fromJSON(opt$`traits-config`, simplifyVector = TRUE)
# cfg can be a list of objects or a named list; try to find by "trait" name (case-insensitive)
rows <- NULL
if (is.data.frame(cfg)) {
  rows <- cfg
} else if (is.list(cfg)) {
  # try to bind rows if it's a list
  try(rows <- data.table::rbindlist(cfg, fill = TRUE), silent = TRUE)
}

if (is.null(rows) || !nrow(rows)) {
  stop("Invalid JSON configuration format. Expected a data frame or list of trait configurations.")
}

trait_row <- rows[tolower(rows$trait) == tolower(opt$trait), ]
if (nrow(trait_row) == 0) {
  stop(sprintf("Trait '%s' not found in JSON configuration. Available traits: %s", 
               opt$trait, paste(rows$trait, collapse=", ")))
}

# Direct mapping from JSON structure (no fallbacks needed)
json_map <- list(
  rsid = as.character(trait_row$rsid_col),
  chr  = as.character(trait_row$chr_col),
  pos  = as.character(trait_row$pos_col),
  a1   = as.character(trait_row$effect_allele),
  a2   = as.character(trait_row$other_allele),
  n    = as.character(trait_row$n_col),
  beta = as.character(trait_row$summary_statistic_name)
  # Add other fields as needed
)

# Validate presence of required fields (allow Z or BETA+SE)
missing_basic <- c()
for (k in c("rsid","chr","pos","a1","a2","n")) if (is.null(json_map[[k]]) || !nzchar(json_map[[k]])) missing_basic <- c(missing_basic, toupper(k))
if (length(missing_basic)) {
  stop(sprintf(
    "Missing required columns in JSON configuration for trait '%s': %s. Required: SNP/RSID, CHR, BP/POS, A1/EffectAllele, A2/OtherAllele, N, and Z or BETA+SE or OR+SE.",
    opt$trait, paste(missing_basic, collapse=", ")
  ))
}
if (is.null(json_map$z) || !nzchar(json_map$z)) {
  # Require either (BETA and SE) or (OR and SE)
  if (
    (is.null(json_map$beta) || !nzchar(json_map$beta) || is.null(json_map$se) || !nzchar(json_map$se)) &&
    (is.null(json_map$or)   || !nzchar(json_map$or)   || is.null(json_map$se) || !nzchar(json_map$se))
  ) {
    stop("Provide Z or both BETA and SE, or OR and SE in JSON configuration for trait '", opt$trait, "'.")
  }
}

# Extract columns safely by display name (case-insensitive)
col_get <- function(colname) if (is.null(colname) || !nzchar(colname)) NULL else dt[[ name_by_key[[ to_key(colname) ]] ]]

RSID <- as.character(col_get(json_map$rsid))
CHR  <- suppressWarnings(as.integer(col_get(json_map$chr)))
POS  <- suppressWarnings(as.integer(col_get(json_map$pos)))
A1   <- toupper(as.character(col_get(json_map$a1)))
A2   <- toupper(as.character(col_get(json_map$a2)))
N    <- suppressWarnings(as.numeric(col_get(json_map$n)))

Z <- NULL
if (!is.null(json_map$z) && nzchar(json_map$z)) {
  Z <- suppressWarnings(as.numeric(col_get(json_map$z)))
} else {
  BETA <- if (!is.null(json_map$beta) && nzchar(json_map$beta)) suppressWarnings(as.numeric(col_get(json_map$beta))) else NULL
  SE   <- if (!is.null(json_map$se)   && nzchar(json_map$se))   suppressWarnings(as.numeric(col_get(json_map$se)))   else NULL
  OR   <- if (!is.null(json_map$or)   && nzchar(json_map$or))   suppressWarnings(as.numeric(col_get(json_map$or)))   else NULL

  if (!is.null(BETA) && !is.null(SE)) {
    Z <- BETA / SE
  } else if (!is.null(OR) && !is.null(SE)) {
    # Guard against non-positive OR
    OR[!is.finite(OR) | OR <= 0] <- NA_real_
    Z <- log(OR) / SE
  } else {
    stop("Could not compute Z: need (BETA and SE) or (OR and SE).")
  }
}

out <- data.table(RSID=RSID, CHR=CHR, POS=POS, A1=A1, A2=A2, N=N, Z=Z)

# Basic QC: drop missing, non-ACGT, identical alleles
out <- out[complete.cases(out)]
out <- out[A1 %in% c("A","C","G","T") & A2 %in% c("A","C","G","T")]
out <- out[A1 != A2]

# Write as TSV.gz (MiXeR accepts .gz)
fwrite(out, opt$output, sep="\t")
cat(sprintf("Trait %s: wrote %s with %d rows\n", opt$trait, opt$output, nrow(out)))