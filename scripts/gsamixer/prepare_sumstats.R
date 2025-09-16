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
             help="Path to GWAS_input.json. If present, column mapping is read for the given trait.") |>
  # Optional explicit column flags override JSON if provided
  add_option("--col-rsid", type="character", default=NULL) |>
  add_option("--col-chr",  type="character", default=NULL) |>
  add_option("--col-pos",  type="character", default=NULL) |>
  add_option("--col-a1",   type="character", default=NULL) |>
  add_option("--col-a2",   type="character", default=NULL) |>
  add_option("--col-n",    type="character", default=NULL) |>
  add_option("--col-z",    type="character", default=NULL) |>
  add_option("--col-beta", type="character", default=NULL) |>
  add_option("--col-se",   type="character", default=NULL) |>
  # NEW: optional OR column
  add_option("--col-or",   type="character", default=NULL) |>
  parse_args(commandArgs(TRUE))

stopifnot(!is.null(opt$trait), !is.null(opt$input), !is.null(opt$output))

# Load GWAS
dt <- fread(opt$input, sep = "\t", header = TRUE, nThread = 1, showProgress = FALSE)

# Helpers
to_key <- function(x) gsub("[^a-z0-9]", "", tolower(x))
cn <- to_key(names(dt))
name_by_key <- setNames(names(dt), cn)

pick_first <- function(keys) {
  hits <- which(cn %in% to_key(keys))
  if (length(hits)) names(dt)[hits[1]] else NA_character_
}

# 1) Start with JSON mapping (if provided)
json_map <- list()
if (!is.null(opt$`traits-config`) && nzchar(opt$`traits-config`)) {
  cfg <- jsonlite::fromJSON(opt$`traits-config`, simplifyVector = TRUE)
  # cfg can be a list of objects or a named list; try to find by "trait" name (case-insensitive)
  rows <- NULL
  if (is.data.frame(cfg)) {
    rows <- cfg
  } else if (is.list(cfg)) {
    # try to bind rows if it’s a list
    try(rows <- data.table::rbindlist(cfg, fill = TRUE), silent = TRUE)
  }
  if (!is.null(rows) && nrow(rows)) {
    trait_row <- rows[tolower(rows$trait) == tolower(opt$trait), ]
    if (nrow(trait_row) == 1) {
      # Accept multiple common field names
      getv <- function(...) {
        for (k in list(...)) {
          if (!is.null(trait_row[[k]]) && nzchar(as.character(trait_row[[k]]))) return(as.character(trait_row[[k]]))
        }
        return(NULL)
      }
      json_map <- list(
        rsid = getv("rsid_col","snp_col","snp","rsid","id"),
        chr  = getv("chr_col","chr","chrom","chromosome"),
        pos  = getv("pos_col","bp_col","pos","bp","position"),
        a1   = getv("a1_col","effect_allele","effectallele","a1"),
        a2   = getv("a2_col","other_allele","otherallele","a2"),
        n    = getv("n_col","n","neff","samplesize"),
        z    = getv("z_col","z"),
        beta = getv("beta_col","beta"),
        se   = getv("se_col","se"),
        # NEW: odds ratio column (if provided in JSON)
        or   = getv("or_col","or","odds_ratio","oddsratio")
      )
    }
  }
}

# 2) Override with explicit CLI flags if present
cli_map <- list(
  rsid = opt$`col-rsid`, chr = opt$`col-chr`, pos = opt$`col-pos`,
  a1   = opt$`col-a1`,   a2  = opt$`col-a2`,  n   = opt$`col-n`,
  z    = opt$`col-z`,    beta= opt$`col-beta`, se  = opt$`col-se`,
  # NEW: CLI override for OR column
  or   = opt$`col-or`
)
for (k in names(cli_map)) {
  if (!is.null(cli_map[[k]]) && nzchar(cli_map[[k]])) json_map[[k]] <- cli_map[[k]]
}

# 3) If still missing, try heuristic detection
heur <- list(
  rsid = pick_first(c("snp","rsid","id")),
  chr  = pick_first(c("chr","chrom","chromosome")),
  pos  = pick_first(c("bp","pos","position")),
  a1   = pick_first(c("a1","effectallele","effect_allele","allele1")),
  a2   = pick_first(c("a2","otherallele","other_allele","allele2")),
  n    = pick_first(c("n","samplesize","neff")),
  z    = pick_first(c("z","zscore","z_stat","zstat","zvalue","z.value")),
  beta = pick_first(c("beta","effect","beta_hat","logor","lnor","log_odds")),
  se   = pick_first(c("se","stderr","se_beta","sebeta")),
  # NEW: detect OR column
  or   = pick_first(c("or","oddsratio","odds_ratio","odds.ratio"))
)
for (k in names(heur)) if (is.null(json_map[[k]]) || !nzchar(json_map[[k]])) json_map[[k]] <- heur[[k]]

# Validate presence of required fields (allow Z or BETA+SE)
missing_basic <- c()
for (k in c("rsid","chr","pos","a1","a2","n")) if (is.null(json_map[[k]]) || !nzchar(json_map[[k]])) missing_basic <- c(missing_basic, toupper(k))
if (length(missing_basic)) {
  stop(sprintf(
    "Missing required columns (from JSON/flags/header): %s. Required: SNP/RSID, CHR, BP/POS, A1/EffectAllele, A2/OtherAllele, N, and Z or BETA+SE or OR+SE.",
    paste(missing_basic, collapse=", ")
  ))
}
if (is.null(json_map$z) || !nzchar(json_map$z)) {
  // Require either (BETA and SE) or (OR and SE)
  if (
    (is.null(json_map$beta) || !nzchar(json_map$beta) || is.null(json_map$se) || !nzchar(json_map$se)) &&
    (is.null(json_map$or)   || !nzchar(json_map$or)   || is.null(json_map$se) || !nzchar(json_map$se))
  ) {
    stop("Provide Z or both BETA and SE, or OR and SE in JSON/flags/headers for GSA-MiXeR.")
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
    // Guard against non-positive OR
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