  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)

# --------- Config (edit or add arg parsing if desired) ----------
input_path <- normalizePath("/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/Malacards")
out_dir    <- "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/Malacards"

is_dir <- function(p) isTRUE(file.info(p)$isdir)

# --------- Read MalaCards ----------
read_malacards <- function(f) {
  lines <- readLines(f, warn = FALSE)
  hdr <- which(grepl("\\bSymbol\\b", lines, ignore.case = TRUE) &
                 grepl("\\bScore\\b",  lines, ignore.case = TRUE))[1]
  if (is.na(hdr)) stop("Header with 'Symbol' and 'Score' not found: ", f)
  dt <- suppressWarnings(tryCatch(
    fread(f, skip = hdr - 1, header = TRUE,fill=T),
    error = function(e) as.data.table(read.csv(f, skip = hdr - 1, header = TRUE, check.names = FALSE))
  ))
  setnames(dt, names(dt), make.names(names(dt)))
  sym_col <- names(dt)[tolower(names(dt)) == "symbol"]
  if (!length(sym_col)) sym_col <- names(dt)[grepl("symbol", names(dt), ignore.case = TRUE)]
  if (!length(sym_col)) stop("Symbol column not found after parsing: ", f)
  dt[, Original_Symbol := get(sym_col[1])]
  dt
}

# --------- SYMBOL -> Ensembl mapping only ----------
cat("Loading HGNC SYMBOL -> Ensembl mapping...\n")
sym_map <- as.data.table(AnnotationDbi::select(
  org.Hs.eg.db,
  keys = keys(org.Hs.eg.db, keytype = "SYMBOL"),
  keytype = "SYMBOL",
  columns = c("SYMBOL","ENSEMBL")
))
sym_map <- sym_map[!is.na(ENSEMBL) & ENSEMBL != ""]
setnames(sym_map, "SYMBOL", "Symbol")
sym_map[, Symbol := toupper(Symbol)]

map_symbols <- function(symbols) {
  symbols <- toupper(symbols)
  u <- unique(symbols)
  hits <- sym_map[Symbol %in% u]
  if (!nrow(hits)) {
    return(data.table(Symbol = u, ENSEMBL = NA_character_, Mapping_Source = "UNMAPPED"))
  }
  agg <- hits[, .(
    n_ids = uniqueN(ENSEMBL),
    ENSEMBL_list = paste(unique(ENSEMBL), collapse = ";"),
    First_ENS = unique(ENSEMBL)[1]
  ), by = Symbol]
  agg[, `:=`(
    ENSEMBL = ifelse(n_ids == 1, First_ENS, NA_character_),
    Mapping_Source = fifelse(n_ids == 0, "UNMAPPED",
                       fifelse(n_ids == 1, "SYMBOL", "AMBIGUOUS"))
  )]
  out <- merge(data.table(Symbol = u),
               agg[, .(Symbol, ENSEMBL, Mapping_Source)],
               all.x = TRUE)
  out[is.na(Mapping_Source), `:=`(ENSEMBL = NA_character_, Mapping_Source = "UNMAPPED")]
  out
}

# --------- Collect files ----------
files <- if (is_dir(input_path)) {
  f <- list.files(input_path, pattern = "\\.csv$", full.names = TRUE)
  if (!length(f)) stop("No CSV files in directory: ", input_path)
  f
} else {
  if (!grepl("\\.csv$", input_path, ignore.case = TRUE)) stop("Input must be a .csv file")
  input_path
}
cat("Found", length(files), "MalaCards file(s).\n")

summary_list <- list()
all_map_rows <- list()

for (f in files) {
  cat("Processing:", basename(f), "...\n")
  dt <- tryCatch(read_malacards(f), error = function(e) { cat("  Skipped (error):", e$message, "\n"); return(NULL) })
  if (is.null(dt)) next
  dt[, Symbol_Clean := toupper(gsub("[^A-Za-z0-9_-]", "", Original_Symbol))]
  map_dt <- map_symbols(dt$Symbol_Clean)
  conv <- merge(dt, map_dt, by.x = "Symbol_Clean", by.y = "Symbol", all.x = TRUE)
  conv[, Input_File := basename(f)]

  out_file <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(f)), "_with_ensembl.csv"))
  if (file.exists(out_file) && !overwrite) {
    cat("  Skipping existing:", out_file, "\n")
  } else {
    fwrite(conv, out_file)
    cat("  Wrote:", out_file, "\n")
  }

  summary_list[[length(summary_list) + 1]] <- data.table(
    File = basename(f),
    Total = nrow(conv),
    Unique_Symbols = uniqueN(conv$Symbol_Clean),
    Mapped = sum(!is.na(conv$ENSEMBL)),
    Unmapped = sum(is.na(conv$ENSEMBL)),
    Ambiguous = sum(conv$Mapping_Source == "AMBIGUOUS", na.rm = TRUE)
  )

  all_map_rows[[length(all_map_rows) + 1]] <- unique(conv[, .(Symbol_Clean, ENSEMBL, Mapping_Source)])
}

if (length(summary_list)) {
  summary_dt <- rbindlist(summary_list)
  fwrite(summary_dt, file.path(out_dir, "malacards_mapping_summary.csv"))
  cat("Summary written.\n")
}

if (length(all_map_rows)) {
  master <- unique(rbindlist(all_map_rows))
  fwrite(master, file.path(out_dir, "malacards_symbol_to_ensembl_master.csv"))
  cat("Master mapping written.\n")
}

cat("Done.\n")