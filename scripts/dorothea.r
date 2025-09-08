library(tidyverse)
library(ggplot2)
library(dplyr)
library(dorothea)
library(GSA)
library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)

map_symbols_to_ensembl <- function(symbols) {
  symbols <- as.character(symbols)
  keys <- unique(na.omit(symbols))
  if (length(keys) == 0) return(rep(NA_character_, length(symbols)))

  mapping <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = keys,
    keytype = "SYMBOL",
    columns = c("ENSEMBL")
  )

  # Keep first ENSEMBL per SYMBOL
  mapping <- as.data.frame(mapping, stringsAsFactors = FALSE)
  mapping <- mapping[!is.na(mapping$ENSEMBL) & mapping$ENSEMBL != "", ]
  mapping <- mapping[!duplicated(mapping$SYMBOL), c("SYMBOL", "ENSEMBL")]

  mapping$ENSEMBL[match(symbols, mapping$SYMBOL)]
}

dat <- as.data.frame(entire_database)

# filter for functionally inferred interactions (not literature curated)
dat_nolit <- dat[dat$is_evidence_curated==FALSE,]

# Convert HGNC symbols in tf/target to Ensembl IDs (adds new columns)
dat_nolit <- dat_nolit %>%
  mutate(
    tf = as.character(tf),
    target = as.character(target),
    tf_ensembl = map_symbols_to_ensembl(tf),
    target_ensembl = map_symbols_to_ensembl(target)
  )

# 1. Load pathway data
dat <- GSA.read.gmt("/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/pathway_db/msigdb/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt")
path.list <- dat$genesets
names(path.list) <- dat$geneset.names

# Convert the path list to a data table
genes_long <- rbindlist(lapply(names(path.list), function(name) {
  data.table(value = path.list[[name]], name = name)
}))

# get number of overlapping genes
length(intersect(unique(genes_long$value),c(unique(dat_nolit$tf_ensembl),unique(dat_nolit$target_ensembl))))

# get number of overlapping genes (TF vs msigdb)
length(intersect(unique(genes_long$value),c(unique(dat_nolit$tf_ensembl))))



# Compute the pairwise pathway score for a given gene set
pairwise_pathway_score <- function(genes, dorothea_df,
                                   keep_conf = c("A","B","C","D","E"),
                                   conf_w = c(A=1, B=0.8, C=0.6, D=0.4, E=0.2),
                                   pair_weight_fn = c("min","geom"),    # how to combine two edges for a shared TF
                                   binary = FALSE                       # TRUE = unweighted (0/1) pairs
){
  pair_weight_fn <- match.arg(pair_weight_fn)

  # Coerce to plain character vector (handles list/list-column/factor)
  G <- genes
  G <- as.character(G)
  G <- toupper(G)
  G <- unique(na.omit(G))
  m <- length(G)
  if (m < 2) return(0)

  D <- dorothea_df
  D$tf     <- toupper(D$tf_ensembl)
  D$target <- toupper(D$target_ensembl)
  D$tf     <- sub("\\.\\d+$", "", D$tf)
  D$target <- sub("\\.\\d+$", "", D$target)
  if (length(keep_conf)) D <- D[D$confidence %in% keep_conf, , drop=FALSE]
  if (!nrow(D)) return(0)

  # map confidence to weights
  wmap <- function(x) unname(conf_w[as.character(x)])
  D$w <- wmap(D$confidence); D$w[is.na(D$w)] <- 0

  # --- (1) Direct TF->target evidence INSIDE the pathway ---
  Ein <- D[D$tf %in% G & D$target %in% G, c("tf","target","w"), drop=FALSE]
  # combine multiple edges by taking the max weight
  if (nrow(Ein)) {
    Ein <- aggregate(Ein$w, by=list(tf=Ein$tf, target=Ein$target),
                     FUN=max)
  } else {
    Ein <- data.frame(tf=character(0), target=character(0), x=numeric(0))
  }

  # quick lookup for direct evidence in either direction
  key <- function(a,b) paste0(a,"||",b)
  direct_map <- new.env(parent=emptyenv())
  if (nrow(Ein)) {
    for (k in seq_len(nrow(Ein))) {
      assign(key(Ein$tf[k], Ein$target[k]), Ein$x[k], envir = direct_map)
    }
  }

  # --- (2) Co-target evidence: pairs of targets sharing ANY TF (TF may be outside the pathway) ---
  Etarget <- D[D$target %in% G & !(D$tf %in% G), c("tf","target","w"), drop=FALSE]
  cotarget_map <- new.env(parent=emptyenv())
  if (nrow(Etarget)) {
    # For each TF, get its targets ∩ G and their edge weights
    split_tf <- split(Etarget, Etarget$tf)
    for (t in names(split_tf)) {
      df <- split_tf[[t]]
      if (nrow(df) < 2) next
      # take max weight for duplicate gene pairs
      df <-df %>% group_by(target) %>% slice_max(w, n = 1, with_ties = FALSE) %>% ungroup()
      tg <- df$target; wg <- df$w
      if (length(tg) < 2) next
      # all unordered pairs among tg
      idx <- utils::combn(seq_along(tg), 2)
      for (c in seq_len(ncol(idx))) {
        i <- idx[1, c]; j <- idx[2, c]
        a <- tg[i]; b <- tg[j]
        # combine the two edges for this TF into one unit u_t
        u_t <- switch(pair_weight_fn,
                      min  = min(wg[i], wg[j]),
                      geom = sqrt(wg[i] * wg[j]))
            # max-accumulate across TFs: keep the strongest single-TF evidence
            k <- if (a < b) key(a,b) else key(b,a)
            prev <- if (exists(k, envir=cotarget_map, inherits=FALSE)) get(k, envir=cotarget_map) else 0
            assign(k, max(prev, u_t), envir=cotarget_map)
      }
    }
  }

  # --- Combine per-pair
  genes <- sort(G)
  n_pairs <- choose(m, 2L)
  acc <- 0
  cnt <- 0
  for (u in 1:(m-1)) {
    for (v in (u+1):m) {
      a <- genes[u]; b <- genes[v]
      # direct in either direction:
      d_ab <- if (exists(key(a,b), envir=direct_map, inherits=FALSE)) get(key(a,b), envir=direct_map) else 0
      d_ba <- if (exists(key(b,a), envir=direct_map, inherits=FALSE)) get(key(b,a), envir=direct_map) else 0
      d_pair <- max(d_ab, d_ba)
      # co-target:
      c_pair <- if (exists(key(a,b), envir=cotarget_map, inherits=FALSE)) get(key(a,b), envir=cotarget_map) else 0
        s_ij <- max(d_pair, c_pair)
      if (binary) s_ij <- as.numeric(s_ij > 0)
      acc <- acc + s_ij
      cnt <- cnt + 1L
    }
  }
  acc / cnt
}

# Apply to all pathways
pairwise_score_all <- function(msigdb_list, dorothea_df,
                               keep_conf = c("A","B","C"),
                               conf_w = c(A=1, B=0.8, C=0.6, D=0.4, E=0.2),
                               pair_weight_fn = c("min","geom"),
                               binary = FALSE) {
  pair_weight_fn <- match.arg(pair_weight_fn)
  data.frame(
    pathway = names(msigdb_list),
    size    = vapply(msigdb_list, length, 1L),
    score   = vapply(msigdb_list, function(gs)
      pairwise_pathway_score(gs, dorothea_df, keep_conf, conf_w, pair_weight_fn, binary),
      numeric(1)
    ),
    row.names = NULL
  )[order(-.$score), ]
}

# -------- Run analysis --------

  message("Scoring pathways with Dorothea...")
  scores <- pairwise_score_all(
    msigdb_list = path.list,
    dorothea_df = dat_nolit,
    keep_conf   = c("A","B","C"),
    pair_weight_fn = "geom",
    binary = FALSE
  )

  out_csv <- "dorothea_pairwise_scores.csv"
  write.csv(scores, out_csv, row.names = FALSE)
  message("Wrote: ", out_csv)
  print(head(scores, 10))












