library(tidyverse)
library(ggplot2)
library(dplyr)
library(dorothea)
library(GSA)
library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(parallel)

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
  
  # Create a gene index map for fast lookups
  gene_to_idx <- setNames(seq_along(G), G)
  
  # Create matrices to store direct evidence and co-target evidence
  direct_evidence <- matrix(0, m, m)
  cotarget_evidence <- matrix(0, m, m)
  
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
  
  # Fill the direct evidence matrix
  if (nrow(Ein) > 0) {
    # Aggregate to get max weight for each TF-target pair
    Ein <- as.data.table(Ein)[, .(x = max(w)), by = .(tf, target)]
    
    # Map to matrix indices and fill the matrix
    for (i in 1:nrow(Ein)) {
      tf_idx <- gene_to_idx[Ein$tf[i]]
      target_idx <- gene_to_idx[Ein$target[i]]
      direct_evidence[tf_idx, target_idx] <- Ein$x[i]
    }
  }
  
  # --- (2) Co-target evidence ---
  Etarget <- D[D$target %in% G & !(D$tf %in% G), c("tf","target","w"), drop=FALSE]
  
  if (nrow(Etarget) > 0) {
    # Process by TF groups more efficiently
    tf_groups <- split(Etarget, Etarget$tf)
    
    for (tf_data in tf_groups) {
      if (nrow(tf_data) < 2) next
      
      # Take max weight per target
      tf_data <- tf_data %>% 
        group_by(target) %>% 
        slice_max(w, n = 1, with_ties = FALSE) %>% 
        ungroup()
      
      targets <- tf_data$target
      weights <- tf_data$w
      
      target_indices <- gene_to_idx[targets]
      n_targets <- length(target_indices)
      
      if (n_targets < 2) next
      
      # Process all pairs for this TF
      for (i in 1:(n_targets-1)) {
        for (j in (i+1):n_targets) {
          idx_i <- target_indices[i]
          idx_j <- target_indices[j]
          
          # Combine weights based on method
          combined_weight <- switch(pair_weight_fn,
                                  min = min(weights[i], weights[j]),
                                  geom = sqrt(weights[i] * weights[j]))
          
          # Update co-target evidence (symmetric matrix)
          cotarget_evidence[idx_i, idx_j] <- max(cotarget_evidence[idx_i, idx_j], combined_weight)
          cotarget_evidence[idx_j, idx_i] <- cotarget_evidence[idx_i, idx_j]
        }
      }
    }
  }
  
  # --- Combine evidence and calculate score ---
  total_score <- 0
  pair_count <- 0
  
  for (i in 1:(m-1)) {
    for (j in (i+1):m) {
      # Get maximum of direct evidence in either direction
      d_pair <- max(direct_evidence[i,j], direct_evidence[j,i])
      
      # Get co-target evidence
      c_pair <- cotarget_evidence[i,j]
      
      # Overall evidence for the pair
      s_ij <- max(d_pair, c_pair)
      if (binary) s_ij <- as.numeric(s_ij > 0)
      
      total_score <- total_score + s_ij
      pair_count <- pair_count + 1
    }
  }
  
  # Return average score
  return(total_score / pair_count)
}

# Apply to all pathways
pairwise_score_all <- function(msigdb_list, dorothea_df,
                               keep_conf = c("A","B","C"),
                               conf_w = c(A=1, B=0.8, C=0.6, D=0.4, E=0.2),
                               pair_weight_fn = c("min","geom"),
                               binary = FALSE) {
  pair_weight_fn <- match.arg(pair_weight_fn)
  n_pathways <- length(msigdb_list)
  results <- data.frame(
    pathway = names(msigdb_list),
    size = integer(n_pathways),
    score = numeric(n_pathways)
  )
  
  cores <- max(1, detectCores() - 1)  # Leave one core free
  results <- mcmapply(function(gs, name) {
    c(size = length(gs),
      score = pairwise_pathway_score(gs, dorothea_df, keep_conf, conf_w, pair_weight_fn, binary))
  }, msigdb_list, names(msigdb_list), mc.cores = cores, SIMPLIFY = FALSE)

  # Convert to data frame
  results_df <- as.data.frame(do.call(rbind, results))
  results_df$pathway <- names(msigdb_list)
  results_df <- results_df[order(-results_df$score), ]
  
  return(results_df)
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












