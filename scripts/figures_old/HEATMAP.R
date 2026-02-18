# Delta-rank correlation heatmap for MAGMA and PRSet
# Rows: traits; Columns: correlation variables (Open Targets, MalaCards, Tissue specificity)
# Only "All Pathways" subset is used

library(data.table)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

# ---- Config ----
trait_map <- list(
  ad     = list(label = "AD",     tissue = "BrainCortex_mean"),
  bmi    = list(label = "BMI",    tissue = "AdiposeSubcutaneous_mean"),
  breast = list(label = "Breast cancer", tissue = "BreastMammaryTissue_mean"),
  cad    = list(label = "CAD",    tissue = "ArteryCoronary_mean"),
  ibd    = list(label = "IBD",    tissue = "ColonTransverse_mean"),
  mdd    = list(label = "MDD",    tissue = "BrainCortex_mean"),
  scz    = list(label = "SCZ",    tissue = "BrainCortex_mean"),
  t2d    = list(label = "T2D",    tissue = "Pancreas_mean")
)

traits <- names(trait_map)
cor_vars <- c("Open Targets", "MalaCards", "DoRothEA", "TissueSpec")
base_dir <- "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results"
delta_rank_dirs <- list(
  "Open Targets" = file.path(base_dir, "delta_rank_ot_correlation"),
  "MalaCards"    = file.path(base_dir, "delta_rank_malacards_correlation"),
  "DoRothEA"      = file.path(base_dir, "delta_rank_dorothea_correlation"),
  "TissueSpec"   = file.path(base_dir, "delta_rank_tissue_correlation")
)

# ---- Helper to read one result ----
get_value <- function(trait, tool_base, cor_type) {
  cor_type_file <- switch(
    cor_type,
    "Open Targets" = "ot",
    "MalaCards"    = "malacards",
    "DoRothEA"      = "dorothea",
    "TissueSpec"   = "tissue"
  )
  file <- file.path(delta_rank_dirs[[cor_type]], trait,
                    sprintf("%s_%s_birewire_delta_rank_%s_correlation_summary.csv",
                            trait, tool_base, cor_type_file))
  if (!file.exists(file)) return(c(NA, NA))
  dt <- fread(file)
  dt <- dt[subset == "All Pathways"]
  if (nrow(dt) != 1) return(c(NA, NA))
  c(dt$spearman_rho_mean, dt$spearman_p_mean)
}

# ---- Build matrices ----
make_matrix <- function(tool_base) {
  rho_mat <- matrix(NA, nrow=length(traits), ncol=length(cor_vars))
  p_mat   <- matrix(NA, nrow=length(traits), ncol=length(cor_vars))
  for (i in seq_along(traits)) {
    for (j in seq_along(cor_vars)) {
      vals <- get_value(traits[i], tool_base, cor_vars[j])
      rho_mat[i, j] <- vals[1]
      p_mat[i, j]   <- vals[2]
    }
  }
  rownames(rho_mat) <- traits
  colnames(rho_mat) <- cor_vars
  rownames(p_mat) <- traits
  colnames(p_mat) <- cor_vars
  list(rho=rho_mat, p=p_mat)
}

magma <- make_matrix("magma")
prset <- make_matrix("prset")

# ---- Custom color scheme ----
my_palette <- colorRampPalette(c("blue", "white", "red"))(99)
breaks <- seq(-0.4, 0.4, length.out=100)

# ---- Mask non-significant correlations ----
mask_non_sig <- function(rho, p, alpha=0.05) {
  out <- rho
  out[p > alpha | is.na(p)] <- NA
  out
}

# ---- Significance annotation function ----
get_sig_asterisks <- function(p_mat) {
  sig_mat <- matrix("", nrow=nrow(p_mat), ncol=ncol(p_mat))
  sig_mat[p_mat < 0.05]    <- "*"
  sig_mat[p_mat < 0.001]   <- "**"
  sig_mat[p_mat < 0.0001]  <- "***"
  sig_mat
}

# ---- Capitalize trait labels ----
trait_labels <- traits
trait_labels <- gsub("^ad$", "AD", trait_labels)
trait_labels <- gsub("^bmi$", "BMI", trait_labels)
trait_labels <- gsub("^breast$", "Breast cancer", trait_labels)
trait_labels <- gsub("^cad$", "CAD", trait_labels)
trait_labels <- gsub("^ibd$", "IBD", trait_labels)
trait_labels <- gsub("^mdd$", "MDD", trait_labels)
trait_labels <- gsub("^scz$", "SCZ", trait_labels)
trait_labels <- gsub("^t2d$", "T2D", trait_labels)

# ---- Make long data frame for ggplot ----
make_long_df <- function(tool_base) {
  df <- data.frame()
  for (i in seq_along(traits)) {
    for (cor_var in cor_vars) {
      vals <- get_value(traits[i], tool_base, cor_var)
      rho <- vals[1]
      p <- vals[2]
      sig <- ""
      if (!is.na(p)) {
        if (p < 0.0001) sig <- "***"
        else if (p < 0.001) sig <- "**"
        else if (p < 0.05) sig <- "*"
      }
      df <- rbind(df, data.frame(Trait=trait_labels[i], Variable=cor_var, Rho=rho, P=p, Sig=sig))
    }
  }
  df
}

magma_df <- make_long_df("magma")
prset_df <- make_long_df("prset")

# ---- Combine data frames and add method column ----
magma_df$Method <- "MAGMA"
prset_df$Method <- "PRSet"
combined_df <- rbind(magma_df, prset_df)

# ---- Order factors for plotting ----
combined_df$Trait <- factor(combined_df$Trait, levels=unique(trait_labels))
combined_df$Variable <- factor(combined_df$Variable, levels=cor_vars)
combined_df$Method <- factor(combined_df$Method, levels=c("MAGMA", "PRSet"))

# ---- Create combined label ----
combined_df$Label <- ifelse(is.na(combined_df$Rho), "", 
                            ifelse(combined_df$Sig == "", sprintf("%.2f", combined_df$Rho), 
                                   paste0(sprintf("%.2f", combined_df$Rho), "\n", combined_df$Sig)))

# ---- Add legend for p-value asterisks ----
asterisk_legend <- data.frame(
  Sig = c("***", "**", "*"),
  Description = c("p < 0.0001", "p < 0.001", "p < 0.05")
)

# ---- Plot combined heatmap ----
heatmap_plot <- ggplot(combined_df, aes(x=Variable, y=Trait, fill=Rho)) +
  geom_tile(color="grey80") +
  geom_text(aes(label=Label), size=4, vjust=0.5, fontface="bold") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, na.value="grey", limits=c(-0.4, 0.4)) +
  facet_grid(Method ~ ., switch="y") +
  labs(title="Delta-Rank Correlations (All Pathways)", x=NULL, y=NULL, fill="Spearman\nrho") +
  theme_minimal(base_size=14) +
  theme(
    strip.text.y.left = element_text(angle=0, face="bold", size=14, color="black"),
    axis.text.x = element_text(angle=45, hjust=1, size=12, color="black"),
    axis.text.y = element_text(size=12, color="black"),
    plot.title = element_text(hjust=0.5, face="bold", size=16, color="black"),
    panel.spacing.y = unit(1, "lines"),
    strip.placement = "outside"
  )

asterisk_table <- tableGrob(
  asterisk_legend,
  rows = NULL,
  theme = ttheme_minimal(
    core = list(fg_params = list(fontface = "bold", fontsize = 12)),
    colhead = list(fg_params = list(fontface = "bold", fontsize = 12))
  )
)

pdf("combined_heatmap_ggplot.pdf", width=8, height=11)
grid.arrange(
  heatmap_plot,
  asterisk_table,
  nrow = 2,
  heights = c(8, 1)
)
dev.off()

# ---- Helper for tissue correlation ----
get_tissue_value <- function(trait, tool_base, tissue_name) {
  trait_dir <- file.path(delta_rank_dirs[["TissueSpec"]], tolower(trait))
  file <- file.path(trait_dir, sprintf("%s_%s_birewire_delta_rank_tissue_correlation_summary.csv",
                                       tolower(trait), tool_base))
  if (!file.exists(file)) return(c(NA, NA))
  dt <- fread(file)
  dt <- dt[subset == "All Pathways" & tissue_metric == tissue_name]
  if (nrow(dt) != 1) return(c(NA, NA))
  c(dt$spearman_rho, dt$spearman_p)
}

# ---- Add tissue results to long data frame ----
add_tissue_results <- function(df, tool_base) {
  for (trait_code in traits) {
    trait_label <- trait_map[[trait_code]]$label
    tissue <- trait_map[[trait_code]]$tissue
    vals <- get_tissue_value(trait_code, tool_base, tissue)
    rho <- vals[1]
    p <- vals[2]
    sig <- ""
    if (!is.na(p)) {
      if (p < 0.0001) sig <- "***"
      else if (p < 0.001) sig <- "**"
      else if (p < 0.05) sig <- "*"
    }
    df <- rbind(df, data.frame(Trait=trait_label, Variable="TissueSpec", Rho=rho, P=p, Sig=sig, stringsAsFactors=FALSE))
  }
  df
}

# ---- Build new long data frames for MAGMA and PRSet ----
magma_df <- make_long_df("magma")
prset_df <- make_long_df("prset")
magma_df <- magma_df[magma_df$Variable != "TissueSpec", ]
prset_df <- prset_df[prset_df$Variable != "TissueSpec", ]
magma_df <- add_tissue_results(magma_df, "magma")
prset_df <- add_tissue_results(prset_df, "prset")

# ---- Combine and plot as before ----
magma_df$Method <- "MAGMA"
prset_df$Method <- "PRSet"
combined_df <- rbind(magma_df, prset_df)

combined_df$Trait <- factor(combined_df$Trait, levels=unique(trait_labels))
combined_df$Variable <- factor(combined_df$Variable, levels=cor_vars)
combined_df$Method <- factor(combined_df$Method, levels=c("MAGMA", "PRSet"))

combined_df$Label <- ifelse(is.na(combined_df$Rho), "", 
                            ifelse(combined_df$Sig == "", sprintf("%.2f", combined_df$Rho), 
                                   paste0(sprintf("%.2f", combined_df$Rho), "\n", combined_df$Sig)))

# ---- Plot as before ----
heatmap_plot <- ggplot(combined_df, aes(x=Variable, y=Trait, fill=Rho)) +
  geom_tile(color="grey80") +
  geom_text(aes(label=Label), size=4, vjust=0.5, fontface="bold") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, na.value="grey", limits=c(-0.4, 0.4)) +
  facet_grid(Method ~ ., switch="y") +
  labs(title="Delta-Rank Correlations (All Pathways)", x=NULL, y=NULL, fill="Spearman\nrho") +
  theme_minimal(base_size=14) +
  theme(
    strip.text.y.left = element_text(angle=0, face="bold", size=14, color="black"),
    axis.text.x = element_text(angle=45, hjust=1, size=12, color="black"),
    axis.text.y = element_text(size=12, color="black"),
    plot.title = element_text(hjust=0.5, face="bold", size=16, color="black"),
    panel.spacing.y = unit(1, "lines"),
    strip.placement = "outside"
  )

asterisk_table <- tableGrob(
  asterisk_legend,
  rows = NULL,
  theme = ttheme_minimal(
    core = list(fg_params = list(fontface = "bold", fontsize = 12)),
    colhead = list(fg_params = list(fontface = "bold", fontsize = 12))
  )
)

pdf("combined_heatmap_ggplot.pdf", width=8, height=11)
grid.arrange(
  heatmap_plot,
  asterisk_table,
  nrow = 2,
  heights = c(8, 1)
)
dev.off()










