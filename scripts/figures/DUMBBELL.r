library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalt)
library(patchwork)

options(dplyr.summarise.inform = FALSE)

# ------------------------------------------------------------------------------
# Config
# ------------------------------------------------------------------------------
SUBSETS <- c("Top 500 Pathways")
METHODS3 <- c("BireWire_EmpPvalStdBeta", "PvalueBeta")

ANALYSIS_FOLDERS <- list(
  Malacards = "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/malacards_correlation",
  DoRothEA = "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/dorothea_correlation",
  OpenTargets = "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/opentargets_correlation",
  TissueSpecificity = "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/tissue_correlation"
)

panel_list_magma <- list()
panel_list_prset <- list()
trait_order <- NULL

parse_and_normalize <- function(f) {
  bn <- basename(f)
  mm <- regmatches(bn, regexec("^(.+?)_([A-Za-z0-9]+)_rank_correlation_summary\\.csv$", bn))[[1]]
  if (length(mm) != 3) return(NULL)
  d <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  d[, trait := mm[2]]
  d[, tool_base := mm[3]]
  d
}

for (analysis_name in names(ANALYSIS_FOLDERS)) {
  files <- list.files(ANALYSIS_FOLDERS[[analysis_name]], pattern = "_rank_correlation_summary\\.csv$", recursive = TRUE, full.names = TRUE)
  dat <- rbindlist(Filter(Negate(is.null), lapply(files, parse_and_normalize)), fill = TRUE)
  if (nrow(dat) == 0) next

  dat <- dat %>%
    mutate(
      tool_base = case_when(
        tool_base == "magma_birewire" ~ "magma",
        tool_base == "birewire" & grepl("_prset$", trait) ~ "prset",
        TRUE ~ tool_base
      ),
      trait = if_else(tool_base == "prset", sub("_prset$", "", trait), trait)
    ) %>%
    filter(tool_base %in% c("magma", "prset"), subset %in% SUBSETS, method %in% METHODS3)

  if (nrow(dat) == 0) next

  plot_dat <- dat %>%
    mutate(
      trait_label = case_when(
        grepl("IBD|inflamm", trait, ignore.case = TRUE) ~ "IBD",
        grepl("\\bCAD\\b|coronary", trait, ignore.case = TRUE) ~ "CAD",
        grepl("\\bAD\\b|alzheimer", trait, ignore.case = TRUE) ~ "AD",
        grepl("\\bT2D\\b|type\\s*2.*diab|t2d", trait, ignore.case = TRUE) ~ "T2D",
        grepl("SCZ|schizo", trait, ignore.case = TRUE) ~ "SCZ",
        grepl("MDD|depress", trait, ignore.case = TRUE) ~ "MDD",
        grepl("breast", trait, ignore.case = TRUE) ~ "Breast cancer",
        TRUE ~ trait
      )
    ) %>%
    select(tool_base, trait_label, method, rho = rank_mean_score_correlation, pval = rank_mean_score_pvalue) %>%
    pivot_wider(names_from = method, values_from = c(rho, pval))

  # Set trait order from first analysis
  if (is.null(trait_order) && any(plot_dat$tool_base == "magma")) {
    trait_order <- plot_dat %>% filter(tool_base == "magma") %>% pull(trait_label)
  }

  annotate_significance <- function(pval) {
    case_when(
      pval < 0.0001 ~ "***",
      pval < 0.001 ~ "**",
      pval < 0.05 ~ "*",
      TRUE ~ ""
    )
  }

  plot_dat <- plot_dat %>%
    mutate(
      star_BireWire = annotate_significance(pval_BireWire_EmpPvalStdBeta),
      star_Pvalue = annotate_significance(pval_PvalueBeta)
    )

  make_panel <- function(df, show_y) {
    ggplot(df, aes(
      x = rho_BireWire_EmpPvalStdBeta,
      xend = rho_PvalueBeta,
      y = factor(trait_label, levels = trait_order)
    )) +
      geom_vline(xintercept = 0, color = "black", linewidth = 1.2) +
      geom_dumbbell(colour = "#a3c4dc", size = 3, colour_x = "#0e668b", colour_xend = "#dd1c77") +
      geom_text(aes(x = rho_BireWire_EmpPvalStdBeta, label = star_BireWire), vjust = -1, size = 5, color = "#0e668b") +
      geom_text(aes(x = rho_PvalueBeta, label = star_Pvalue), vjust = -1, size = 5, color = "#dd1c77") +
      labs(
        y = if (show_y) "Trait" else NULL,
        x = "Spearman rho",
        title = analysis_name
      ) +
      xlim(min(df$rho_BireWire_EmpPvalStdBeta, df$rho_PvalueBeta, na.rm = TRUE),
           max(df$rho_BireWire_EmpPvalStdBeta, df$rho_PvalueBeta, na.rm = TRUE)) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.y = if (show_y) element_text(size = 10, color = "black") else element_blank(),
        axis.title.y = if (show_y) element_text(color = "black") else element_blank(),
        axis.title.x = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5)
      )
  }

  magma_df <- plot_dat %>% filter(tool_base == "magma")
  prset_df <- plot_dat %>% filter(tool_base == "prset")

  if (nrow(magma_df) > 0) panel_list_magma[[analysis_name]] <- make_panel(magma_df, show_y = analysis_name == names(ANALYSIS_FOLDERS)[1])
  if (nrow(prset_df) > 0) panel_list_prset[[analysis_name]] <- make_panel(prset_df, show_y = FALSE)
}

# Remove empty panels
panel_list_magma <- panel_list_magma[!vapply(panel_list_magma, is.null, logical(1))]
panel_list_prset <- panel_list_prset[!vapply(panel_list_prset, is.null, logical(1))]

# Combine and save
p_combined <- wrap_plots(
  wrap_plots(panel_list_magma, ncol = 1),
  wrap_plots(panel_list_prset, ncol = 1),
  ncol = 2
)

ggsave(
  filename = "dumbbell_birewire_rawp_top500_panel_multiscore.pdf",
  plot = p_combined, width = 12, height = 4.2 * length(panel_list_magma), device = cairo_pdf
)
