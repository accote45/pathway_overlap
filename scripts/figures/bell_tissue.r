library(data.table)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------------------------
# Config
# ------------------------------------------------------------------------------
SUBSETS <- c("Top 500 Pathways")
METHODS3 <- c("BireWire_EmpPvalStdBeta", "PvalueBeta")
MALACARDS_FOLDER <- "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/tissue_correlation"

trait_map <- list(
  ad     = list(label = "AD",     tissue = "BrainCortex"),
  bmi    = list(label = "BMI",    tissue = "AdiposeSubcutaneous"),
  breast = list(label = "Breast cancer", tissue = "BreastMammaryTissue"),
  cad    = list(label = "CAD",    tissue = "ArteryCoronary"),
  ibd    = list(label = "IBD",    tissue = "ColonTransverse"),
  mdd    = list(label = "MDD",    tissue = "BrainCortex"),
  scz    = list(label = "SCZ",    tissue = "BrainCortex"),
  t2d    = list(label = "T2D",    tissue = "Pancreas"),
  alkaline_phosphatase = list(label = "Alkaline phosphatase", tissue = "Liver"),
  eosinophill_percentage = list(label = "Eosinophil %", tissue = "WholeBlood"),
  hdl_cholesterol = list(label = "HDL cholesterol", tissue = "Liver"),
  lipoprotein_a = list(label = "Lipoprotein A", tissue = "Liver"),
  lymphocyte_count = list(label = "Lymphocyte count", tissue = "WholeBlood"),
  mean_platelet_thrombocyte_volume = list(label = "Mean platelet volume", tissue = "WholeBlood"),
  monocyte_percentage = list(label = "Monocyte %", tissue = "WholeBlood"),
  platelet_crit = list(label = "Platelet crit", tissue = "WholeBlood"),
  urea = list(label = "Urea", tissue = "KidneyCortex"),
  vitamin_d = list(label = "Vitamin D", tissue = "SkinNotSunExposedSuprapubic")
)

# Read and parse files
parse_and_normalize <- function(f) {
  bn <- basename(f)
  mm <- regmatches(bn, regexec("^(.+?)_([A-Za-z0-9]+)_correlation_summary\\.csv$", bn))[[1]]
  if (length(mm) != 3) return(NULL)
  d <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  d[, trait := mm[2]]
  d[, tool_base := mm[3]]
  d
}

files <- list.files(MALACARDS_FOLDER, pattern = "_correlation_summary\\.csv$", recursive = TRUE, full.names = TRUE)
dat <- rbindlist(Filter(Negate(is.null), lapply(files, parse_and_normalize)), fill = TRUE)
if (nrow(dat) == 0) stop("No data found for tissue analysis.")

# Filter
dat <- dat %>%
  filter(subset %in% SUBSETS,
    method %in% METHODS3
  )

if (nrow(dat) == 0) stop("No filtered data for tissue analysis.")

# Extract trait and tool_base from the trait column
dat <- dat %>%
  mutate(
    trait_extracted = tolower(gsub("(_magma|_prset_birewire)$", "", trait)),
    tool_base_extracted = case_when(
      grepl("_magma$", trait) ~ "magma",
      grepl("_prset_birewire$", trait) ~ "prset",
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    trait = trait_extracted,
    tool_base = tool_base_extracted,
    everything(),
    -trait,        # remove old trait column
    -tool_base     # remove old tool_base column
  )

# ------------------------------------------------------------------------------
# Filter each trait for its relevant tissue
# ------------------------------------------------------------------------------
dat <- dat %>%
  mutate(
    tissue_metric_clean = gsub("_mean$", "", tissue_metric),
    relevant_tissue = sapply(trait, function(x) {
      if (!is.null(trait_map[[x]])) trait_map[[x]]$tissue else NA
    })
  ) %>%
  filter(tissue_metric_clean == relevant_tissue)

if (nrow(dat) == 0) stop("No data after tissue filtering.")

# ------------------------------------------------------------------------------
# Reshape and clean
# ------------------------------------------------------------------------------
## remove na rows
dat <- dat %>% filter(!is.na(kendall_tau) & !is.na(kendall_pvalue))

plot_dat <- dat %>%
  mutate(
    trait_label = trait
  ) %>%
  select(tool_base, trait_label, method, rho = spearman_rho, pval = correlation_pvalue) %>%
  pivot_wider(names_from = method, values_from = c(rho, pval))

# Convert comma-separated p-values to numeric (first value only)
plot_dat <- plot_dat %>%
  mutate(
    pval_BireWire_EmpPvalStdBeta_num = as.numeric(sapply(strsplit(as.character(pval_BireWire_EmpPvalStdBeta), ","), `[`, 1)),
    pval_PvalueBeta_num = as.numeric(sapply(strsplit(as.character(pval_PvalueBeta), ","), `[`, 1))
  )

# Annotate significance
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
    star_BireWire = annotate_significance(pval_BireWire_EmpPvalStdBeta_num),
    star_Pvalue = annotate_significance(pval_PvalueBeta_num)
  )

# Print for troubleshooting
print(plot_dat)

# Uncomment below to plot after verifying data
library(ggplot2)
library(ggalt)
library(patchwork)

trait_order <- plot_dat$trait_label

# Map trait_label to simplified trait names
plot_dat <- plot_dat %>%
  mutate(
    trait_simple = case_when(
      grepl("CAD|cad", trait_label) ~ "CAD",
      grepl("IBD|ibd", trait_label) ~ "IBD",
      grepl("Breast", trait_label, ignore.case = TRUE) ~ "Breast cancer",
      grepl("T2D|type2|diab", trait_label, ignore.case = TRUE) ~ "T2D",
      grepl("SCZ|schizo", trait_label, ignore.case = TRUE) ~ "SCZ",
      grepl("BMI|bmi", trait_label, ignore.case = TRUE) ~ "BMI",
      grepl("AD|alzheimer", trait_label, ignore.case = TRUE) ~ "AD",
      grepl("MDD|depress", trait_label, ignore.case = TRUE) ~ "MDD",
      TRUE ~ trait_label
    )
  )

trait_order <- unique(plot_dat$trait_simple)

# Subset to traits available in both magma and prset
traits_magma <- plot_dat %>% filter(tool_base == "magma") %>% pull(trait_simple)
traits_prset <- plot_dat %>% filter(tool_base == "prset") %>% pull(trait_simple)
traits_both <- intersect(traits_magma, traits_prset)

plot_dat <- plot_dat %>% filter(trait_simple %in% traits_both)

# Now split data
magma_df <- plot_dat %>% filter(grepl("magma", tool_base))
prset_df <- plot_dat %>% filter(grepl("prset", tool_base))

# Get shared x-axis limits
all_rho <- c(
  as.numeric(unlist(strsplit(as.character(plot_dat$rho_BireWire_EmpPvalStdBeta), ","))),
  as.numeric(unlist(strsplit(as.character(plot_dat$rho_PvalueBeta), ",")))
)
x_min <- min(all_rho, na.rm = TRUE)
x_max <- max(all_rho, na.rm = TRUE)

# Panel function with y = trait_simple
make_panel <- function(df, title, show_y) {
  ggplot(df, aes(
    x = rho_BireWire_EmpPvalStdBeta,
    xend = rho_PvalueBeta,
    y = factor(trait_simple, levels = trait_order)
  )) +
    geom_vline(xintercept = 0, color = "black", linewidth = 1.2) +
    geom_dumbbell(colour = "#a3c4dc", size = 3, colour_x = "#0e668b", colour_xend = "#dd1c77") +
    geom_text(aes(x = rho_BireWire_EmpPvalStdBeta, label = star_BireWire), vjust = -1, size = 5, color = "#0e668b") +
    geom_text(aes(x = rho_PvalueBeta, label = star_Pvalue), vjust = -1, size = 5, color = "#dd1c77") +
    labs(y = if (show_y) "Trait" else NULL, x = "Spearman rho", title = title) +
    xlim(x_min, x_max) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = if (show_y) element_text(size = 10, color = "black") else element_blank(),
      axis.title.y = if (show_y) element_text(color = "black") else element_blank(),
      axis.title.x = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5)
    )
}

# Create panels
panel_magma <- make_panel(magma_df, "MAGMA", show_y = TRUE)
panel_prset <- make_panel(prset_df, "PRSet", show_y = FALSE)

# Combine panels
p <- panel_magma + panel_prset + plot_layout(ncol = 2)

ggsave("dumbbell_tissue_twopanel.pdf", plot = p, width = 12, height = 4.2)

# Add mapped labels and tissues for downstream use
plot_dat <- plot_dat %>%
  mutate(
    trait_label = tolower(gsub("_magma_birewire|_birewire|_prset", "", trait_label)),
    trait_mapped_label = sapply(trait_label, function(x) if (!is.null(trait_map[[x]])) trait_map[[x]]$label else x),
    trait_mapped_tissue = sapply(trait_label, function(x) if (!is.null(trait_map[[x]])) trait_map[[x]]$tissue else NA)
  )