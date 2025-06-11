# T2D = MONDO_0005148
# CAD = EFO_0001645
# BMI = EFO_0004340
# AD = MONDO_0004975
# Depression = MONDO_0002050
# SCZ = MONDO_0005090
# IBD = EFO_0000555
# breast cancer = MONDO_0007254
# HDL = EFO_0004612
# lymphocyte = EFO_0004587
# platelet crit = EFO_0007985
# alkaline phosphatase = EFO_0004533

library(jsonlite)
library(tidyverse)
library(data.table)
library(GSA)

# Define trait mapping
trait_mapping <- list(
  "t2d" = "MONDO_0005148",
  "cad" = "EFO_0001645", 
  "bmi" = "EFO_0004340",
  "ad" = "MONDO_0004975",
  "major_depression" = "MONDO_0002050",
  "scz" = "MONDO_0005090",
  "ibd" = "EFO_0000555",
  "breast" = "MONDO_0007254",
  "hdl_cholesterol" = "EFO_0004612",
  "lymphocyte_count" = "EFO_0004587",
  "platelet_crit" = "EFO_0007985",
  "alkaline_phosphatase" = "EFO_0004533"
)

# Set parameters
background <- "pathwaydb_enrichment_msigdbgenes"
method <- "magma"

setwd('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/results/drugtarget_test/associationByDatatypeDirect')

# Updated function to process open targets JSON files for specific disease ID
process_json_files <- function(files, disease_id) {
  cat("Processing JSON files for disease ID:", disease_id, "\n")
  
  dat <- lapply(files, function(file) {
    json_data <- stream_in(file(file), verbose = FALSE)
    json_data <- as.data.frame(json_data)
    json_data %>%
      filter(diseaseId == disease_id & datatypeId %in% c('animal_model', 'known_drug', 'literature'))
  })
  
  combined_dat <- do.call(rbind, dat)
  
  if(nrow(combined_dat) == 0) {
    cat("Warning: No data found for disease ID:", disease_id, "\n")
    return(NULL)
  }
  
  return(combined_dat)
}

# Read gene set data once (used for all traits)
dat <- GSA.read.gmt('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/pathway_db/msigdb/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt')
path.list <- dat$genesets
names(path.list) <- dat$geneset.names

# Convert the path list to a data table
genes_long <- rbindlist(lapply(names(path.list), function(name) {
  data.table(value = path.list[[name]], name = name)
}))
genes_long <- as.data.frame(genes_long)
genes_long$targetId <- genes_long$value

# Get JSON files
files <- list.files(pattern="*.json")

# Function to read real results
read_real_results <- function(file_path) {
  if(!file.exists(file_path)) {
    cat("Warning: File not found:", file_path, "\n")
    return(NULL)
  }
  
  data <- read.table(file_path, header = TRUE)
  data$name <- data$FULL_NAME
  data$Zscore <- data$BETA/data$SE
  data$Zscore_N <- (data$BETA / data$SE) / data$NGENES
  data$Zscore_rank <- rank(-data$Zscore)
  data$ZscoreN_rank <- rank(-data$Zscore_N)
  data <- data %>% arrange(EMPIRICAL_P, desc(BETA)) %>% mutate(EmpPBeta_rank = row_number())  
  data <- data %>% arrange(P, desc(BETA)) %>% mutate(PBeta_rank = row_number())  
  return(data)
}

# Main processing loop for all traits
all_results <- list()

for (trait_name in names(trait_mapping)) {
  cat("\n=== Processing trait:", trait_name, "===\n")
  
  disease_id <- trait_mapping[[trait_name]]
  
  # Process JSON files for this trait
  master <- process_json_files(files, disease_id)
  
  if(is.null(master)) {
    cat("Skipping trait", trait_name, "- no OpenTargets data found\n")
    next
  }
  
  # Select max evidence count for each target
  master2 <- master %>%
    group_by(targetId) %>%
    slice_max(evidenceCount, with_ties=FALSE) %>%
    ungroup()
  
  # Read real gene set enrichment results
  birewire_path <- paste0('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/results/', 
                         background, '/', method, '_real/', trait_name, 
                         '/setreal.empP.birewire.gsa.out')
  keeppath_path <- paste0('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/results/', 
                         background, '/', method, '_real/', trait_name, 
                         '/setreal.empP.keeppathsize.gsa.out')
  
  real_birewire <- read_real_results(birewire_path)
  real_keeppath <- read_real_results(keeppath_path)
  
  if(is.null(real_birewire) || is.null(real_keeppath)) {
    cat("Skipping trait", trait_name, "- MAGMA results not found\n")
    next
  }
  
  # Merge OpenTargets data with pathway genes
  masterfin <- merge(genes_long, master2, by="targetId", all.x=TRUE)
  masterfin[is.na(masterfin)] <- 0
  masterfin <- masterfin %>% 
    group_by(name) %>% 
    summarise(avg_score = mean(score, na.rm = TRUE)) %>% 
    as.data.frame()
  
  # Merge with real data
  merged_birewire <- merge(masterfin, real_birewire, by = "name")
  merged_keeppath <- merge(masterfin, real_keeppath, by = "name")
  
  # Calculate ranks
  merged_birewire$targetscore_rank <- rank(-merged_birewire$avg_score)
  merged_keeppath$targetscore_rank <- rank(-merged_keeppath$avg_score)
  
  # Store results
  all_results[[trait_name]] <- list(
    birewire = merged_birewire,
    keeppath = merged_keeppath,
    opentargets_count = nrow(master2)
  )
  
  # Perform correlation tests for this trait
  cat("Correlation tests for", trait_name, ":\n")
  
  cor_zscore <- cor.test(merged_birewire$Zscore_rank, merged_birewire$targetscore_rank, method = "kendall")
  cor_zscore_n <- cor.test(merged_birewire$ZscoreN_rank, merged_birewire$targetscore_rank, method = "kendall")
  
  cat("Z-score rank correlation (tau):", round(cor_zscore$estimate, 3), 
      "p-value:", format.pval(cor_zscore$p.value), "\n")
  cat("Z-score/N rank correlation (tau):", round(cor_zscore_n$estimate, 3), 
      "p-value:", format.pval(cor_zscore_n$p.value), "\n")
}

# Summary across all traits
cat("\n=== Summary across all traits ===\n")
cat("Successfully processed traits:", paste(names(all_results), collapse = ", "), "\n")
cat("Total traits processed:", length(all_results), "\n")

# Create summary plot for all traits
if(length(all_results) > 0) {
  # Combine data from all traits for summary visualization
  summary_data <- rbindlist(lapply(names(all_results), function(trait) {
    birewire_data <- all_results[[trait]]$birewire
    birewire_data$trait <- trait
    birewire_data$method <- "birewire"
    
    keeppath_data <- all_results[[trait]]$keeppath  
    keeppath_data$trait <- trait
    keeppath_data$method <- "keeppathsize"
    
    rbind(birewire_data[, c("trait", "method", "avg_score", "Zscore", "P", "EMPIRICAL_P")],
          keeppath_data[, c("trait", "method", "avg_score", "Zscore", "P", "EMPIRICAL_P")])
  }))
  
  # Create summary plots
  library(patchwork)
  pdf('opentargets_summary_all_traits.pdf', width=15, height=10)
  
  p1 <- ggplot(summary_data, aes(x=P, y=avg_score, color=trait)) + 
    geom_point(alpha=0.6) + 
    facet_wrap(~method) +
    theme_classic() + 
    geom_smooth(method="lm", se=TRUE) + 
    xlab('P-value') + 
    ylab('Open Targets score') +
    ggtitle('P-value vs OpenTargets Score by Method')
  
  p2 <- ggplot(summary_data, aes(x=EMPIRICAL_P, y=avg_score, color=trait)) + 
    geom_point(alpha=0.6) + 
    facet_wrap(~method) +
    theme_classic() + 
    geom_smooth(method="lm", se=TRUE) + 
    xlab('Empirical P-value') + 
    ylab('Open Targets score') +
    ggtitle('Empirical P-value vs OpenTargets Score by Method')
  
  print(p1 / p2)
  dev.off()
}








