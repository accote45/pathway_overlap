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

# Function to read real results and merge with pre-calculated empirical p-values
read_real_results <- function(file_path, trait_name, random_method) {
  if(!file.exists(file_path)) {
    cat("Warning: File not found:", file_path, "\n")
    return(NULL)
  }
  
  # Read real MAGMA results
  data <- read.table(file_path, header = TRUE)
  data$name <- data$FULL_NAME
  data$Zscore <- data$BETA/data$SE
  data$Zscore_N <- (data$BETA / data$SE) / data$NGENES
  
  # Read pre-calculated empirical p-values
  emp_p_file <- paste0('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_random/',
                     'birewire', '/msigdbgenes/', trait_name, '/', 
                     trait_name, '_empirical_pvalues.csv')
  
  if(file.exists(emp_p_file)) {
    cat("Reading empirical p-values from:", emp_p_file, "\n")
    emp_p_data <- read.csv(emp_p_file)
    
    # Filter for the specific randomization method we want
    emp_p_data_filtered <- emp_p_data[emp_p_data$method == random_method, ]
    
    if(nrow(emp_p_data_filtered) > 0) {
      # Merge empirical p-values with real results
      data <- merge(data, emp_p_data_filtered[, !colnames(emp_p_data_filtered) %in% c("method")], 
                    by.x = "FULL_NAME", by.y = "FULL_NAME", all.x = TRUE)
    } else {
      cat("Warning: No empirical p-values found for method", random_method, "\n")
      cat("Using nominal p-values instead\n")
      data$empirical_pvalue <- data$P
      data$more_significant_count <- NA
      data$total_random <- NA
    }
    
    # Fill missing empirical p-values if any
    if(any(is.na(data$empirical_pvalue))) {
      cat("Warning: Some pathways missing empirical p-values, using nominal p-values for those\n")
      data$empirical_pvalue[is.na(data$empirical_pvalue)] <- data$P[is.na(data$empirical_pvalue)]
    }
  } else {
    cat("Warning: Empirical p-values file not found:", emp_p_file, "\n")
    cat("Using nominal p-values instead\n")
    data$empirical_pvalue <- data$P
    data$more_significant_count <- NA
    data$total_random <- NA
  }
  
  # Calculate ranks
  data$Zscore_rank <- rank(-data$Zscore)
  data$ZscoreN_rank <- rank(-data$Zscore_N)
  data <- data %>% arrange(empirical_pvalue, desc(BETA)) %>% mutate(EmpPBeta_rank = row_number())  
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
  
  # Read real gene set enrichment results with empirical p-values
  birewire_path <- paste0('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/results/', 
                         background, '/', method, '_real/', trait_name, 
                         '/setreal.empP.birewire.gsa.out')
  keeppath_path <- paste0('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/results/', 
                         background, '/', method, '_real/', trait_name, 
                         '/setreal.empP.keeppathsize.gsa.out')
  
  real_birewire <- read_real_results(birewire_path, trait_name, "birewire")
  real_keeppath <- read_real_results(keeppath_path, trait_name, "keeppathsize")
  
  if(is.null(real_birewire) || is.null(real_keeppath)) {
    cat("Skipping trait", trait_name, "- MAGMA results not found\n")
    next
  }
  
  # Merge OpenTargets data with pathway genes
  masterfin <- merge(genes_long, master2, by="targetId", all.x=TRUE)
  masterfin[is.na(masterfin)] <- 0
  
  # Calculate pathway-level scores, accounting for pathway size
  masterfin <- masterfin %>% 
    group_by(name) %>% 
    summarise(
      avg_score = mean(score, na.rm = TRUE),
      num_genes = n(),
      num_with_evidence = sum(score > 0, na.rm = TRUE),
      coverage = sum(score > 0, na.rm = TRUE) / n()
    ) %>% 
    as.data.frame()
  
  # Merge with real data
  merged_birewire <- merge(masterfin, real_birewire, by = "name")
  merged_keeppath <- merge(masterfin, real_keeppath, by = "name")
  
  # Calculate ranks
  merged_birewire$targetscore_rank <- rank(-merged_birewire$avg_score)
  merged_keeppath$targetscore_rank <- rank(-merged_keeppath$avg_score)
  
  # Calculate size-adjusted scores
  merged_birewire$coverage_adj_score <- merged_birewire$avg_score / sqrt(merged_birewire$coverage)
  merged_keeppath$coverage_adj_score <- merged_keeppath$avg_score / sqrt(merged_keeppath$coverage)
  
  # Store results
  all_results[[trait_name]] <- list(
    birewire = merged_birewire,
    keeppath = merged_keeppath,
    opentargets_count = nrow(master2)
  )
  
  # Perform correlation tests for this trait
  cat("Correlation tests for", trait_name, ":\n")
  
  # Nominal p-value correlations
  cat("=== Nominal P-value correlations ===\n")
  cor_zscore <- cor.test(merged_birewire$Zscore_rank, merged_birewire$targetscore_rank, method = "kendall")
  cor_zscore_n <- cor.test(merged_birewire$ZscoreN_rank, merged_birewire$targetscore_rank, method = "kendall")
  cor_pbeta <- cor.test(merged_birewire$PBeta_rank, merged_birewire$targetscore_rank, method = "kendall")
  
  cat("Z-score rank correlation (tau):", round(cor_zscore$estimate, 3), 
      "p-value:", format.pval(cor_zscore$p.value), "\n")
  cat("Z-score/N rank correlation (tau):", round(cor_zscore_n$estimate, 3), 
      "p-value:", format.pval(cor_zscore_n$p.value), "\n")
  cat("P+Beta rank correlation (tau):", round(cor_pbeta$estimate, 3), 
      "p-value:", format.pval(cor_pbeta$p.value), "\n")
  
  # Empirical p-value correlations
  cat("=== Empirical P-value correlations ===\n")
  cor_emppbeta <- cor.test(merged_birewire$EmpPBeta_rank, merged_birewire$targetscore_rank, method = "kendall")
  
  cat("Empirical P+Beta rank correlation (tau):", round(cor_emppbeta$estimate, 3), 
      "p-value:", format.pval(cor_emppbeta$p.value), "\n")
  
  # Size-stratified analysis for this trait
  cat("=== Size-stratified analysis ===\n")
  # Create size bins for stratification
  merged_birewire$size_bin <- cut(merged_birewire$NGENES, 
                                breaks = c(0, 50, 100, 200, 500, Inf),
                                labels = c("<50", "50-100", "100-200", "200-500", ">500"))
  
  # Correlation by size bin
  size_cors <- merged_birewire %>%
    group_by(size_bin) %>%
    summarize(
      n = n(),
      nom_p_cor = cor(-log10(P), avg_score, method = "spearman"),
      emp_p_cor = cor(-log10(empirical_pvalue), avg_score, method = "spearman")
    ) %>%
    filter(n >= 5)  # Only include bins with enough data
  
  print(size_cors)
}

# Summary across all traits
cat("\n=== Summary across all traits ===\n")
cat("Successfully processed traits:", paste(names(all_results), collapse = ", "), "\n")
cat("Total traits processed:", length(all_results), "\n")

# Create summary plots for all traits
if(length(all_results) > 0) {
  # Combine data from all traits for summary visualization
  summary_data <- rbindlist(lapply(names(all_results), function(trait) {
    birewire_data <- all_results[[trait]]$birewire
    birewire_data$trait <- trait
    birewire_data$method <- "birewire"
    
    keeppath_data <- all_results[[trait]]$keeppath  
    keeppath_data$trait <- trait
    keeppath_data$method <- "keeppathsize"
    
    rbind(birewire_data[, c("trait", "method", "avg_score", "coverage_adj_score", "Zscore", "P", "empirical_pvalue", "NGENES")],
          keeppath_data[, c("trait", "method", "avg_score", "coverage_adj_score", "Zscore", "P", "empirical_pvalue", "NGENES")])
  }))
  
  # Create summary plots comparing nominal and empirical p-values
  library(patchwork)
  pdf('opentargets_summary_all_traits.pdf', width=15, height=12)
  
  # Plot 1: Nominal P vs OpenTargets score
  p1 <- ggplot(summary_data, aes(x=-log10(P), y=avg_score, color=trait)) + 
    geom_point(alpha=0.6, aes(size=NGENES)) + 
    facet_wrap(~method) +
    theme_classic() + 
    geom_smooth(method="lm", se=TRUE, color="black") + 
    xlab('-log10(Nominal P-value)') + 
    ylab('Open Targets score') +
    ggtitle('Nominal P-value vs OpenTargets Score by Method') +
    theme(legend.position="bottom")
  
  # Plot 2: Empirical P vs OpenTargets score  
  p2 <- ggplot(summary_data, aes(x=-log10(empirical_pvalue), y=avg_score, color=trait)) + 
    geom_point(alpha=0.6, aes(size=NGENES)) + 
    facet_wrap(~method) +
    theme_classic() + 
    geom_smooth(method="lm", se=TRUE, color="black") + 
    xlab('-log10(Empirical P-value)') + 
    ylab('Open Targets score') +
    ggtitle('Empirical P-value vs OpenTargets Score by Method') +
    theme(legend.position="bottom")

  # Plot 3: Size-adjusted analysis
  p3 <- ggplot(summary_data, aes(x=-log10(P), y=coverage_adj_score, color=trait)) + 
    geom_point(alpha=0.6) + 
    facet_wrap(~method) +
    theme_classic() + 
    geom_smooth(method="lm", se=TRUE, color="black") + 
    xlab('-log10(Nominal P-value)') + 
    ylab('Coverage-adjusted Open Targets score') +
    ggtitle('Nominal P-value vs Coverage-adjusted OpenTargets Score')
  
  p4 <- ggplot(summary_data, aes(x=-log10(empirical_pvalue), y=coverage_adj_score, color=trait)) + 
    geom_point(alpha=0.6) + 
    facet_wrap(~method) +
    theme_classic() + 
    geom_smooth(method="lm", se=TRUE, color="black") + 
    xlab('-log10(Empirical P-value)') + 
    ylab('Coverage-adjusted Open Targets score') +
    ggtitle('Empirical P-value vs Coverage-adjusted OpenTargets Score')
  
  # Size-stratified analysis
  summary_data$size_bin <- cut(summary_data$NGENES, 
                              breaks = c(0, 50, 100, 200, 500, Inf),
                              labels = c("<50", "50-100", "100-200", "200-500", ">500"))
  
  p5 <- ggplot(summary_data, aes(x=-log10(P), y=avg_score, color=trait)) + 
    geom_point(alpha=0.6) + 
    facet_grid(method ~ size_bin) +
    theme_classic() + 
    geom_smooth(method="lm", se=FALSE, color="black") + 
    xlab('-log10(Nominal P-value)') + 
    ylab('Open Targets score') +
    ggtitle('Nominal P-value vs OpenTargets Score by Pathway Size') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  p6 <- ggplot(summary_data, aes(x=-log10(empirical_pvalue), y=avg_score, color=trait)) + 
    geom_point(alpha=0.6) + 
    facet_grid(method ~ size_bin) +
    theme_classic() + 
    geom_smooth(method="lm", se=FALSE, color="black") + 
    xlab('-log10(Empirical P-value)') + 
    ylab('Open Targets score') +
    ggtitle('Empirical P-value vs OpenTargets Score by Pathway Size') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  # Print all plots
  print(p1 / p2)
  print(p3 / p4)
  print(p5)
  print(p6)
  dev.off()
  
  # Calculate overall correlation improvement
  cat("\n=== Overall correlation improvement with empirical p-values ===\n")
  improvement_stats <- summary_data %>%
    group_by(trait, method) %>%
    summarize(
      nom_p_cor = cor(-log10(P), avg_score, method = "spearman", use = "pairwise.complete.obs"),
      emp_p_cor = cor(-log10(empirical_pvalue), avg_score, method = "spearman", use = "pairwise.complete.obs"),
      improvement = emp_p_cor - nom_p_cor,
      n_pathways = n()
    ) %>%
    arrange(desc(improvement))
  
  print(improvement_stats)
  
  # Write summary tables to files
  write.csv(improvement_stats, "improvement_by_trait_method.csv", row.names = FALSE)
  
  # Create size-stratified correlation table
  size_strat <- summary_data %>%
    group_by(method, size_bin) %>%
    summarize(
      nom_p_cor = cor(-log10(P), avg_score, method = "spearman", use = "pairwise.complete.obs"),
      emp_p_cor = cor(-log10(empirical_pvalue), avg_score, method = "spearman", use = "pairwise.complete.obs"),
      improvement = emp_p_cor - nom_p_cor,
      n_pathways = n()
    ) %>%
    arrange(method, size_bin)
  
  write.csv(size_strat, "improvement_by_size_bin.csv", row.names = FALSE)
}








