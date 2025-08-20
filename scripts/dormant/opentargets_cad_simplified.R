# Simplified script focusing only on CAD trait

library(jsonlite)
library(tidyverse)
library(data.table)
library(GSA)
library(patchwork)

# Define CAD disease ID
cad_id <- "EFO_0001645"

# Set parameters
background <- "pathwaydb_enrichment_msigdbgenes"
method <- "magma"

# Set working directory
setwd('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/results/drugtarget_test/associationByDatatypeDirect')

# Read gene set data
cat("Reading gene set data...\n")
dat <- GSA.read.gmt('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/pathway_db/msigdb/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt')
path.list <- dat$genesets
names(path.list) <- dat$geneset.names

# Convert the path list to a data table
genes_long <- rbindlist(lapply(names(path.list), function(name) {
  data.table(value = path.list[[name]], name = name)
}))
genes_long$targetId <- genes_long$value

# Process OpenTargets JSON files for CAD
cat("Processing OpenTargets data for CAD...\n")
files <- list.files(pattern="*.json")

# Simple function to read JSON files
process_cad_data <- function(files) {
  cat("Reading JSON files for CAD...\n")
  
  dat <- lapply(files, function(file) {
    tryCatch({
      json_data <- stream_in(file(file), verbose = FALSE)
      json_data <- as.data.frame(json_data)
      json_data %>%
        filter(diseaseId == cad_id & datatypeId %in% c('animal_model', 'known_drug', 'literature'))
    }, error = function(e) {
      cat("Error reading file:", file, "\n")
      return(NULL)
    })
  })
  
  combined_dat <- do.call(rbind, dat[!sapply(dat, is.null)])
  
  if(nrow(combined_dat) == 0) {
    stop("No OpenTargets data found for CAD")
  }
  
  return(combined_dat)
}

# Read CAD OpenTargets data
cad_data <- process_cad_data(files)

# Select max evidence count for each target
cad_targets <- cad_data %>%
  group_by(targetId) %>%
  slice_max(evidenceCount, with_ties=FALSE) %>%
  ungroup()

# Read real MAGMA results for CAD
cat("Reading MAGMA results for CAD...\n")
real_results <- paste0('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/',method,'_real/cad/cad_real_set.gsa.out')

# Read MAGMA results and empirical p-values
read_cad_results <- function(file_path, random_method) {
  # Read real MAGMA results
  data <- read.table(file_path, header = TRUE)
  data$name <- data$FULL_NAME
  data$Zscore <- data$BETA/data$SE
  data$Zscore_N <- (data$BETA / data$SE) / data$NGENES
  
  # Read pre-calculated empirical p-values
  emp_p_file <- paste0('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_random/birewire/msigdbgenes/cad/cad_empirical_pvalues.csv')
  
  emp_p_data <- read.csv(emp_p_file)
  
  # Filter for the specific randomization method
  emp_p_data_filtered <- emp_p_data[emp_p_data$method == random_method, ]
  
  # Merge empirical p-values with real results
  data <- merge(data, emp_p_data_filtered[, !colnames(emp_p_data_filtered) %in% c("method")], 
                by.x = "FULL_NAME", by.y = "FULL_NAME", all.x = TRUE)
  
  return(data)
}

# Get CAD results for both methods
cad_birewire <- read_cad_results(real_results, "birewire")
cad_keeppath <- read_cad_results(real_results, "keeppathsize")

# Merge OpenTargets data with pathway genes
cat("Merging OpenTargets data with pathways...\n")
masterfin <- merge(genes_long, cad_targets, by="targetId", all.x=TRUE)
masterfin[is.na(masterfin)] <- 0

# Calculate pathway-level scores
masterfin <- masterfin %>% 
  group_by(name) %>% 
  summarise(
    avg_score = mean(score, na.rm = TRUE),
    num_genes = n(),
    num_with_evidence = sum(score > 0, na.rm = TRUE),
    coverage = sum(score > 0, na.rm = TRUE) / n()
  ) %>% 
  as.data.frame()

# Merge with MAGMA results
merged_birewire <- merge(masterfin, cad_birewire, by = "name")
merged_keeppath <- merge(masterfin, cad_keeppath, by = "name")

# Calculate coverage-adjusted scores
merged_birewire$coverage_adj_score <- merged_birewire$avg_score / sqrt(merged_birewire$coverage)
merged_keeppath$coverage_adj_score <- merged_keeppath$avg_score / sqrt(merged_keeppath$coverage)

# Add size bins
merged_birewire$size_bin <- cut(merged_birewire$NGENES, 
                               breaks = c(0, 50, 100, 200, 500, Inf),
                               labels = c("<50", "50-100", "100-200", "200-500", ">500"))
merged_keeppath$size_bin <- cut(merged_keeppath$NGENES, 
                               breaks = c(0, 50, 100, 200, 500, Inf),
                               labels = c("<50", "50-100", "100-200", "200-500", ">500"))

# Combine data for plotting
plot_data <- rbind(
  transform(merged_birewire, method = "birewire"),
  transform(merged_keeppath, method = "keeppathsize")
)

# Calculate correlations
cat("\n=== Correlation statistics ===\n")
birewire_nom_cor <- cor(-log10(merged_birewire$P), merged_birewire$avg_score, method = "spearman")
birewire_emp_cor <- cor(-log10(merged_birewire$empirical_pvalue), merged_birewire$avg_score, method = "spearman")
keeppath_nom_cor <- cor(-log10(merged_keeppath$P), merged_keeppath$avg_score, method = "spearman")
keeppath_emp_cor <- cor(-log10(merged_keeppath$empirical_pvalue), merged_keeppath$avg_score, method = "spearman")

cat("Birewire - Nominal P correlation:", round(birewire_nom_cor, 3), "\n")
cat("Birewire - Empirical P correlation:", round(birewire_emp_cor, 3), "\n")
cat("Birewire - Improvement:", round(birewire_emp_cor - birewire_nom_cor, 3), "\n\n")

cat("Keeppathsize - Nominal P correlation:", round(keeppath_nom_cor, 3), "\n")
cat("Keeppathsize - Empirical P correlation:", round(keeppath_emp_cor, 3), "\n")
cat("Keeppathsize - Improvement:", round(keeppath_emp_cor - keeppath_nom_cor, 3), "\n")

# Size-stratified analysis
cat("\n=== Size-stratified analysis ===\n")
size_stats <- plot_data %>%
  group_by(method, size_bin) %>%
  summarize(
    n = n(),
    nom_p_cor = cor(-log10(P), avg_score, method = "spearman", use = "pairwise.complete.obs"),
    emp_p_cor = cor(-log10(empirical_pvalue), avg_score, method = "spearman", use = "pairwise.complete.obs"),
    improvement = emp_p_cor - nom_p_cor
  ) %>%
  filter(n >= 5)  # Only include bins with enough data

print(size_stats)

# Create comparison plots
cat("\nCreating plots...\n")
pdf('cad_opentargets_analysis.pdf', width=12, height=10)

# Plot 1: Nominal vs Empirical P-values with OpenTargets score
p1 <- ggplot(plot_data, aes(x=-log10(P), y=avg_score, color=method)) + 
  geom_point(alpha=0.6, aes(size=NGENES)) + 
  geom_smooth(method="lm", se=TRUE) + 
  labs(x='-log10(Nominal P-value)', 
       y='Open Targets score',
       title='CAD: Nominal P-value vs OpenTargets Score')

p2 <- ggplot(plot_data, aes(x=-log10(empirical_pvalue), y=avg_score, color=method)) + 
  geom_point(alpha=0.6, aes(size=NGENES)) + 
  geom_smooth(method="lm", se=TRUE) + 
  labs(x='-log10(Empirical P-value)', 
       y='Open Targets score',
       title='CAD: Empirical P-value vs OpenTargets Score')

# Plot by pathway size
p3 <- ggplot(plot_data, aes(x=-log10(P), y=avg_score, color=method)) + 
  geom_point(alpha=0.6) + 
  geom_smooth(method="lm", se=FALSE) + 
  facet_wrap(~ size_bin) +
  labs(x='-log10(Nominal P-value)', y='Open Targets score',
       title='CAD: Nominal P-value vs OpenTargets Score by Pathway Size')

p4 <- ggplot(plot_data, aes(x=-log10(empirical_pvalue), y=avg_score, color=method)) + 
  geom_point(alpha=0.6) + 
  geom_smooth(method="lm", se=FALSE) + 
  facet_wrap(~ size_bin) +
  labs(x='-log10(Empirical P-value)', y='Open Targets score',
       title='CAD: Empirical P-value vs OpenTargets Score by Pathway Size')

# Print all plots
print(p1 / p2)
print(p3)
print(p4)
dev.off()

# Write results to CSV
write.csv(plot_data, "cad_opentargets_with_empirical_pvalues.csv", row.names = FALSE)
write.csv(size_stats, "cad_size_stratified_correlations.csv", row.names = FALSE)

cat("\nAnalysis complete. Results saved to cad_opentargets_analysis.pdf\n")

# After your existing size-stratified analysis, add this more granular analysis:

# Create finer-grained size bins for more detailed analysis
cat("\n=== Fine-grained size bin analysis ===\n")

# Define more specific size bins
plot_data$fine_size_bin <- cut(plot_data$NGENES, 
                              breaks = c(0, 20, 50, 75, 100, 125, 150, 175, 200, 250, 300, 400, 500, Inf),
                              labels = c("<20", "20-50", "51-75", "76-100", "101-125", "126-150", 
                                        "151-175", "176-200", "201-250", "251-300", "301-400", 
                                        "401-500", ">500"))

# Calculate correlations for these finer bins
fine_size_stats <- plot_data %>%
  group_by(method, fine_size_bin) %>%
  summarize(
    n = n(),
    nom_p_cor = cor(-log10(P), avg_score, method = "spearman", use = "pairwise.complete.obs"),
    emp_p_cor = cor(-log10(empirical_pvalue), avg_score, method = "spearman", use = "pairwise.complete.obs"),
    improvement = emp_p_cor - nom_p_cor,
    median_ngenes = median(NGENES)
  ) %>%
  filter(n >= 3)  # Include bins with at least 3 pathways

print(fine_size_stats)

# Create a plot showing improvement by pathway size bin
p5 <- ggplot(fine_size_stats, aes(x=fine_size_bin, y=improvement, fill=method)) + 
  geom_bar(stat="identity", position="dodge") +
  geom_text(aes(label=n), position=position_dodge(width=0.9), vjust=-0.5, size=3) +
  labs(x='Pathway Size Bin', y='Correlation Improvement (Empirical - Nominal)',
       title='CAD: Correlation Improvement by Fine-Grained Pathway Size') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Create detailed plots for specific size ranges of interest
# Example: focused analysis on pathways with 175-225 genes
focused_range_pathways <- plot_data %>%
  filter(NGENES >= 175 & NGENES <= 225)

if(nrow(focused_range_pathways) >= 10) {  # Only create if we have enough data points
  cat("\n=== Focused analysis on pathways with 175-225 genes ===\n")
  cat("Number of pathways in this range:", nrow(focused_range_pathways), "\n")
  
  # Calculate correlations for this specific range
  focused_cor <- focused_range_pathways %>%
    group_by(method) %>%
    summarize(
      n = n(),
      nom_p_cor = cor(-log10(P), avg_score, method = "spearman"),
      emp_p_cor = cor(-log10(empirical_pvalue), avg_score, method = "spearman"),
      improvement = emp_p_cor - nom_p_cor
    )
  
  print(focused_cor)
  
  # Create plots for these pathways
  p6 <- ggplot(focused_range_pathways, aes(x=-log10(P), y=avg_score, color=method)) + 
    geom_point(alpha=0.8, aes(size=NGENES)) + 
    geom_smooth(method="lm", se=TRUE) + 
    labs(x='-log10(Nominal P-value)', y='Open Targets score',
         title='CAD: Nominal P-value vs OpenTargets Score (175-225 genes)')
  
  p7 <- ggplot(focused_range_pathways, aes(x=-log10(empirical_pvalue), y=avg_score, color=method)) + 
    geom_point(alpha=0.8, aes(size=NGENES)) + 
    geom_smooth(method="lm", se=TRUE) + 
    labs(x='-log10(Empirical P-value)', y='Open Targets score',
         title='CAD: Empirical P-value vs OpenTargets Score (175-225 genes)')
}

# Update the PDF with these additional plots
cat("\nUpdating plots with detailed size analysis...\n")
pdf('cad_opentargets_analysis_detailed.pdf', width=15, height=12)

# Print original plots
print(p1 / p2)
print(p3)
print(p4)

# Print new plots
print(p5)
if(exists("p6") && exists("p7")) {
  print(p6 / p7)
}

# Create a detailed heatmap of improvement by bin size
size_matrix <- matrix(NA, nrow = length(unique(fine_size_stats$method)), 
                     ncol = length(unique(fine_size_stats$fine_size_bin)))
rownames(size_matrix) <- unique(fine_size_stats$method)
colnames(size_matrix) <- unique(fine_size_stats$fine_size_bin)

for(i in 1:nrow(fine_size_stats)) {
  row <- which(rownames(size_matrix) == fine_size_stats$method[i])
  col <- which(colnames(size_matrix) == fine_size_stats$fine_size_bin[i])
  size_matrix[row, col] <- fine_size_stats$improvement[i]
}

# Plot heatmap
library(reshape2)
melted_matrix <- melt(size_matrix)
colnames(melted_matrix) <- c("Method", "Size_Bin", "Improvement")

p8 <- ggplot(melted_matrix, aes(x=Size_Bin, y=Method, fill=Improvement)) + 
  geom_tile() + 
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title="Improvement in Correlation by Method and Pathway Size",
       x="Pathway Size Bin", y="Randomization Method")

print(p8)

# Create a scatter plot showing improvement vs. pathway size
improvement_by_size <- plot_data %>%
  group_by(method) %>%
  mutate(size_group = cut(NGENES, breaks=10)) %>%
  group_by(method, size_group) %>%
  summarize(
    mean_size = mean(NGENES),
    mean_improvement = mean(abs(-log10(empirical_pvalue) - -log10(P))),
    n = n()
  )

p9 <- ggplot(improvement_by_size, aes(x=mean_size, y=mean_improvement, color=method, size=n)) +
  geom_point(alpha=0.7) +
  geom_smooth(method="loess", se=TRUE) +
  labs(x="Mean Pathway Size", y="Mean Absolute Difference in -log10(p)",
       title="Impact of Pathway Size on P-value Differences") +
  theme_minimal()

print(p9)

dev.off()

# Write detailed results to CSV files
write.csv(fine_size_stats, "cad_fine_grained_size_correlations.csv", row.names = FALSE)

# Create additional focused analyses for different size ranges
size_ranges <- list(
  "small" = c(10, 50),
  "medium" = c(90, 110),
  "medium_large" = c(190, 210),
  "large" = c(400, 500)
)

# Function to analyze a specific size range
analyze_size_range <- function(data, min_size, max_size, label) {
  subset <- data %>% filter(NGENES >= min_size & NGENES <= max_size)
  
  if(nrow(subset) < 5) {
    cat("\nNot enough pathways in size range", min_size, "-", max_size, "genes\n")
    return(NULL)
  }
  
  cat("\n=== Analysis for pathways with", min_size, "-", max_size, "genes ===\n")
  cat("Number of pathways:", nrow(subset), "\n")
  
  # Calculate correlations
  result <- subset %>%
    group_by(method) %>%
    summarize(
      n = n(),
      median_size = median(NGENES),
      nom_p_cor = cor(-log10(P), avg_score, method = "spearman"),
      emp_p_cor = cor(-log10(empirical_pvalue), avg_score, method = "spearman"),
      improvement = emp_p_cor - nom_p_cor
    )
  
  print(result)
  write.csv(result, paste0("cad_size_range_", label, ".csv"), row.names = FALSE)
  
  return(result)
}

# Run analyses for each size range
size_range_results <- lapply(names(size_ranges), function(label) {
  range <- size_ranges[[label]]
  analyze_size_range(plot_data, range[1], range[2], label)
})

cat("\nDetailed size analysis complete. Results saved to files.\n")