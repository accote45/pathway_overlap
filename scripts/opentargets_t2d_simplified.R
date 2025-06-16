# Simplified script focusing only on T2D trait

library(jsonlite)
library(tidyverse)
library(data.table)
library(GSA)
library(patchwork)

# Define T2D disease ID
t2d_id <- "MONDO_0005148"

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

# Process OpenTargets JSON files for T2D
cat("Processing OpenTargets data for T2D...\n")
files <- list.files(pattern="*.json")

# Simple function to read JSON files
process_t2d_data <- function(files) {
  cat("Reading JSON files for T2D...\n")
  
  dat <- lapply(files, function(file) {
    tryCatch({
      json_data <- stream_in(file(file), verbose = FALSE)
      json_data <- as.data.frame(json_data)
      json_data %>%
        filter(diseaseId == t2d_id & datatypeId %in% c('animal_model', 'known_drug', 'literature'))
    }, error = function(e) {
      cat("Error reading file:", file, "\n")
      return(NULL)
    })
  })
  
  combined_dat <- do.call(rbind, dat[!sapply(dat, is.null)])
  
  if(nrow(combined_dat) == 0) {
    stop("No OpenTargets data found for T2D")
  }
  
  return(combined_dat)
}

# Read T2D OpenTargets data
t2d_data <- process_t2d_data(files)

# Select max evidence count for each target
t2d_targets <- t2d_data %>%
  group_by(targetId) %>%
  slice_max(evidenceCount, with_ties=FALSE) %>%
  ungroup()

# Read real MAGMA results for T2D
cat("Reading MAGMA results for T2D...\n")
real_results <- paste0('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/',method,'_real/t2d/t2d_real_set.gsa.out')

# Read MAGMA results and empirical p-values
read_t2d_results <- function(file_path, random_method) {
  # Read real MAGMA results
  data <- read.table(file_path, header = TRUE)
  data$name <- data$FULL_NAME
  data$Zscore <- data$BETA/data$SE
  data$Zscore_N <- (data$BETA / data$SE) / data$NGENES
  
  # Read pre-calculated empirical p-values
  emp_p_file <- paste0('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_random/birewire/msigdbgenes/t2d/t2d_empirical_pvalues.csv')
  
  if(file.exists(emp_p_file)) {
    emp_p_data <- read.csv(emp_p_file)
    
    # Filter for the specific randomization method
    emp_p_data_filtered <- emp_p_data[emp_p_data$method == random_method, ]
    
    if(nrow(emp_p_data_filtered) > 0) {
      # Merge empirical p-values with real results
      data <- merge(data, emp_p_data_filtered[, !colnames(emp_p_data_filtered) %in% c("method")], 
                  by.x = "FULL_NAME", by.y = "FULL_NAME", all.x = TRUE)
      
      # Fill missing empirical p-values if any
      if(any(is.na(data$empirical_pvalue))) {
        cat("Warning: Some pathways missing empirical p-values, using nominal p-values for those\n")
        data$empirical_pvalue[is.na(data$empirical_pvalue)] <- data$P[is.na(data$empirical_pvalue)]
      }
    } else {
      cat("Warning: No empirical p-values found for method", random_method, "\n")
      cat("Using nominal p-values instead\n")
      data$empirical_pvalue <- data$P
      data$more_significant_count <- NA
      data$total_random <- NA
    }
  } else {
    cat("Warning: Empirical p-values file not found:", emp_p_file, "\n")
    cat("Using nominal p-values instead\n")
    data$empirical_pvalue <- data$P
    data$more_significant_count <- NA
    data$total_random <- NA
  }
  
  return(data)
}

# Get T2D results for both methods
t2d_birewire <- read_t2d_results(real_results, "birewire")
t2d_keeppath <- read_t2d_results(real_results, "keeppathsize")

# Merge OpenTargets data with pathway genes
cat("Merging OpenTargets data with pathways...\n")
masterfin <- merge(genes_long, t2d_targets, by="targetId", all.x=TRUE)
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
merged_birewire <- merge(masterfin, t2d_birewire, by = "name")
merged_keeppath <- merge(masterfin, t2d_keeppath, by = "name")

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
cat("\n=== Correlation statistics for T2D ===\n")
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
pdf('t2d_opentargets_analysis.pdf', width=12, height=10)

# Plot 1: Nominal vs Empirical P-values with OpenTargets score
p1 <- ggplot(plot_data, aes(x=-log10(P), y=avg_score, color=method)) + 
  geom_point(alpha=0.6, aes(size=NGENES)) + 
  geom_smooth(method="lm", se=TRUE) + 
  labs(x='-log10(Nominal P-value)', 
       y='Open Targets score',
       title='T2D: Nominal P-value vs OpenTargets Score')

p2 <- ggplot(plot_data, aes(x=-log10(empirical_pvalue), y=avg_score, color=method)) + 
  geom_point(alpha=0.6, aes(size=NGENES)) + 
  geom_smooth(method="lm", se=TRUE) + 
  labs(x='-log10(Empirical P-value)', 
       y='Open Targets score',
       title='T2D: Empirical P-value vs OpenTargets Score')

# Plot 3: Analysis with coverage-adjusted scores
p3 <- ggplot(plot_data, aes(x=-log10(P), y=coverage_adj_score, color=method)) + 
  geom_point(alpha=0.6, aes(size=NGENES)) + 
  geom_smooth(method="lm", se=TRUE) + 
  labs(x='-log10(Nominal P-value)', 
       y='Coverage-adjusted score',
       title='T2D: Nominal P-value vs Coverage-adjusted OpenTargets Score')

p4 <- ggplot(plot_data, aes(x=-log10(empirical_pvalue), y=coverage_adj_score, color=method)) + 
  geom_point(alpha=0.6, aes(size=NGENES)) + 
  geom_smooth(method="lm", se=TRUE) + 
  labs(x='-log10(Empirical P-value)', 
       y='Coverage-adjusted score',
       title='T2D: Empirical P-value vs Coverage-adjusted OpenTargets Score')

# Plot by pathway size
p5 <- ggplot(plot_data, aes(x=-log10(P), y=avg_score, color=method)) + 
  geom_point(alpha=0.6) + 
  geom_smooth(method="lm", se=FALSE) + 
  facet_wrap(~ size_bin) +
  labs(x='-log10(Nominal P-value)', y='Open Targets score',
       title='T2D: Nominal P-value vs OpenTargets Score by Pathway Size')

p6 <- ggplot(plot_data, aes(x=-log10(empirical_pvalue), y=avg_score, color=method)) + 
  geom_point(alpha=0.6) + 
  geom_smooth(method="lm", se=FALSE) + 
  facet_wrap(~ size_bin) +
  labs(x='-log10(Empirical P-value)', y='Open Targets score',
       title='T2D: Empirical P-value vs OpenTargets Score by Pathway Size')

# Print all plots
print(p1 / p2)
print(p3 / p4)
print(p5)
print(p6)

# Add scatter plot comparing nominal vs empirical p-values
p7 <- ggplot(plot_data, aes(x=-log10(P), y=-log10(empirical_pvalue), color=size_bin)) + 
  geom_point(alpha=0.8) +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  labs(x='-log10(Nominal P-value)', 
       y='-log10(Empirical P-value)',
       title='T2D: Comparison of Nominal vs Empirical P-values') +
  theme_minimal() +
  facet_wrap(~method)

print(p7)

# Show top pathways based on empirical p-values
top_pathways <- plot_data %>%
  filter(method == "birewire") %>%
  arrange(empirical_pvalue) %>%
  head(20) %>%
  select(name, empirical_pvalue, P, avg_score, coverage, num_genes, num_with_evidence)

# Plot highlighting top pathways
p8 <- ggplot(plot_data %>% filter(method == "birewire"), 
             aes(x=-log10(empirical_pvalue), y=avg_score)) +
  geom_point(alpha=0.3, color="grey") +
  geom_point(data=plot_data %>% filter(method == "birewire") %>% 
               arrange(empirical_pvalue) %>% head(20),
             aes(size=coverage), color="red") +
  geom_text(data=plot_data %>% filter(method == "birewire") %>% 
              arrange(empirical_pvalue) %>% head(10),
            aes(label=name), size=3, hjust=-0.1, vjust=0) +
  labs(title="Top T2D pathways by empirical p-value",
       x="-log10(Empirical P-value)", y="OpenTargets Score") +
  theme_minimal() +
  theme(legend.position="bottom")

print(p8)

dev.off()

# Write results to CSV
write.csv(plot_data, "t2d_opentargets_with_empirical_pvalues.csv", row.names = FALSE)
write.csv(size_stats, "t2d_size_stratified_correlations.csv", row.names = FALSE)
write.csv(top_pathways, "t2d_top_pathways.csv", row.names = FALSE)

cat("\nAnalysis complete. Results saved to t2d_opentargets_analysis.pdf\n")