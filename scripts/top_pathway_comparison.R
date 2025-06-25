# ===== Top Pathway Prioritization Analysis =====
# This script compares the biological relevance of top-ranked pathways 
# identified by different randomization methods

library(tidyverse)
library(ggplot2)
library(patchwork)
library(data.table)
library(jsonlite)
library(GSA)

# Set parameters
disease_name <- "CAD"
disease_id <- "EFO_0001645"  # CAD
sig_threshold <- 0.05

# Set working directory - adjust as needed
setwd('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/results/drugtarget_test/associationByDatatypeDirect')

# ===== Load data =====
cat("Loading pathway and OpenTargets data...\n")

# Read gene set data
dat <- GSA.read.gmt('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/pathway_db/msigdb/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt')
path.list <- dat$genesets
names(path.list) <- dat$geneset.names

# Convert the path list to a data table
genes_long <- rbindlist(lapply(names(path.list), function(name) {
  data.table(value = path.list[[name]], name = name)
}))
genes_long$targetId <- genes_long$value

# Process OpenTargets JSON files for disease
cat("Processing OpenTargets data...\n")
files <- list.files(pattern="*.json")

# Simple function to read JSON files
process_disease_data <- function(files, disease_id) {
  cat("Reading JSON files...\n")
  
  dat <- lapply(files, function(file) {
    tryCatch({
      json_data <- stream_in(file(file), verbose = FALSE)
      json_data <- as.data.frame(json_data)
      json_data %>%
        filter(diseaseId == disease_id & datatypeId %in% c('animal_model', 'known_drug', 'literature'))
    }, error = function(e) {
      cat("Error reading file:", file, "\n")
      return(NULL)
    })
  })
  
  combined_dat <- do.call(rbind, dat[!sapply(dat, is.null)])
  
  if(nrow(combined_dat) == 0) {
    stop(paste("No OpenTargets data found for disease ID:", disease_id))
  }
  
  return(combined_dat)
}

# Read disease OpenTargets data
disease_data <- process_disease_data(files, disease_id)

# Select max evidence count for each target
disease_targets <- disease_data %>%
  group_by(targetId) %>%
  slice_max(evidenceCount, with_ties=FALSE) %>%
  ungroup()

# Read MAGMA results with different randomization methods
cat("Reading MAGMA results...\n")
real_results <- paste0('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_real/cad/cad_real_set.gsa.out')

# Read empirical p-values file - both methods are in the same file
emp_p_file <- '/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_random/birewire/msigdbgenes/cad/cad_empirical_pvalues.csv'
cat("Reading empirical p-values from:", emp_p_file, "\n")
emp_p_data <- read.csv(emp_p_file)

# Extract results for each method
birewire_emp_p <- emp_p_data %>% filter(method == "birewire")
keeppath_emp_p <- emp_p_data %>% filter(method == "keeppathsize")

# Check if we have data for both methods
if(nrow(birewire_emp_p) == 0) {
  stop("No birewire results found in empirical p-values file")
}
if(nrow(keeppath_emp_p) == 0) {
  stop("No keeppathsize results found in empirical p-values file")
}

cat("Found", nrow(birewire_emp_p), "BireWire results and", 
    nrow(keeppath_emp_p), "KeepPathSize results\n")

# Read real MAGMA results
magma_results <- read.table(real_results, header = TRUE)
magma_results$name <- magma_results$FULL_NAME

# Merge with empirical p-values
birewire_results <- merge(
  magma_results, 
  birewire_emp_p %>% select(-method),
  by.x = "FULL_NAME", by.y = "FULL_NAME"
)

keeppath_results <- merge(
  magma_results, 
  keeppath_emp_p %>% select(-method),
  by.x = "FULL_NAME", by.y = "FULL_NAME"  
)

# Merge OpenTargets data with pathway genes
masterfin <- merge(genes_long, disease_targets, by="targetId", all.x=TRUE)
masterfin[is.na(masterfin)] <- 0

# Calculate pathway-level scores
pathway_scores <- masterfin %>% 
  group_by(name) %>% 
  summarise(
    avg_score = mean(score, na.rm = TRUE),
    num_genes = n(),
    num_with_evidence = sum(score > 0, na.rm = TRUE),
    coverage = sum(score > 0, na.rm = TRUE) / n(),
    evidence_density = num_with_evidence / num_genes,
    max_score = max(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE)
  ) %>% 
  as.data.frame()

# ===== Top Pathway Analysis =====
cat("\n===== Top Pathway Analysis =====\n")

# Function to analyze top N pathways
analyze_top_pathways <- function(birewire_data, keeppath_data, pathway_scores, N_values=c(10, 20, 50, 100)) {
  results <- list()
  
  for(N in N_values) {
    cat(paste("\nAnalyzing top", N, "pathways...\n"))
    
    # Get top N pathways by each method
    top_birewire <- birewire_data %>% 
      arrange(empirical_pvalue) %>% 
      head(N)
      
    top_keeppath <- keeppath_data %>% 
      arrange(empirical_pvalue) %>% 
      head(N)
    
    # Pathway overlap
    common_pathways <- intersect(top_birewire$name, top_keeppath$name)
    birewire_only <- setdiff(top_birewire$name, top_keeppath$name)
    keeppath_only <- setdiff(top_keeppath$name, top_birewire$name)
    
    cat("Pathway overlap:\n")
    cat("- Common pathways:", length(common_pathways), "\n")
    cat("- Unique to BireWire:", length(birewire_only), "\n")
    cat("- Unique to KeepPathSize:", length(keeppath_only), "\n")
    
    # Calculate evidence metrics for each set
    get_evidence_metrics <- function(pathway_names) {
      pathway_data <- pathway_scores %>% filter(name %in% pathway_names)
      if(nrow(pathway_data) == 0) {
        return(list(
          n = 0,
          mean_size = NA,
          mean_ot_score = NA,
          median_ot_score = NA,
          mean_evidence_density = NA,
          median_evidence_density = NA,
          pct_with_evidence = NA,
          mean_coverage = NA,
          mean_max_score = NA
        ))
      }
      return(list(
        n = nrow(pathway_data),
        mean_size = mean(pathway_data$num_genes),
        mean_ot_score = mean(pathway_data$avg_score),
        median_ot_score = median(pathway_data$avg_score),
        mean_evidence_density = mean(pathway_data$evidence_density),
        median_evidence_density = median(pathway_data$evidence_density),
        pct_with_evidence = 100 * mean(pathway_data$num_with_evidence > 0),
        mean_coverage = mean(pathway_data$coverage),
        mean_max_score = mean(pathway_data$max_score)
      ))
    }
    
    # Get metrics for each group
    birewire_metrics <- get_evidence_metrics(top_birewire$name)
    keeppath_metrics <- get_evidence_metrics(top_keeppath$name)
    common_metrics <- get_evidence_metrics(common_pathways)
    birewire_only_metrics <- get_evidence_metrics(birewire_only)
    keeppath_only_metrics <- get_evidence_metrics(keeppath_only)
    
    # Create result entry
    results[[as.character(N)]] <- list(
      N = N,
      common_count = length(common_pathways),
      birewire_only_count = length(birewire_only),
      keeppath_only_count = length(keeppath_only),
      
      # Size metrics
      birewire_mean_size = birewire_metrics$mean_size,
      keeppath_mean_size = keeppath_metrics$mean_size,
      birewire_only_mean_size = birewire_only_metrics$mean_size,
      keeppath_only_mean_size = keeppath_only_metrics$mean_size,
      
      # Evidence metrics
      birewire_mean_ot_score = birewire_metrics$mean_ot_score,
      keeppath_mean_ot_score = keeppath_metrics$mean_ot_score,
      birewire_only_mean_ot_score = birewire_only_metrics$mean_ot_score,
      keeppath_only_mean_ot_score = keeppath_only_metrics$mean_ot_score,
      
      # Evidence density
      birewire_mean_evidence_density = birewire_metrics$mean_evidence_density,
      keeppath_mean_evidence_density = keeppath_metrics$mean_evidence_density,
      birewire_only_mean_evidence_density = birewire_only_metrics$mean_evidence_density,
      keeppath_only_mean_evidence_density = keeppath_only_metrics$mean_evidence_density,
      
      # Max score metrics (strongest evidence in pathway)
      birewire_mean_max_score = birewire_metrics$mean_max_score,
      keeppath_mean_max_score = keeppath_metrics$mean_max_score,
      birewire_only_mean_max_score = birewire_only_metrics$mean_max_score,
      keeppath_only_mean_max_score = keeppath_only_metrics$mean_max_score,
      
      # Pathway names for further analysis
      common_pathways = common_pathways,
      birewire_only = birewire_only,
      keeppath_only = keeppath_only
    )
    
    # Print evidence comparison
    cat("\nEvidence comparison:\n")
    cat("BireWire top", N, "pathways:\n")
    cat("- Mean OpenTargets score:", round(birewire_metrics$mean_ot_score, 4), "\n")
    cat("- Mean evidence density:", round(birewire_metrics$mean_evidence_density, 4), "\n")
    
    cat("KeepPathSize top", N, "pathways:\n")
    cat("- Mean OpenTargets score:", round(keeppath_metrics$mean_ot_score, 4), "\n")
    cat("- Mean evidence density:", round(keeppath_metrics$mean_evidence_density, 4), "\n")
    
    if(length(birewire_only) > 0) {
      cat("\nUnique to BireWire:\n")
      cat("- Mean OpenTargets score:", round(birewire_only_metrics$mean_ot_score, 4), "\n")
      cat("- Mean evidence density:", round(birewire_only_metrics$mean_evidence_density, 4), "\n")
    }
    
    if(length(keeppath_only) > 0) {
      cat("Unique to KeepPathSize:\n")
      cat("- Mean OpenTargets score:", round(keeppath_only_metrics$mean_ot_score, 4), "\n")
      cat("- Mean evidence density:", round(keeppath_only_metrics$mean_evidence_density, 4), "\n")
    }
  }
  
  return(results)
}

# Analyze top pathways for different N values
top_pathway_results <- analyze_top_pathways(
  birewire_results, 
  keeppath_results, 
  pathway_scores, 
  N_values=c(10, 20, 50, 100)
)

# Convert results to data frame for plotting and saving
results_df <- do.call(rbind, lapply(names(top_pathway_results), function(N) {
  result <- top_pathway_results[[N]]
  data.frame(
    N = result$N,
    common_count = result$common_count,
    birewire_only_count = result$birewire_only_count,
    keeppath_only_count = result$keeppath_only_count,
    
    birewire_mean_size = result$birewire_mean_size,
    keeppath_mean_size = result$keeppath_mean_size,
    
    birewire_mean_ot_score = result$birewire_mean_ot_score,
    keeppath_mean_ot_score = result$keeppath_mean_ot_score,
    
    birewire_mean_evidence_density = result$birewire_mean_evidence_density,
    keeppath_mean_evidence_density = result$keeppath_mean_evidence_density,
    
    birewire_only_mean_ot_score = result$birewire_only_mean_ot_score,
    keeppath_only_mean_ot_score = result$keeppath_only_mean_ot_score,
    
    birewire_only_mean_evidence_density = result$birewire_only_mean_evidence_density,
    keeppath_only_mean_evidence_density = result$keeppath_only_mean_evidence_density,
    
    birewire_mean_max_score = result$birewire_mean_max_score,
    keeppath_mean_max_score = result$keeppath_mean_max_score
  )
}))

# ===== Create Visualizations =====
cat("\nCreating visualizations...\n")

# Prepare long-format data for cleaner bar plots
evidence_long <- results_df %>%
  pivot_longer(cols=c(birewire_mean_ot_score, keeppath_mean_ot_score),
               names_to="method", 
               values_to="mean_ot_score") %>%
  mutate(method = ifelse(method == "birewire_mean_ot_score", "BireWire", "KeepPathSize"))

density_long <- results_df %>%
  pivot_longer(cols=c(birewire_mean_evidence_density, keeppath_mean_evidence_density),
               names_to="method", 
               values_to="mean_evidence_density") %>%
  mutate(method = ifelse(method == "birewire_mean_evidence_density", "BireWire", "KeepPathSize"))

max_score_long <- results_df %>%
  pivot_longer(cols=c(birewire_mean_max_score, keeppath_mean_max_score),
               names_to="method", 
               values_to="mean_max_score") %>%
  mutate(method = ifelse(method == "birewire_mean_max_score", "BireWire", "KeepPathSize"))

# Create better bar plots
p1 <- ggplot(evidence_long, aes(x=as.factor(N), y=mean_ot_score, fill=method)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=c("BireWire"="blue", "KeepPathSize"="red")) +
  labs(title=paste(disease_name, "- Mean OpenTargets Evidence in Top Pathways"),
       x="Top N Pathways", y="Mean OpenTargets Score", fill="Method") +
  theme_minimal()

p2 <- ggplot(density_long, aes(x=as.factor(N), y=mean_evidence_density, fill=method)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=c("BireWire"="blue", "KeepPathSize"="red")) +
  labs(title=paste(disease_name, "- Evidence Density in Top Pathways"),
       x="Top N Pathways", y="Mean Evidence Density", fill="Method") +
  theme_minimal()

p3 <- ggplot(max_score_long, aes(x=as.factor(N), y=mean_max_score, fill=method)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=c("BireWire"="blue", "KeepPathSize"="red")) +
  labs(title=paste(disease_name, "- Maximum Evidence in Top Pathways"),
       x="Top N Pathways", y="Mean Maximum OpenTargets Score", fill="Method") +
  theme_minimal()

# Analysis of unique pathways
unique_evidence_data <- results_df %>%
  pivot_longer(cols=c(birewire_only_mean_ot_score, keeppath_only_mean_ot_score),
               names_to="method", 
               values_to="unique_mean_score") %>%
  mutate(method = ifelse(method == "birewire_only_mean_ot_score", "BireWire Unique", "KeepPathSize Unique"))

unique_density_data <- results_df %>%
  pivot_longer(cols=c(birewire_only_mean_evidence_density, keeppath_only_mean_evidence_density),
               names_to="method", 
               values_to="unique_density") %>%
  mutate(method = ifelse(method == "birewire_only_mean_evidence_density", "BireWire Unique", "KeepPathSize Unique"))

p4 <- ggplot(unique_evidence_data, aes(x=as.factor(N), y=unique_mean_score, fill=method)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=c("BireWire Unique"="darkblue", "KeepPathSize Unique"="darkred")) +
  labs(title=paste(disease_name, "- Evidence in Method-Specific Pathways"),
       x="Top N Pathways", y="Mean OpenTargets Score", fill="Method") +
  theme_minimal()

p5 <- ggplot(unique_density_data, aes(x=as.factor(N), y=unique_density, fill=method)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=c("BireWire Unique"="darkblue", "KeepPathSize Unique"="darkred")) +
  labs(title=paste(disease_name, "- Evidence Density in Method-Specific Pathways"),
       x="Top N Pathways", y="Mean Evidence Density", fill="Method") +
  theme_minimal()

# Summary plot showing which method performs better
p6 <- ggplot(results_df) +
  geom_segment(aes(x=as.factor(N), xend=as.factor(N), 
                   y=0, yend=birewire_mean_evidence_density - keeppath_mean_evidence_density),
               color=ifelse(results_df$birewire_mean_evidence_density > results_df$keeppath_mean_evidence_density, 
                           "blue", "red"),
               size=1.5) +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(title=paste(disease_name, "- Evidence Advantage by Method"),
       subtitle="Positive values favor BireWire, negative values favor KeepPathSize",
       x="Top N Pathways", y="Evidence Density Difference (BireWire - KeepPathSize)") +
  theme_minimal()

# Create overlap analysis
overlap_data <- data.frame(
  N = results_df$N,
  Common = results_df$common_count,
  BireWire_Only = results_df$birewire_only_count,
  KeepPathSize_Only = results_df$keeppath_only_count
) %>%
  pivot_longer(cols=c(Common, BireWire_Only, KeepPathSize_Only),
               names_to="Category", values_to="Count")

p7 <- ggplot(overlap_data, aes(x=as.factor(N), y=Count, fill=Category)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("Common"="purple", "BireWire_Only"="blue", "KeepPathSize_Only"="red")) +
  labs(title=paste(disease_name, "- Overlap Between Top Pathways"),
       x="Top N Pathways", y="Number of Pathways", fill="") +
  theme_minimal()

# Save all plots
pdf(paste0(tolower(disease_name), '_top_pathway_comparison.pdf'), width=10, height=12)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
dev.off()

# ===== Individual Pathway Analysis =====
cat("\n===== Individual Pathway Analysis =====\n")

# Examine top 20 pathways from each method
top_n <- 20
top_birewire <- birewire_results %>% 
  arrange(empirical_pvalue) %>% 
  head(top_n) %>%
  select(name, NGENES, P, empirical_pvalue)

top_keeppath <- keeppath_results %>% 
  arrange(empirical_pvalue) %>% 
  head(top_n) %>%
  select(name, NGENES, P, empirical_pvalue)

# Add OpenTargets evidence
top_birewire <- merge(top_birewire, pathway_scores, by="name")
top_keeppath <- merge(top_keeppath, pathway_scores, by="name")

# Create detailed tables
top_birewire <- top_birewire %>%
  arrange(empirical_pvalue) %>%
  mutate(rank = row_number()) %>%
  select(rank, name, NGENES, empirical_pvalue, P, avg_score, evidence_density, num_with_evidence, max_score)

top_keeppath <- top_keeppath %>%
  arrange(empirical_pvalue) %>%
  mutate(rank = row_number()) %>%
  select(rank, name, NGENES, empirical_pvalue, P, avg_score, evidence_density, num_with_evidence, max_score)

# Write detailed pathway tables to CSV
write.csv(top_birewire, paste0(tolower(disease_name), '_top_birewire_pathways.csv'), row.names=FALSE)
write.csv(top_keeppath, paste0(tolower(disease_name), '_top_keeppath_pathways.csv'), row.names=FALSE)
write.csv(results_df, paste0(tolower(disease_name), '_top_pathway_summary.csv'), row.names=FALSE)

# ===== Statistical Tests =====
cat("\n===== Statistical Tests =====\n")

# Perform paired t-tests to compare methods
# Test 1: Compare mean OT scores in top pathways
t_test_scores <- t.test(
  results_df$birewire_mean_ot_score,
  results_df$keeppath_mean_ot_score,
  paired = TRUE
)

# Test 2: Compare evidence density in top pathways
t_test_density <- t.test(
  results_df$birewire_mean_evidence_density,
  results_df$keeppath_mean_evidence_density,
  paired = TRUE
)

# Test 3: Compare max scores in top pathways
t_test_max <- t.test(
  results_df$birewire_mean_max_score,
  results_df$keeppath_mean_max_score,
  paired = TRUE
)

# Print test results
cat("\nPaired t-test comparing mean OpenTargets scores:\n")
cat("t =", t_test_scores$statistic, ", p-value =", t_test_scores$p.value, "\n")
cat("BireWire mean:", mean(results_df$birewire_mean_ot_score), "\n")
cat("KeepPathSize mean:", mean(results_df$keeppath_mean_ot_score), "\n")

cat("\nPaired t-test comparing evidence density:\n")
cat("t =", t_test_density$statistic, ", p-value =", t_test_density$p.value, "\n")
cat("BireWire mean:", mean(results_df$birewire_mean_evidence_density), "\n")
cat("KeepPathSize mean:", mean(results_df$keeppath_mean_evidence_density), "\n")

cat("\nPaired t-test comparing maximum scores:\n")
cat("t =", t_test_max$statistic, ", p-value =", t_test_max$p.value, "\n")
cat("BireWire mean:", mean(results_df$birewire_mean_max_score), "\n")
cat("KeepPathSize mean:", mean(results_df$keeppath_mean_max_score), "\n")

# ===== Final Summary =====
cat("\n===== Final Summary =====\n")
better_method <- ifelse(mean(results_df$birewire_mean_evidence_density) > mean(results_df$keeppath_mean_evidence_density),
                      "BireWire", "KeepPathSize")

cat("Overall better method for top pathway prioritization:", better_method, "\n")
cat("Evidence density advantage:", 
    abs(mean(results_df$birewire_mean_evidence_density) - mean(results_df$keeppath_mean_evidence_density)), "\n")
cat("Analysis complete. Results saved to CSV files and PDF.\n")

# ===== Pathway Size Analysis =====
cat("\n===== Pathway Size Analysis =====\n")

# Calculate average sizes for each method's top pathways
size_analysis <- results_df %>%
  select(N, birewire_mean_size, keeppath_mean_size) %>%
  mutate(
    size_difference = birewire_mean_size - keeppath_mean_size,
    size_ratio = birewire_mean_size / keeppath_mean_size,
    percent_difference = 100 * (birewire_mean_size - keeppath_mean_size) / keeppath_mean_size
  )

# Print size comparison
cat("Average pathway sizes:\n")
for(i in 1:nrow(size_analysis)) {
  cat(paste0("Top ", size_analysis$N[i], " pathways:\n"))
  cat("- BireWire: ", round(size_analysis$birewire_mean_size[i], 1), " genes\n", sep="")
  cat("- KeepPathSize: ", round(size_analysis$keeppath_mean_size[i], 1), " genes\n", sep="")
  cat("- Difference: ", round(size_analysis$size_difference[i], 1), " genes (", 
      ifelse(size_analysis$size_difference[i] > 0, "+", ""), 
      round(size_analysis$percent_difference[i], 1), "%)\n\n", sep="")
}

# Size distribution of top pathways
# Extract the actual sizes for the top N pathways from each method
gather_top_n_sizes <- function(results, pathway_scores, N) {
  top_birewire <- results %>% 
    arrange(empirical_pvalue) %>% 
    head(N) %>%
    select(name, NGENES)
  
  top_keeppath <- keeppath_results %>% 
    arrange(empirical_pvalue) %>% 
    head(N) %>%
    select(name, NGENES)
  
  birewire_sizes <- data.frame(
    method = "BireWire",
    N = N,
    name = top_birewire$name,
    size = top_birewire$NGENES
  )
  
  keeppath_sizes <- data.frame(
    method = "KeepPathSize",
    N = N,
    name = top_keeppath$name,
    size = top_keeppath$NGENES
  )
  
  return(rbind(birewire_sizes, keeppath_sizes))
}

# Gather size data for different N values
all_sizes <- do.call(rbind, lapply(c(10, 20, 50, 100), function(n) {
  gather_top_n_sizes(birewire_results, keeppath_results, n)
}))

# Statistical test for size differences
size_tests <- data.frame()

for(n_value in c(10, 20, 50, 100)) {
  # Extract sizes for each method at this N value
  birewire_sizes <- all_sizes$size[all_sizes$method == "BireWire" & all_sizes$N == n_value]
  keeppath_sizes <- all_sizes$size[all_sizes$method == "KeepPathSize" & all_sizes$N == n_value]
  
  # Perform Wilcoxon rank sum test (Mann-Whitney U test)
  wilcox_result <- wilcox.test(birewire_sizes, keeppath_sizes)
  
  size_tests <- rbind(size_tests, data.frame(
    N = n_value,
    birewire_median = median(birewire_sizes),
    keeppath_median = median(keeppath_sizes),
    p_value = wilcox_result$p.value
  ))
}

cat("\nStatistical tests for pathway size differences:\n")
print(size_tests)

# Create visualizations
# Boxplot of pathway sizes
p8 <- ggplot(all_sizes, aes(x=as.factor(N), y=size, fill=method)) +
  geom_boxplot() +
  scale_fill_manual(values=c("BireWire"="blue", "KeepPathSize"="red")) +
  labs(title=paste(disease_name, "- Size Distribution of Top Pathways"),
       x="Top N Pathways", y="Pathway Size (# of Genes)", fill="Method") +
  theme_minimal()

# Violin plot for more detailed size distribution
p9 <- ggplot(all_sizes, aes(x=as.factor(N), y=size, fill=method)) +
  geom_violin(position=position_dodge(width=0.8), alpha=0.7) +
  geom_jitter(aes(color=method), position=position_dodge(width=0.8), alpha=0.5, size=1) +
  scale_fill_manual(values=c("BireWire"="blue", "KeepPathSize"="red")) +
  scale_color_manual(values=c("BireWire"="darkblue", "KeepPathSize"="darkred")) +
  labs(title=paste(disease_name, "- Detailed Size Distribution of Top Pathways"),
       x="Top N Pathways", y="Pathway Size (# of Genes)", fill="Method") +
  theme_minimal() +
  guides(color="none")

# Add to existing PDF
pdf(paste0(tolower(disease_name), '_top_pathway_comparison.pdf'), width=10, height=12)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
print(p8)
print(p9)

# Additional plot: Size difference in relation to evidence advantage
p10 <- ggplot(size_analysis, aes(x=size_difference, 
                                y=results_df$birewire_mean_evidence_density - results_df$keeppath_mean_evidence_density)) +
  geom_point(aes(size=N), color="purple") +
  geom_text(aes(label=N), vjust=-0.7) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  labs(title=paste(disease_name, "- Relationship Between Size Difference and Evidence Advantage"),
       x="Size Difference (BireWire - KeepPathSize)",
       y="Evidence Density Advantage (BireWire - KeepPathSize)") +
  theme_minimal()
print(p10)
dev.off()

# Create a size-specific CSV file
write.csv(size_analysis, paste0(tolower(disease_name), '_pathway_size_analysis.csv'), row.names=FALSE)
write.csv(all_sizes, paste0(tolower(disease_name), '_detailed_pathway_sizes.csv'), row.names=FALSE)
write.csv(size_tests, paste0(tolower(disease_name), '_size_statistical_tests.csv'), row.names=FALSE)

# Add a section comparing top pathways by size category
# Define size categories
all_sizes$size_category <- cut(all_sizes$size, 
                              breaks=c(0, 50, 100, 200, 500, Inf),
                              labels=c("<50", "50-100", "100-200", "200-500", ">500"))

# Count pathways in each size bin
size_distribution <- all_sizes %>%
  group_by(method, N, size_category) %>%
  summarize(count = n()) %>%
  spread(size_category, count, fill=0)

# Print distribution
cat("\nSize distribution of top pathways:\n")
print(size_distribution)

# Create stacked bar plot of size distribution
size_dist_long <- all_sizes %>%
  group_by(method, N, size_category) %>%
  summarize(count = n()) %>%
  ungroup()

p11 <- ggplot(size_dist_long, aes(x=as.factor(N), y=count, fill=size_category)) +
  geom_bar(stat="identity") +
  facet_wrap(~method) +
  scale_fill_brewer(palette="Set3") +
  labs(title=paste(disease_name, "- Size Distribution of Top Pathways"),
       x="Top N Pathways", y="Number of Pathways", fill="Size Category") +
  theme_minimal()

# Add to the PDF
pdf(paste0(tolower(disease_name), '_pathway_size_distribution.pdf'), width=10, height=6)
print(p11)
dev.off()

# Add to the summary results
write.csv(size_dist_long, paste0(tolower(disease_name), '_pathway_size_distribution.csv'), row.names=FALSE)

cat("\nPathway size analysis complete. Results added to PDFs and CSV files.\n")

# ===== Size-Matched Pathway Analysis =====
cat("\n===== Size-Matched Pathway Analysis =====\n")

# Function to perform analysis on size-matched pathways
analyze_size_matched_pathways <- function(birewire_results, keeppath_results, pathway_scores, 
                                         size_bins=c(0, 50, 100, 200, 500, Inf)) {
  
  # Create size categories
  birewire_results$size_bin <- cut(birewire_results$NGENES, breaks=size_bins)
  keeppath_results$size_bin <- cut(keeppath_results$NGENES, breaks=size_bins)
  
  # Results container
  all_results <- data.frame()
  
  # Analyze each size bin separately
  for(bin in levels(birewire_results$size_bin)) {
    cat("\nAnalyzing size bin:", bin, "\n")
    
    # Get pathways in this size bin
    bin_birewire <- birewire_results[birewire_results$size_bin == bin, ]
    bin_keeppath <- keeppath_results[keeppath_results$size_bin == bin, ]
    
    # Get top N for each method within this size bin
    n_values <- c(5, 10, 20)
    n_values <- n_values[n_values <= min(nrow(bin_birewire), nrow(bin_keeppath))/2]
    
    if(length(n_values) == 0) {
      cat("Not enough pathways in bin", bin, "for analysis\n")
      next
    }
    
    for(n in n_values) {
      # Get top n pathways in this size bin
      top_birewire <- bin_birewire %>% arrange(empirical_pvalue) %>% head(n)
      top_keeppath <- bin_keeppath %>% arrange(empirical_pvalue) %>% head(n)
      
      # Get evidence metrics
      bw_evidence <- pathway_scores %>% 
        filter(name %in% top_birewire$name) %>%
        summarize(
          mean_size = mean(num_genes),
          mean_score = mean(avg_score),
          evidence_density = mean(evidence_density),
          max_score = mean(max_score)
        )
      
      kp_evidence <- pathway_scores %>% 
        filter(name %in% top_keeppath$name) %>%
        summarize(
          mean_size = mean(num_genes),
          mean_score = mean(avg_score),
          evidence_density = mean(evidence_density),
          max_score = mean(max_score)
        )
      
      # Store results
      bin_result <- data.frame(
        size_bin = bin,
        n_pathways = n,
        mean_size_birewire = bw_evidence$mean_size,
        mean_size_keeppath = kp_evidence$mean_size,
        
        mean_score_birewire = bw_evidence$mean_score,
        mean_score_keeppath = kp_evidence$mean_score,
        score_advantage = bw_evidence$mean_score - kp_evidence$mean_score,
        
        density_birewire = bw_evidence$evidence_density, 
        density_keeppath = kp_evidence$evidence_density,
        density_advantage = bw_evidence$evidence_density - kp_evidence$evidence_density,
        
        max_score_birewire = bw_evidence$max_score,
        max_score_keeppath = kp_evidence$max_score,
        max_advantage = bw_evidence$max_score - kp_evidence$max_score
      )
      
      all_results <- rbind(all_results, bin_result)
      
      # Print results
      cat("Top", n, "pathways in size bin", bin, ":\n")
      cat("- BireWire mean score:", round(bw_evidence$mean_score, 4), 
          ", evidence density:", round(bw_evidence$evidence_density, 4), "\n")
      cat("- KeepPathSize mean score:", round(kp_evidence$mean_score, 4), 
          ", evidence density:", round(kp_evidence$evidence_density, 4), "\n")
    }
  }
  
  return(all_results)
}

# Run size-matched analysis
size_matched_results <- analyze_size_matched_pathways(birewire_results, keeppath_results, pathway_scores)

# Visualize results
if(nrow(size_matched_results) > 0) {
  # Plot advantage by size bin - USING ORIGINAL CODE WITHOUT FIXING X-AXIS
  size_matched_long <- size_matched_results %>%
    select(size_bin, n_pathways, score_advantage, density_advantage, max_advantage) %>%
    pivot_longer(cols=c(score_advantage, density_advantage, max_advantage),
                 names_to="metric", values_to="advantage")
  
  p_size_matched <- ggplot(size_matched_long, 
                          aes(x=size_bin, y=advantage, fill=metric)) +
    geom_bar(stat="identity", position=position_dodge()) +
    facet_wrap(~n_pathways, labeller=labeller(n_pathways=function(x) paste0("Top ", x, " Pathways"))) +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_fill_brewer(palette="Set1", 
                    labels=c("density_advantage"="Evidence Density", 
                             "max_advantage"="Max Evidence", 
                             "score_advantage"="Mean Score")) +
    labs(title="BireWire Advantage by Size Bin (Controlling for Pathway Size)",
         x="Pathway Size Bin", y="BireWire Advantage", fill="Metric") +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  
  # Save results
  pdf(paste0(tolower(disease_name), '_size_matched_analysis.pdf'), width=10, height=7)
  print(p_size_matched)
  dev.off()
  
  # Save to CSV
  write.csv(size_matched_results, paste0(tolower(disease_name), '_size_matched_results.csv'), row.names=FALSE)
}













