# ===== Top Pathway Prioritization Analysis =====
# This script compares the biological relevance of top-ranked pathways 
# identified by different randomization methods

library(tidyverse)
library(ggplot2)
library(patchwork)
library(data.table)
library(jsonlite)
library(GSA)
library(MatchIt)

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
    mean_score = mean(score, na.rm = TRUE),
    num_genes = n(),
    num_with_evidence = sum(score > 0, na.rm = TRUE),
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
        mean_ot_score = mean(pathway_data$mean_score),
        median_ot_score = median(pathway_data$mean_score),
        mean_evidence_density = mean(pathway_data$evidence_density),
        median_evidence_density = median(pathway_data$evidence_density),
        pct_with_evidence = 100 * mean(pathway_data$num_with_evidence > 0),
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
write.csv(top_birewire, paste0(tolower(disease_name), '_top_birewire_pathways.csv'), row.names=FALSE,quote=F)
write.csv(top_keeppath, paste0(tolower(disease_name), '_top_keeppath_pathways.csv'), row.names=FALSE,quote=F)
write.csv(results_df, paste0(tolower(disease_name), '_top_pathway_summary.csv'), row.names=FALSE,quote=F)


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

# Add to the summary results
write.csv(size_distribution, paste0(tolower(disease_name), '_pathway_size_distribution.csv'), row.names=FALSE)




# ===== Propensity Score Matching Analysis =====
cat("\n===== Propensity Score Matching Analysis =====\n")
# Function to perform propensity score matching
perform_psm_analysis <- function(birewire_results, keeppath_results, pathway_scores) {
  # Prepare data for matching
  # First, identify all unique pathways from both methods
  all_pathways <- unique(c(birewire_results$name, keeppath_results$name))
  
  # Create a dataset with method assignment and pathway size
  matching_data <- data.frame(
    name = all_pathways,
    in_birewire_top100 = all_pathways %in% head(birewire_results[order(birewire_results$empirical_pvalue),]$name, 20),
    in_keeppath_top100 = all_pathways %in% head(keeppath_results[order(keeppath_results$empirical_pvalue),]$name, 20)
  )
  
  # Add pathway size information
  matching_data$size <- NA
  for(i in 1:nrow(matching_data)) {
    idx <- which(birewire_results$name == matching_data$name[i])
    if(length(idx) > 0) {
      matching_data$size[i] <- birewire_results$NGENES[idx[1]]
    } else {
      idx <- which(keeppath_results$name == matching_data$name[i])
      if(length(idx) > 0) {
        matching_data$size[i] <- keeppath_results$NGENES[idx[1]]
      }
    }
  }
  
  # Remove pathways with missing size
  matching_data <- matching_data[!is.na(matching_data$size),]
  
  # Add OpenTargets evidence information
  matching_data <- merge(matching_data, pathway_scores[,c("name", "mean_score")], by="name", all.x=TRUE)
  
  # Perform size matching for BireWire top pathways
  tryCatch({
    m.out1 <- matchit(in_birewire_top100 ~ size, data = matching_data, 
                     method = "nearest", ratio = 1)
    
    matched_data1 <- match.data(m.out1)
    
    # Compare OT scores
    bw_density <- mean(matched_data1$mean_score[matched_data1$in_birewire_top100], na.rm=TRUE)
    other_density <- mean(matched_data1$mean_score[!matched_data1$in_birewire_top100], na.rm=TRUE)
    
    cat("Size-matched comparison (BireWire top 100 vs. size-matched pathways):\n")
    cat("- Average evidence density in BireWire top 100:", round(bw_density, 4), "\n")
    cat("- Average evidence density in size-matched pathways:", round(other_density, 4), "\n")
    cat("- Advantage:", round(bw_density - other_density, 4), 
        "(", ifelse(bw_density > other_density, "BireWire better", "Other better"), ")\n\n")
    
    # Perform matching for KeepPathSize top pathways
    m.out2 <- matchit(in_keeppath_top100 ~ size, data = matching_data, 
                     method = "nearest", ratio = 1)
    
    matched_data2 <- match.data(m.out2)
    
    # Compare OT scores
    kp_density <- mean(matched_data2$mean_score[matched_data2$in_keeppath_top100], na.rm=TRUE)
    other_density2 <- mean(matched_data2$mean_score[!matched_data2$in_keeppath_top100], na.rm=TRUE)
    
    cat("Size-matched comparison (KeepPathSize top 100 vs. size-matched pathways):\n")
    cat("- Average OT score in KeepPathSize top 100:", round(kp_density, 4), "\n")
    cat("- Average OT score in size-matched pathways:", round(other_density2, 4), "\n")
    cat("- Advantage:", round(kp_density - other_density2, 4),
        "(", ifelse(kp_density > other_density2, "KeepPathSize better", "Other better"), ")\n")
    
    return(list(
      birewire_match = matched_data1,
      keeppath_match = matched_data2
    ))
  }, error = function(e) {
    cat("Error in pathway size matching:", e$message, "\n")
    cat("Consider installing the MatchIt package: install.packages('MatchIt')\n")
    return(NULL)
  })
}

# Use the function
psm_results <- perform_psm_analysis(birewire_results, keeppath_results, pathway_scores)

# If successful, save results
if(!is.null(psm_results)) {
  lapply(names(psm_results), function(name) {
    write.csv(psm_results[[name]], 
             paste0(tolower(disease_name), '_psm_', name, '.csv'), 
             row.names=FALSE,quote=F)
  })
}


# ===== Enhanced Size Matching Analysis =====
cat("\n===== Enhanced Size Matching Analysis =====\n")

# Function to perform size matching with multiple values of N
perform_enhanced_matching <- function(birewire_results, keeppath_results, pathway_scores, n_values = c(10, 20, 50, 100)) {
  all_results <- list()
  
  for (n in n_values) {
    cat(paste("\nAnalyzing top", n, "pathways...\n"))
    
    # Prepare data for matching
    # First, identify all unique pathways from both methods
    all_pathways <- unique(c(birewire_results$name, keeppath_results$name))
    
    # Create a dataset with method assignment and pathway size
    matching_data <- data.frame(
      name = all_pathways,
      in_birewire_top = all_pathways %in% head(birewire_results[order(birewire_results$empirical_pvalue),]$name, n),
      in_keeppath_top = all_pathways %in% head(keeppath_results[order(keeppath_results$empirical_pvalue),]$name, n)
    )
    
    # Add pathway size information
    matching_data$size <- NA
    for(i in 1:nrow(matching_data)) {
      idx <- which(birewire_results$name == matching_data$name[i])
      if(length(idx) > 0) {
        matching_data$size[i] <- birewire_results$NGENES[idx[1]]
      } else {
        idx <- which(keeppath_results$name == matching_data$name[i])
        if(length(idx) > 0) {
          matching_data$size[i] <- keeppath_results$NGENES[idx[1]]
        }
      }
    }
    
    # Remove pathways with missing size
    matching_data <- matching_data[!is.na(matching_data$size),]
    
    # Add OpenTargets evidence information - include all metrics
    matching_data <- merge(matching_data, pathway_scores, by="name", all.x=TRUE)
    
    # BireWire analysis
    tryCatch({
      # Perform size matching for BireWire top pathways
      m.out1 <- matchit(in_birewire_top ~ size, data = matching_data, 
                       method = "nearest", ratio = 1)
      
      matched_data1 <- match.data(m.out1)
      
      # Compare all OT metrics
      metrics <- c("mean_score", "evidence_density", "max_score", "median_score")
      bw_comparison <- data.frame()
      
      for(metric in metrics) {
        if(metric %in% colnames(matched_data1)) {
          bw_value <- mean(matched_data1[[metric]][matched_data1$in_birewire_top], na.rm=TRUE)
          other_value <- mean(matched_data1[[metric]][!matched_data1$in_birewire_top], na.rm=TRUE)
          
          bw_comparison <- rbind(bw_comparison, data.frame(
            metric = metric,
            birewire_value = bw_value,
            control_value = other_value,
            difference = bw_value - other_value,
            percent_diff = 100 * (bw_value - other_value) / ifelse(other_value == 0, 1, other_value)
          ))
        }
      }
      
      # Calculate size balance
      bw_size_avg <- mean(matched_data1$size[matched_data1$in_birewire_top])
      other_size_avg <- mean(matched_data1$size[!matched_data1$in_birewire_top])
      
      cat("Size-matched comparison (BireWire top", n, "vs. size-matched pathways):\n")
      cat("- Average size in BireWire top:", round(bw_size_avg, 1), 
          "vs. matched pathways:", round(other_size_avg, 1), "\n")
      
      # Print metric comparisons
      for(i in 1:nrow(bw_comparison)) {
        metric_name <- bw_comparison$metric[i]
        bw_val <- bw_comparison$birewire_value[i]
        other_val <- bw_comparison$control_value[i]
        diff <- bw_comparison$difference[i]
        
        cat("- Average", metric_name, "in BireWire top", n, ":", round(bw_val, 4), 
            "vs. matched pathways:", round(other_val, 4), 
            "(diff:", round(diff, 4), ",", 
            ifelse(diff > 0, "BireWire better)", "Control better)"), "\n")
      }
      
      # KeepPathSize analysis
      m.out2 <- matchit(in_keeppath_top ~ size, data = matching_data, 
                       method = "nearest", ratio = 1)
      
      matched_data2 <- match.data(m.out2)
      
      # Compare all OT metrics for KeepPathSize
      kp_comparison <- data.frame()
      
      for(metric in metrics) {
        if(metric %in% colnames(matched_data2)) {
          kp_value <- mean(matched_data2[[metric]][matched_data2$in_keeppath_top], na.rm=TRUE)
          other_value <- mean(matched_data2[[metric]][!matched_data2$in_keeppath_top], na.rm=TRUE)
          
          kp_comparison <- rbind(kp_comparison, data.frame(
            metric = metric,
            keeppath_value = kp_value,
            control_value = other_value,
            difference = kp_value - other_value,
            percent_diff = 100 * (kp_value - other_value) / ifelse(other_value == 0, 1, other_value)
          ))
        }
      }
      
      # Calculate size balance for KeepPathSize
      kp_size_avg <- mean(matched_data2$size[matched_data2$in_keeppath_top])
      other_size_avg2 <- mean(matched_data2$size[!matched_data2$in_keeppath_top])
      
      cat("\nSize-matched comparison (KeepPathSize top", n, "vs. size-matched pathways):\n")
      cat("- Average size in KeepPathSize top:", round(kp_size_avg, 1), 
          "vs. matched pathways:", round(other_size_avg2, 1), "\n")
      
      # Print metric comparisons for KeepPathSize
      for(i in 1:nrow(kp_comparison)) {
        metric_name <- kp_comparison$metric[i]
        kp_val <- kp_comparison$keeppath_value[i]
        other_val <- kp_comparison$control_value[i]
        diff <- kp_comparison$difference[i]
        
        cat("- Average", metric_name, "in KeepPathSize top", n, ":", round(kp_val, 4), 
            "vs. matched pathways:", round(other_val, 4), 
            "(diff:", round(diff, 4), ",", 
            ifelse(diff > 0, "KeepPathSize better)", "Control better)"), "\n")
      }
      
      # Add N to the matched datasets for plotting
      matched_data1$N <- n
      matched_data2$N <- n
      
      # Add method identifier for later comparison
      matched_data1$method <- "BireWire"
      matched_data2$method <- "KeepPathSize"
      
      # Store results
      all_results[[paste0("birewire_", n)]] <- matched_data1
      all_results[[paste0("keeppath_", n)]] <- matched_data2
      all_results[[paste0("bw_metrics_", n)]] <- bw_comparison
      all_results[[paste0("kp_metrics_", n)]] <- kp_comparison
      
      # Create boxplots for this N value
      create_comparison_plots(matched_data1, matched_data2, n, disease_name)
      
    }, error = function(e) {
      cat("Error in matching for N =", n, ":", e$message, "\n")
    })
  }
  
  # Combine all matched data for overall comparison
  bw_combined <- do.call(rbind, all_results[grep("^birewire_\\d+$", names(all_results))])
  kp_combined <- do.call(rbind, all_results[grep("^keeppath_\\d+$", names(all_results))])
  
  all_results[["birewire_combined"]] <- bw_combined
  all_results[["keeppath_combined"]] <- kp_combined
  
  # Create combined visualizations
  create_combined_plots(bw_combined, kp_combined, disease_name)
  
  return(all_results)
}

# Helper function to create comparison plots for a specific N value
create_comparison_plots <- function(bw_data, kp_data, n, disease_name) {
  # BireWire plots
  pdf(paste0(tolower(disease_name), "_birewire_top", n, "_matching_plots.pdf"), width=12, height=8)
  
  # Mean score comparison
  p1 <- ggplot(bw_data, aes(x=factor(in_birewire_top, labels=c("Control", "BireWire")), 
                            y=mean_score)) +
    geom_boxplot(fill="lightblue", width=0.5) +
    geom_jitter(width=0.2, alpha=0.6) +
    stat_summary(fun=mean, geom="point", shape=23, size=4, fill="red") +
    labs(title=paste("OpenTargets Mean Score - BireWire Top", n, "vs Size-Matched Controls"),
         x="Group", y="Mean OpenTargets Score") +
    theme_minimal()
  
  print(p1)
  
  # Evidence density comparison
  p2 <- ggplot(bw_data, aes(x=factor(in_birewire_top, labels=c("Control", "BireWire")), 
                           y=evidence_density)) +
    geom_boxplot(fill="lightblue", width=0.5) +
    geom_jitter(width=0.2, alpha=0.6) +
    stat_summary(fun=mean, geom="point", shape=23, size=4, fill="red") +
    labs(title=paste("Evidence Density - BireWire Top", n, "vs Size-Matched Controls"),
         x="Group", y="Evidence Density") +
    theme_minimal()
  
  print(p2)
  
  # Size verification - to confirm matching worked
  p3 <- ggplot(bw_data, aes(x=factor(in_birewire_top, labels=c("Control", "BireWire")), 
                           y=size)) +
    geom_boxplot(fill="lightblue", width=0.5) +
    geom_jitter(width=0.2, alpha=0.6) +
    stat_summary(fun=mean, geom="point", shape=23, size=4, fill="red") +
    labs(title=paste("Pathway Size - BireWire Top", n, "vs Size-Matched Controls"),
         x="Group", y="Pathway Size (# genes)") +
    theme_minimal()
  
  print(p3)
  
  dev.off()
  
  # KeepPathSize plots
  pdf(paste0(tolower(disease_name), "_keeppath_top", n, "_matching_plots.pdf"), width=12, height=8)
  
  # Mean score comparison
  p4 <- ggplot(kp_data, aes(x=factor(in_keeppath_top, labels=c("Control", "KeepPathSize")), 
                            y=mean_score)) +
    geom_boxplot(fill="lightpink", width=0.5) +
    geom_jitter(width=0.2, alpha=0.6) +
    stat_summary(fun=mean, geom="point", shape=23, size=4, fill="red") +
    labs(title=paste("OpenTargets Mean Score - KeepPathSize Top", n, "vs Size-Matched Controls"),
         x="Group", y="Mean OpenTargets Score") +
    theme_minimal()
  
  print(p4)
  
  # Evidence density comparison
  p5 <- ggplot(kp_data, aes(x=factor(in_keeppath_top, labels=c("Control", "KeepPathSize")), 
                           y=evidence_density)) +
    geom_boxplot(fill="lightpink", width=0.5) +
    geom_jitter(width=0.2, alpha=0.6) +
    stat_summary(fun=mean, geom="point", shape=23, size=4, fill="red") +
    labs(title=paste("Evidence Density - KeepPathSize Top", n, "vs Size-Matched Controls"),
         x="Group", y="Evidence Density") +
    theme_minimal()
  
  print(p5)
  
  # Size verification - to confirm matching worked
  p6 <- ggplot(kp_data, aes(x=factor(in_keeppath_top, labels=c("Control", "KeepPathSize")), 
                           y=size)) +
    geom_boxplot(fill="lightpink", width=0.5) +
    geom_jitter(width=0.2, alpha=0.6) +
    stat_summary(fun=mean, geom="point", shape=23, size=4, fill="red") +
    labs(title=paste("Pathway Size - KeepPathSize Top", n, "vs Size-Matched Controls"),
         x="Group", y="Pathway Size (# genes)") +
    theme_minimal()
  
  print(p6)
  
  dev.off()
}

# Helper function to create combined plots across all N values
create_combined_plots <- function(bw_combined, kp_combined, disease_name) {
  pdf(paste0(tolower(disease_name), "_combined_matching_analysis.pdf"), width=14, height=10)
  
  # Reshape BireWire data for faceted plot
  bw_plot_data <- bw_combined %>%
    mutate(Group = factor(in_birewire_top, 
                        labels = c("Size-Matched Control", "BireWire Top Pathways"))) %>%
    select(N, Group, name, size, mean_score, evidence_density, max_score) %>%
    mutate(N = factor(N))
  
  # Reshape KeepPathSize data for faceted plot
  kp_plot_data <- kp_combined %>%
    mutate(Group = factor(in_keeppath_top, 
                        labels = c("Size-Matched Control", "KeepPathSize Top Pathways"))) %>%
    select(N, Group, name, size, mean_score, evidence_density, max_score) %>%
    mutate(N = factor(N))
  
  # Combined mean score plots - BireWire
  p1 <- ggplot(bw_plot_data, aes(x=Group, y=mean_score, fill=Group)) +
    geom_boxplot(alpha=0.7) +
    geom_jitter(width=0.2, alpha=0.4) +
    stat_summary(fun=mean, geom="point", shape=23, size=3, fill="black") +
    facet_wrap(~ N, scales = "free_y", labeller = labeller(N = function(x) paste0("Top ", x, " Pathways"))) +
    scale_fill_manual(values=c("Size-Matched Control" = "gray80", "BireWire Top Pathways" = "lightblue")) +
    labs(title="BireWire: OpenTargets Mean Score Comparison Across Different N Values",
         x="", y="Mean OpenTargets Score") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Combined mean score plots - KeepPathSize
  p2 <- ggplot(kp_plot_data, aes(x=Group, y=mean_score, fill=Group)) +
    geom_boxplot(alpha=0.7) +
    geom_jitter(width=0.2, alpha=0.4) +
    stat_summary(fun=mean, geom="point", shape=23, size=3, fill="black") +
    facet_wrap(~ N, scales = "free_y", labeller = labeller(N = function(x) paste0("Top ", x, " Pathways"))) +
    scale_fill_manual(values=c("Size-Matched Control" = "gray80", "KeepPathSize Top Pathways" = "lightpink")) +
    labs(title="KeepPathSize: OpenTargets Mean Score Comparison Across Different N Values",
         x="", y="Mean OpenTargets Score") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Evidence density plots
  p3 <- ggplot(bw_plot_data, aes(x=Group, y=evidence_density, fill=Group)) +
    geom_boxplot(alpha=0.7) +
    geom_jitter(width=0.2, alpha=0.4) +
    stat_summary(fun=mean, geom="point", shape=23, size=3, fill="black") +
    facet_wrap(~ N, scales = "free_y", labeller = labeller(N = function(x) paste0("Top ", x, " Pathways"))) +
    scale_fill_manual(values=c("Size-Matched Control" = "gray80", "BireWire Top Pathways" = "lightblue")) +
    labs(title="BireWire: Evidence Density Comparison Across Different N Values",
         x="", y="Evidence Density") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  p4 <- ggplot(kp_plot_data, aes(x=Group, y=evidence_density, fill=Group)) +
    geom_boxplot(alpha=0.7) +
    geom_jitter(width=0.2, alpha=0.4) +
    stat_summary(fun=mean, geom="point", shape=23, size=3, fill="black") +
    facet_wrap(~ N, scales = "free_y", labeller = labeller(N = function(x) paste0("Top ", x, " Pathways"))) +
    scale_fill_manual(values=c("Size-Matched Control" = "gray80", "KeepPathSize Top Pathways" = "lightpink")) +
    labs(title="KeepPathSize: Evidence Density Comparison Across Different N Values",
         x="", y="Evidence Density") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Direct comparison between methods
  # Create dataset with advantage metrics
  bw_advantage <- bw_plot_data %>%
    group_by(N) %>%
    summarise(
      mean_score_top = mean(mean_score[Group == "BireWire Top Pathways"], na.rm=TRUE),
      mean_score_control = mean(mean_score[Group == "Size-Matched Control"], na.rm=TRUE),
      mean_score_advantage = mean_score_top - mean_score_control,
      evidence_density_top = mean(evidence_density[Group == "BireWire Top Pathways"], na.rm=TRUE),
      evidence_density_control = mean(evidence_density[Group == "Size-Matched Control"], na.rm=TRUE),
      evidence_density_advantage = evidence_density_top - evidence_density_control
    ) %>%
    mutate(method = "BireWire")
  
  kp_advantage <- kp_plot_data %>%
    group_by(N) %>%
    summarise(
      mean_score_top = mean(mean_score[Group == "KeepPathSize Top Pathways"], na.rm=TRUE),
      mean_score_control = mean(mean_score[Group == "Size-Matched Control"], na.rm=TRUE),
      mean_score_advantage = mean_score_top - mean_score_control,
      evidence_density_top = mean(evidence_density[Group == "KeepPathSize Top Pathways"], na.rm=TRUE),
      evidence_density_control = mean(evidence_density[Group == "Size-Matched Control"], na.rm=TRUE),
      evidence_density_advantage = evidence_density_top - evidence_density_control
    ) %>%
    mutate(method = "KeepPathSize")
  
  advantage_data <- rbind(bw_advantage, kp_advantage)
  
  # Create advantage comparison plot
  p5 <- ggplot(advantage_data, aes(x=N, y=mean_score_advantage, fill=method)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_text(aes(label=sprintf("%.3f", mean_score_advantage)), 
              position=position_dodge(width=0.9), vjust=-0.5) +
    scale_fill_manual(values=c("BireWire"="blue", "KeepPathSize"="red")) +
    labs(title="Mean Score Advantage (Top Pathways vs. Size-Matched Controls)",
         x="Top N Pathways", y="Score Difference") +
    theme_minimal()
  
  p6 <- ggplot(advantage_data, aes(x=N, y=evidence_density_advantage, fill=method)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_text(aes(label=sprintf("%.3f", evidence_density_advantage)), 
              position=position_dodge(width=0.9), vjust=-0.5) +
    scale_fill_manual(values=c("BireWire"="blue", "KeepPathSize"="red")) +
    labs(title="Evidence Density Advantage (Top Pathways vs. Size-Matched Controls)",
         x="Top N Pathways", y="Density Difference") +
    theme_minimal()
  
  # Print all plots
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  print(p6)
  
  dev.off()
  
  # Create a summary CSV
  summary_data <- advantage_data %>%
    select(N, method, mean_score_advantage, evidence_density_advantage) %>%
    pivot_wider(
      names_from = method,
      values_from = c(mean_score_advantage, evidence_density_advantage)
    ) %>%
    mutate(
      mean_score_better_method = ifelse(mean_score_advantage_BireWire > mean_score_advantage_KeepPathSize, 
                                      "BireWire", "KeepPathSize"),
      mean_score_advantage_ratio = mean_score_advantage_BireWire / mean_score_advantage_KeepPathSize,
      evidence_better_method = ifelse(evidence_density_advantage_BireWire > evidence_density_advantage_KeepPathSize, 
                                     "BireWire", "KeepPathSize"),
      evidence_advantage_ratio = evidence_density_advantage_BireWire / evidence_density_advantage_KeepPathSize
    )
  
  write.csv(summary_data, paste0(tolower(disease_name), "_matching_advantage_summary.csv"), 
           row.names = FALSE, quote = FALSE)
  
  write.csv(advantage_data, paste0(tolower(disease_name), "_matching_raw_advantage.csv"), 
           row.names = FALSE, quote = FALSE)
}

# Run the enhanced matching analysis
matching_results <- perform_enhanced_matching(
  birewire_results, 
  keeppath_results, 
  pathway_scores, 
  n_values = c(10, 20, 50, 100)
)

# Save all matching results to CSV files
for(name in names(matching_results)) {
  write.csv(matching_results[[name]], 
           paste0(tolower(disease_name), "_size_matched_", name, ".csv"),
           row.names = FALSE, quote = FALSE)
}

# Create a statistical test table
stat_test_results <- data.frame()

for(n in c(20, 50, 100)) {
  bw_data <- matching_results[[paste0("birewire_", n)]]
  kp_data <- matching_results[[paste0("keeppath_", n)]]
  
  # BireWire t-test
  bw_test_mean_score <- t.test(
    bw_data$mean_score[bw_data$in_birewire_top],
    bw_data$mean_score[!bw_data$in_birewire_top]
  )
  
  bw_test_evidence <- t.test(
    bw_data$evidence_density[bw_data$in_birewire_top],
    bw_data$evidence_density[!bw_data$in_birewire_top]
  )
  
  # KeepPathSize t-test
  kp_test_mean_score <- t.test(
    kp_data$mean_score[kp_data$in_keeppath_top],
    kp_data$mean_score[!kp_data$in_keeppath_top]
  )
  
  kp_test_evidence <- t.test(
    kp_data$evidence_density[kp_data$in_keeppath_top],
    kp_data$evidence_density[!kp_data$in_keeppath_top]
  )
  
  # Add to results table
  stat_test_results <- rbind(stat_test_results, data.frame(
    N = n,
    Method = "BireWire",
    Metric = "Mean Score",
    t_value = bw_test_mean_score$statistic,
    p_value = bw_test_mean_score$p.value,
    significant = bw_test_mean_score$p.value < 0.05
  ))
  
  stat_test_results <- rbind(stat_test_results, data.frame(
    N = n,
    Method = "BireWire",
    Metric = "Evidence Density",
    t_value = bw_test_evidence$statistic,
    p_value = bw_test_evidence$p.value,
    significant = bw_test_evidence$p.value < 0.05
  ))
  
  stat_test_results <- rbind(stat_test_results, data.frame(
    N = n,
    Method = "KeepPathSize",
    Metric = "Mean Score",
    t_value = kp_test_mean_score$statistic,
    p_value = kp_test_mean_score$p.value,
    significant = kp_test_mean_score$p.value < 0.05
  ))
  
  stat_test_results <- rbind(stat_test_results, data.frame(
    N = n,
    Method = "KeepPathSize",
    Metric = "Evidence Density",
    t_value = kp_test_evidence$statistic,
    p_value = kp_test_evidence$p.value,
    significant = kp_test_evidence$p.value < 0.05
  ))
}

# Save statistical test results
write.csv(stat_test_results, paste0(tolower(disease_name), "_matching_statistical_tests.csv"), 
         row.names = FALSE, quote = FALSE)

cat("\nEnhanced size matching analysis complete!\n")













