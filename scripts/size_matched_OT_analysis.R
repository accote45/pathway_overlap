# ===== Top Pathway Prioritization Analysis =====
# This script compares the biological relevance of top-ranked pathways 
# identified by different randomization methods

library(plyr)
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
  
  combined_dat <- do.call(base::rbind, dat[!sapply(dat, is.null)])
  
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
    mean_score = mean(score, na.rm = TRUE),
    num_genes = n(),
    num_with_evidence = sum(score > 0, na.rm = TRUE),
    evidence_density = num_with_evidence / num_genes,
    max_score = max(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE)
  ) %>% 
  as.data.frame()


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
          
          bw_comparison <- base::rbind(bw_comparison, data.frame(
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
          
          kp_comparison <- base::rbind(kp_comparison, data.frame(
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
  
# Convert to data frames before combining (avoiding MatchIt's rbind method)
  bw_names <- grep("^birewire_\\d+$", names(all_results), value = TRUE)
  kp_names <- grep("^keeppath_\\d+$", names(all_results), value = TRUE)
  
  # Convert each matched object to a data frame before combining
  bw_combined <- do.call(rbind, lapply(all_results[bw_names], as.data.frame))
  kp_combined <- do.call(rbind, lapply(all_results[kp_names], as.data.frame))
  
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

# Helper function to create combined plots across all N values - with grouped boxplots
create_combined_plots <- function(bw_combined, kp_combined, disease_name) {
  pdf(paste0(tolower(disease_name), "_combined_matching_analysis.pdf"), width=14, height=10)
  
  # Define a larger text size for all plots
  bigger_text <- theme(
    axis.title = element_text(size=18, face="bold"),
    axis.text = element_text(size=16),
    plot.title = element_text(size=20, face="bold"),
    legend.text = element_text(size=16),
    legend.title = element_text(size=16)
  )
  
  # Define new method labels
  new_bw_label <- "Fix gene frequency and pathway size"
  new_kp_label <- "Fix pathway size"
  
  # Reshape BireWire data for grouped boxplots with new labels
  bw_plot_data <- bw_combined %>%
    mutate(Group = factor(in_birewire_top, 
                        labels = c("Control", new_bw_label)),
           N = factor(N, levels = sort(unique(N)), 
                     labels = paste0("Top ", sort(unique(N))))) %>%
    select(N, Group, name, size, mean_score, evidence_density, max_score)
  
  # Reshape KeepPathSize data for grouped boxplots with new labels
  kp_plot_data <- kp_combined %>%
    mutate(Group = factor(in_keeppath_top, 
                        labels = c("Control", new_kp_label)),
           N = factor(N, levels = sort(unique(N)), 
                     labels = paste0("Top ", sort(unique(N))))) %>%
    select(N, Group, name, size, mean_score, evidence_density, max_score)
  
  # Calculate t-test results for each N value (keep for data export, but don't show on plots)
  n_values <- levels(bw_plot_data$N)
  
  # Create data frames to store significance results
  bw_sig_mean <- data.frame(N=n_values, p_value=numeric(length(n_values)), 
                           significance=character(length(n_values)), stringsAsFactors=FALSE)
  bw_sig_evidence <- data.frame(N=n_values, p_value=numeric(length(n_values)), 
                              significance=character(length(n_values)), stringsAsFactors=FALSE)
  kp_sig_mean <- data.frame(N=n_values, p_value=numeric(length(n_values)), 
                           significance=character(length(n_values)), stringsAsFactors=FALSE)
  kp_sig_evidence <- data.frame(N=n_values, p_value=numeric(length(n_values)), 
                              significance=character(length(n_values)), stringsAsFactors=FALSE)
  
  # Calculate significance for each N value
  for(i in 1:length(n_values)) {
    n_label <- n_values[i]
    
    # BireWire mean score
    bw_mean_data <- bw_plot_data %>% filter(N == n_label)
    bw_test_mean <- t.test(mean_score ~ Group, data = bw_mean_data)
    bw_sig_mean$p_value[i] <- bw_test_mean$p.value
    bw_sig_mean$significance[i] <- ifelse(bw_test_mean$p.value < 0.05, "*", "")
    
    # BireWire evidence density
    bw_test_evidence <- t.test(evidence_density ~ Group, data = bw_mean_data)
    bw_sig_evidence$p_value[i] <- bw_test_evidence$p.value
    bw_sig_evidence$significance[i] <- ifelse(bw_test_evidence$p.value < 0.05, "*", "")
    
    # KeepPathSize mean score
    kp_mean_data <- kp_plot_data %>% filter(N == n_label)
    kp_test_mean <- t.test(mean_score ~ Group, data = kp_mean_data)
    kp_sig_mean$p_value[i] <- kp_test_mean$p.value
    kp_sig_mean$significance[i] <- ifelse(kp_test_mean$p.value < 0.05, "*", "")
    
    # KeepPathSize evidence density
    kp_test_evidence <- t.test(evidence_density ~ Group, data = kp_mean_data)
    kp_sig_evidence$p_value[i] <- kp_test_evidence$p.value
    kp_sig_evidence$significance[i] <- ifelse(kp_test_evidence$p.value < 0.05, "*", "")
  }
  
  # BireWire Mean Score - Grouped boxplot (without significance annotation)
  p1 <- ggplot(bw_plot_data, aes(x=N, y=mean_score, fill=Group)) +
    geom_boxplot(position=position_dodge(width=0.8), width=0.7, alpha=0.8, outlier.shape=NA) +
    geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8), 
              alpha=0.6, size=1) +
    scale_fill_manual(values=c("Control" = "gray80", new_bw_label = "lightblue")) +
    labs(title=paste0(new_bw_label, ": OpenTargets Mean Score Across Different N Values"),
         x="Top N Pathways", y="Mean OpenTargets Score") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.title = element_text(face="bold")) +
    bigger_text +
    guides(fill=guide_legend(title="Null Model"))
  
  # KeepPathSize Mean Score - Grouped boxplot (without significance annotation)
  p2 <- ggplot(kp_plot_data, aes(x=N, y=mean_score, fill=Group)) +
    geom_boxplot(position=position_dodge(width=0.8), width=0.7, alpha=0.8, outlier.shape=NA) +
    geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8), 
              alpha=0.6, size=1) +
    scale_fill_manual(values=c("Control" = "gray80", new_kp_label = "lightpink")) +
    labs(title=paste0(new_kp_label, ": OpenTargets Mean Score Across Different N Values"),
         x="Top N Pathways", y="Mean OpenTargets Score") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.title = element_text(face="bold")) +
    bigger_text +
    guides(fill=guide_legend(title="Null Model"))
  
  # BireWire Evidence Density - Grouped boxplot (without significance annotation)
  p3 <- ggplot(bw_plot_data, aes(x=N, y=evidence_density, fill=Group)) +
    geom_boxplot(position=position_dodge(width=0.8), width=0.7, alpha=0.8, outlier.shape=NA) +
    geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8), 
              alpha=0.6, size=1) +
    scale_fill_manual(values=c("Control" = "gray80", new_bw_label = "lightblue")) +
    labs(title=paste0(new_bw_label, ": Evidence Density Across Different N Values"),
         x="Top N Pathways", y="Evidence Density") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.title = element_text(face="bold")) +
    bigger_text +
    guides(fill=guide_legend(title="Null Model"))
  
  # KeepPathSize Evidence Density - Grouped boxplot (without significance annotation)
  p4 <- ggplot(kp_plot_data, aes(x=N, y=evidence_density, fill=Group)) +
    geom_boxplot(position=position_dodge(width=0.8), width=0.7, alpha=0.8, outlier.shape=NA) +
    geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8), 
              alpha=0.6, size=1) +
    scale_fill_manual(values=c("Control" = "gray80", new_kp_label = "lightpink")) +
    labs(title=paste0(new_kp_label, ": Evidence Density Across Different N Values"),
         x="Top N Pathways", y="Evidence Density") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.title = element_text(face="bold")) +
    bigger_text +
    guides(fill=guide_legend(title="Null Model"))
  
  # Direct comparison - Create dataset with advantage metrics
  bw_advantage <- bw_plot_data %>%
    group_by(N) %>%
    summarise(
      mean_score_top = mean(mean_score[Group == new_bw_label], na.rm=TRUE),
      mean_score_control = mean(mean_score[Group == "Control"], na.rm=TRUE),
      mean_score_advantage = mean_score_top - mean_score_control,
      evidence_density_top = mean(evidence_density[Group == new_bw_label], na.rm=TRUE),
      evidence_density_control = mean(evidence_density[Group == "Control"], na.rm=TRUE),
      evidence_density_advantage = evidence_density_top - evidence_density_control
    ) %>%
    mutate(method = new_bw_label)
  
  kp_advantage <- kp_plot_data %>%
    group_by(N) %>%
    summarise(
      mean_score_top = mean(mean_score[Group == new_kp_label], na.rm=TRUE),
      mean_score_control = mean(mean_score[Group == "Control"], na.rm=TRUE),
      mean_score_advantage = mean_score_top - mean_score_control,
      evidence_density_top = mean(evidence_density[Group == new_kp_label], na.rm=TRUE),
      evidence_density_control = mean(evidence_density[Group == "Control"], na.rm=TRUE),
      evidence_density_advantage = evidence_density_top - evidence_density_control
    ) %>%
    mutate(method = new_kp_label)
  
  # Combine advantage data and add significance information
  advantage_data <- base::rbind(bw_advantage, kp_advantage)
  
  # Add significance information to advantage data (for data export only)
  advantage_data <- advantage_data %>%
    mutate(
      mean_score_significance = case_when(
        method == new_bw_label ~ bw_sig_mean$significance[match(N, bw_sig_mean$N)],
        method == new_kp_label ~ kp_sig_mean$significance[match(N, kp_sig_mean$N)],
        TRUE ~ ""
      ),
      evidence_density_significance = case_when(
        method == new_bw_label ~ bw_sig_evidence$significance[match(N, bw_sig_evidence$N)],
        method == new_kp_label ~ kp_sig_evidence$significance[match(N, kp_sig_evidence$N)],
        TRUE ~ ""
      )
    )
  
  # Create advantage comparison plots with larger text (without significance annotation on plot)
  p5 <- ggplot(advantage_data, aes(x=N, y=mean_score_advantage, fill=method)) +
    geom_bar(stat="identity", position=position_dodge(), width=0.7) +
    geom_text(aes(label=sprintf("%.3f", mean_score_advantage)), 
              position=position_dodge(width=0.9), vjust=-0.5, size=5) +
    scale_fill_manual(values=c(new_bw_label="blue", new_kp_label="red")) +
    labs(title="Mean Score Advantage (Top Pathways vs. Size-Matched Controls)",
         x="Top N Pathways", y="Score Difference") +
    theme_minimal() +
    bigger_text +
    theme(legend.position = "bottom",
          legend.title = element_text(face="bold")) +
    guides(fill=guide_legend(title="Null Model"))
  
  p6 <- ggplot(advantage_data, aes(x=N, y=evidence_density_advantage, fill=method)) +
    geom_bar(stat="identity", position=position_dodge(), width=0.7) +
    geom_text(aes(label=sprintf("%.3f", evidence_density_advantage)), 
              position=position_dodge(width=0.9), vjust=-0.5, size=5) +
    scale_fill_manual(values=c(new_bw_label="blue", new_kp_label="red")) +
    labs(title="Evidence Density Advantage (Top Pathways vs. Size-Matched Controls)",
         x="Top N Pathways", y="Density Difference") +
    theme_minimal() +
    bigger_text +
    theme(legend.position = "bottom",
          legend.title = element_text(face="bold")) +
    guides(fill=guide_legend(title="Null Model"))
  
  # Print only the individual method plots and advantage plots
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  print(p6)
  
  dev.off()
  
  # Update method names in the CSV files too
  method_mapping <- c("BireWire" = new_bw_label, "KeepPathSize" = new_kp_label)
  
  # Create a summary CSV with significance information
  summary_data <- advantage_data %>%
    select(N, method, mean_score_advantage, mean_score_significance, 
           evidence_density_advantage, evidence_density_significance) %>%
    pivot_wider(
      names_from = method,
      values_from = c(mean_score_advantage, evidence_density_advantage,
                     mean_score_significance, evidence_density_significance)
    )
  
  # Create column names using new method names
  new_colnames <- colnames(summary_data)
  for (old_name in names(method_mapping)) {
    new_colnames <- gsub(old_name, gsub(" ", "_", method_mapping[old_name]), new_colnames)
  }
  colnames(summary_data) <- new_colnames
  
  # Add comparison columns
  summary_data <- summary_data %>%
    mutate(
      mean_score_better_method = case_when(
        .[[paste0("mean_score_advantage_", gsub(" ", "_", new_bw_label))]] > 
          .[[paste0("mean_score_advantage_", gsub(" ", "_", new_kp_label))]] ~ new_bw_label,
        TRUE ~ new_kp_label
      ),
      mean_score_advantage_ratio = .[[paste0("mean_score_advantage_", gsub(" ", "_", new_bw_label))]] / 
                                   .[[paste0("mean_score_advantage_", gsub(" ", "_", new_kp_label))]],
      evidence_better_method = case_when(
        .[[paste0("evidence_density_advantage_", gsub(" ", "_", new_bw_label))]] > 
          .[[paste0("evidence_density_advantage_", gsub(" ", "_", new_kp_label))]] ~ new_bw_label,
        TRUE ~ new_kp_label
      ),
      evidence_advantage_ratio = .[[paste0("evidence_density_advantage_", gsub(" ", "_", new_bw_label))]] / 
                               .[[paste0("evidence_density_advantage_", gsub(" ", "_", new_kp_label))]]
    )
  
  write.csv(summary_data, paste0(tolower(disease_name), "_matching_advantage_summary.csv"), 
           row.names = FALSE, quote = FALSE)
  
  # Update method names in the raw advantage data
  advantage_data$method <- plyr::mapvalues(advantage_data$method, 
                                         from = c("BireWire", "KeepPathSize"),
                                         to = c(new_bw_label, new_kp_label))
  
  write.csv(advantage_data, paste0(tolower(disease_name), "_matching_raw_advantage.csv"), 
           row.names = FALSE, quote = FALSE)
  
  # Also save the significance testing results separately
  sig_results <- bind_rows(
    bw_sig_mean %>% mutate(Method = new_bw_label, Metric = "Mean Score"),
    bw_sig_evidence %>% mutate(Method = new_bw_label, Metric = "Evidence Density"),
    kp_sig_mean %>% mutate(Method = new_kp_label, Metric = "Mean Score"),
    kp_sig_evidence %>% mutate(Method = new_kp_label, Metric = "Evidence Density")
  )
  
  write.csv(sig_results, paste0(tolower(disease_name), "_significance_tests.csv"), 
           row.names = FALSE, quote = FALSE)
}


# Execute the size-matched analysis with the loaded data
cat("Starting size-matched analysis...\n")

# Define the N values you want to analyze
n_values <- c(10, 20, 50, 100)  # you can adjust this list

# Execute the analysis
results <- perform_enhanced_matching(
  birewire_results, 
  keeppath_results, 
  pathway_scores, 
  n_values
)

cat("Analysis complete.\n")
