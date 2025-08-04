# ===== Top Pathway Prioritization Analysis - Visualizations =====
# Enhanced to support all pathway ranking methods with default colors
library(tidyverse)
library(ggplot2)
library(patchwork)
library(data.table)
library(gridExtra)

# Process command line arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 4) {
  stop("Usage: Rscript size_matched_OT_viz.R <trait> <tool_base> <detailed_advantage_file> <sig_threshold> <data_dir>")
}

trait <- args[1]
tool_base <- args[2]
detailed_advantage_file <- args[3]
sig_threshold <- as.numeric(args[4])
data_dir <- args[5]

# Load advantage data
cat("Loading advantage data from:", detailed_advantage_file, "\n")
advantage_data <- read.csv(detailed_advantage_file)
advantage_data$N <- as.factor(as.numeric(advantage_data$N))

# Extract n_values from advantage data
n_values <- unique(as.numeric(as.character(advantage_data$N)))
cat("Found N values:", paste(n_values, collapse=", "), "\n")

# Method labels for plotting
method_labels <- list(
  BireWire_empP = "BireWire",
  KeepPathSize_empP = "KeepPathSize",
  RawP = "P-value Only",
  PvalueBeta = "P-value & Beta",
  BireWire_EmpPvalStdBeta = "BiReWire: Empirical p-value & Std Beta",  # New method
  KeepPathSize_EmpPvalStdBeta = "KeepPathSize: Empirical p-value & Std Beta"  # New method
)

# Update method names in advantage data to use descriptive labels
cat("Original methods in advantage_data:", paste(unique(advantage_data$method), collapse=", "), "\n")

# Transform method names to descriptive labels
advantage_data <- advantage_data %>%
  mutate(method = case_when(
    method == "BireWire_empP" ~ method_labels$BireWire_empP,
    method == "KeepPathSize_empP" ~ method_labels$KeepPathSize_empP,
    method == "RawP" ~ method_labels$RawP,
    method == "PvalueBeta" ~ method_labels$PvalueBeta,
    method == "BireWire_EmpPvalStdBeta" ~ method_labels$BireWire_EmpPvalStdBeta,
    method == "KeepPathSize_EmpPvalStdBeta" ~ method_labels$KeepPathSize_EmpPvalStdBeta,
    TRUE ~ method
  ))

cat("Transformed methods in advantage_data:", paste(unique(advantage_data$method), collapse=", "), "\n")

# === PART 1: CREATE ADVANTAGE PLOTS ===
cat("\nCreating advantage plots...\n")

# Add significance indicators with stricter threshold
advantage_data <- advantage_data %>%
  mutate(
    # Create significance indicator for mean score - only if p < 0.00625
    score_sig_label = ifelse(p_value_mean_score < 0.00625, "*", ""),
    
    # Create significance indicator for evidence density - only if p < 0.00625
    density_sig_label = ifelse(p_value_evidence_density < 0.00625, "*", "")
  )

# Open PDF for advantage plots
pdf(paste0(trait, "_combined_advantage_plots.pdf"), width=10, height=8)

# === Create Mean Score Advantage Plot with significance markers ===
p1 <- ggplot(advantage_data, aes(x=N, y=mean_score_advantage, fill=method)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.7) +
  geom_text(aes(label=sprintf("%.3f", mean_score_advantage)), 
            position=position_dodge2(width=0.9,padding=0.5), vjust=-0.5, size=4) +
  # Add significance markers only where p<0.00625
  geom_text(aes(label=score_sig_label, group=method),
            position=position_dodge2(width=0.9,padding=0.5),
            vjust=-1.8, size=6, fontface="bold") +
  # Use default ggplot colors
  labs(title="Mean Score Advantage (Top Pathways vs. Size-Matched Controls)", 
       x="Top N Pathways", y="Score Difference",
       subtitle="* p<0.00625 (Bonferroni corrected)") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face="bold"),
    axis.title = element_text(size=14, face="bold"),
    axis.text = element_text(size=12),
    plot.title = element_text(size=16, face="bold"),
    plot.subtitle = element_text(size=10, face="italic"),
    legend.text = element_text(size=12)
  ) +
  guides(fill=guide_legend(title="Null Model"))
print(p1)

# === Create Evidence Density Advantage Plot with significance markers ===
p2 <- ggplot(advantage_data, aes(x=N, y=evidence_density_advantage, fill=method)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.7) +
  geom_text(aes(label=sprintf("%.3f", evidence_density_advantage)), 
            position=position_dodge2(width=0.9,padding=0.5), vjust=-0.5, size=4) +
  # Add significance markers only where p<0.00625
  geom_text(aes(label=density_sig_label, group=method),
            position=position_dodge2(width=0.9,padding=0.5),
            vjust=-1.8, size=6, fontface="bold") +
  # Use default ggplot colors
  labs(title="Evidence Density Advantage (Top Pathways vs. Size-Matched Controls)", 
       x="Top N Pathways", y="Density Difference",
       subtitle="* p<0.00625 (Bonferroni corrected)") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face="bold"),
    axis.title = element_text(size=14, face="bold"),
    axis.text = element_text(size=12),
    plot.title = element_text(size=16, face="bold"),
    plot.subtitle = element_text(size=10, face="italic"),
    legend.text = element_text(size=12)
  ) +
  guides(fill=guide_legend(title="Null Model"))
print(p2)

dev.off()
cat("Created advantage plots:", trait, "_combined_advantage_plots.pdf\n")

# === PART 2: CREATE COMBINED BOXPLOTS ===
cat("\nCreating combined boxplots...\n")

# Load matched data files for all methods
all_bw_empp_files <- list.files(path = data_dir, 
                          pattern = paste0(trait, "_", tool_base, "_n\\d+_birewire_empp.csv"), 
                          full.names = TRUE)
all_kp_empp_files <- list.files(path = data_dir, 
                          pattern = paste0(trait, "_", tool_base, "_n\\d+_keeppathsize_empp.csv"), 
                          full.names = TRUE)
all_rawp_files <- list.files(path = data_dir, 
                            pattern = paste0(trait, "_", tool_base, "_n\\d+_rawp_raw_p.csv"), 
                            full.names = TRUE)
all_pvalbeta_files <- list.files(path = data_dir, 
                               pattern = paste0(trait, "_", tool_base, "_n\\d+_pvaluebeta_p_beta.csv"), 
                               full.names = TRUE)
all_bw_emppstdbeta_files <- list.files(path = data_dir, 
                               pattern = paste0(trait, "_", tool_base, "_n\\d+_birewire_emppvalstdbeta_p_beta.csv"), 
                               full.names = TRUE)
all_kp_emppstdbeta_files <- list.files(path = data_dir, 
                               pattern = paste0(trait, "_", tool_base, "_n\\d+_keeppathsize_emppvalstdbeta_p_beta.csv"), 
                               full.names = TRUE)

# Check if we have any data files
if(length(all_bw_empp_files) + length(all_kp_empp_files) + length(all_rawp_files) + length(all_pvalbeta_files) + length(all_bw_emppstdbeta_files) + length(all_kp_emppstdbeta_files) == 0) {
  cat("Warning: No matched data files found for any method\n")
} else {
  # Initialize combined data frames for each method
  combined_bw_empp_data <- data.frame()
  combined_kp_empp_data <- data.frame()
  combined_rawp_data <- data.frame()
  combined_pvalbeta_data <- data.frame()
  combined_bw_emppstdbeta_data <- data.frame()
  combined_kp_emppstdbeta_data <- data.frame()

  for(n in n_values) {
    cat("\nProcessing data for N =", n, "...\n")
    
    # Find files for this n value across all methods
    bw_empp_file <- grep(paste0("_n", n, "_birewire_empp"), all_bw_empp_files, value = TRUE)
    kp_empp_file <- grep(paste0("_n", n, "_keeppath"), all_kp_empp_files, value = TRUE)
    rawp_file <- grep(paste0("_n", n, "_rawp"), all_rawp_files, value = TRUE)
    pvalbeta_file <- grep(paste0("_n", n, "_pvaluebeta"), all_pvalbeta_files, value = TRUE)
    bw_emppstdbeta_file <- grep(paste0("_n", n, "_birewire_emppvalstdbeta"), all_bw_emppstdbeta_files, value = TRUE)
    kp_emppstdbeta_file <- grep(paste0("_n", n, "_keeppathsize_emppvalstdbeta"), all_kp_emppstdbeta_files, value = TRUE)

    # Process BireWire Empirical P-Value data if available
    if(length(bw_empp_file) > 0) {
      cat("Processing BireWire Empirical P-Value data from", basename(bw_empp_file[1]), "\n")
      bw_empp_data <- read.csv(bw_empp_file[1])
      
      # Check if in_birewire_top exists
      if("in_birewire_top" %in% colnames(bw_empp_data)) {
        # Add N and Group columns
        bw_empp_data$N <- paste0("Top ", n)
        bw_empp_data$Group <- factor(bw_empp_data$in_birewire_top, 
                               labels = c("Control", method_labels$BireWire_empP))
        bw_empp_data$Method_Short <- "BireWire_empP"
        combined_bw_empp_data <- rbind(combined_bw_empp_data, bw_empp_data)
        
        # Calculate stats for debugging
        random_scores <- bw_empp_data$mean_score[bw_empp_data$in_birewire_top]
        control_scores <- bw_empp_data$mean_score[!bw_empp_data$in_birewire_top]
        t_test <- t.test(random_scores, control_scores)
        cat(method_labels$BireWire_empP, "advantage:", 
            mean(random_scores, na.rm=TRUE) - mean(control_scores, na.rm=TRUE), 
            "(p =", t_test$p.value, ")\n")
      } else {
        cat("Warning: in_birewire_top column not found in BireWire data\n")
      }
    }
    
    # Process KeepPathSize Empirical P-Value data if available
    if(length(kp_empp_file) > 0) {
      cat("Processing KeepPathSize Empirical P-Value data from", basename(kp_empp_file[1]), "\n")
      kp_empp_data <- read.csv(kp_empp_file[1])
      
      # Check if in_keeppath_top exists
      if("in_keeppath_top" %in% colnames(kp_empp_data)) {
        # Add N and Group columns
        kp_empp_data$N <- paste0("Top ", n)
        kp_empp_data$Group <- factor(kp_empp_data$in_keeppath_top, 
                              labels = c("Control", method_labels$KeepPathSize_empP))
        kp_empp_data$Method_Short <- "KeepPathSize_empP"
        combined_kp_empp_data <- rbind(combined_kp_empp_data, kp_empp_data)
        
        # Calculate stats for debugging
        random_scores <- kp_empp_data$mean_score[kp_empp_data$in_keeppath_top]
        control_scores <- kp_empp_data$mean_score[!kp_empp_data$in_keeppath_top]
        t_test <- t.test(random_scores, control_scores)
        cat(method_labels$KeepPathSize_empP, "advantage:", 
            mean(random_scores, na.rm=TRUE) - mean(control_scores, na.rm=TRUE), 
            "(p =", t_test$p.value, ")\n")
      } else {
        cat("Warning: in_keeppath_top column not found in KeepPathSize data\n")
      }
    }
    
    # Process Raw P-value data if available
    if(length(rawp_file) > 0) {
      cat("Processing Raw P-value data from", basename(rawp_file[1]), "\n")
      rawp_data <- read.csv(rawp_file[1])
      
      # Check if in_raw_p_top exists
      if("in_raw_p_top" %in% colnames(rawp_data)) {
        # Add N and Group columns
        rawp_data$N <- paste0("Top ", n)
        rawp_data$Group <- factor(rawp_data$in_raw_p_top, labels = c("Control", method_labels$RawP))
        rawp_data$Method_Short <- "RawP"
        combined_rawp_data <- rbind(combined_rawp_data, rawp_data)
        
        # Calculate stats for debugging
        random_scores <- rawp_data$mean_score[rawp_data$in_raw_p_top]
        control_scores <- rawp_data$mean_score[!rawp_data$in_raw_p_top]
        t_test <- t.test(random_scores, control_scores)
        cat(method_labels$RawP, "advantage:", 
            mean(random_scores, na.rm=TRUE) - mean(control_scores, na.rm=TRUE), 
            "(p =", t_test$p.value, ")\n")
      } else {
        cat("Warning: in_raw_p_top column not found in Raw P-value data\n")
      }
    }
    
    # Process P-Value + Beta data if available
    if(length(pvalbeta_file) > 0) {
      cat("Processing P-Value + Beta data from", basename(pvalbeta_file[1]), "\n")
      pvalbeta_data <- read.csv(pvalbeta_file[1])
      
      # Check if in_p_beta_top exists
      if("in_p_beta_top" %in% colnames(pvalbeta_data)) {
        # Add N and Group columns
        pvalbeta_data$N <- paste0("Top ", n)
        pvalbeta_data$Group <- factor(pvalbeta_data$in_p_beta_top, 
                                 labels = c("Control", method_labels$PvalueBeta))
        pvalbeta_data$Method_Short <- "PvalueBeta"
        combined_pvalbeta_data <- rbind(combined_pvalbeta_data, pvalbeta_data)
        
        # Calculate stats for debugging
        random_scores <- pvalbeta_data$mean_score[pvalbeta_data$in_p_beta_top]
        control_scores <- pvalbeta_data$mean_score[!pvalbeta_data$in_p_beta_top]
        t_test <- t.test(random_scores, control_scores)
        cat(method_labels$PvalueBeta, "advantage:", 
            mean(random_scores, na.rm=TRUE) - mean(control_scores, na.rm=TRUE), 
            "(p =", t_test$p.value, ")\n")
      } else {
        cat("Warning: in_p_beta_top column not found in P-Value + Beta data\n")
      }
    }
    
    # Process BireWire Empirical P-Value + Std Beta data if available
    if(length(bw_emppstdbeta_file) > 0) {
      cat("Processing BireWire Empirical P-Value + Std Beta data from", basename(bw_emppstdbeta_file[1]), "\n")
      bw_emppstdbeta_data <- read.csv(bw_emppstdbeta_file[1])
      
      # Check if in_p_beta_top exists (common column name across beta methods)
      if("in_p_beta_top" %in% colnames(bw_emppstdbeta_data)) {
        # Add N and Group columns
        bw_emppstdbeta_data$N <- paste0("Top ", n)
        bw_emppstdbeta_data$Group <- factor(bw_emppstdbeta_data$in_p_beta_top, 
                                       labels = c("Control", method_labels$BireWire_EmpPvalStdBeta))
        bw_emppstdbeta_data$Method_Short <- "BW_EmpPval_StdBeta"
        combined_bw_emppstdbeta_data <- rbind(combined_bw_emppstdbeta_data, bw_emppstdbeta_data)
        
        # Calculate stats for debugging
        random_scores <- bw_emppstdbeta_data$mean_score[bw_emppstdbeta_data$in_p_beta_top]
        control_scores <- bw_emppstdbeta_data$mean_score[!bw_emppstdbeta_data$in_p_beta_top]
        t_test <- t.test(random_scores, control_scores)
        cat(method_labels$BireWire_EmpPvalStdBeta, "advantage:", 
            mean(random_scores, na.rm=TRUE) - mean(control_scores, na.rm=TRUE), 
            "(p =", t_test$p.value, ")\n")
      } else {
        cat("Warning: in_p_beta_top column not found in BireWire Empirical P-Value + Std Beta data\n")
      }
    }
    
    # Process KeepPathSize Empirical P-Value + Std Beta data if available
    if(length(kp_emppstdbeta_file) > 0) {
      cat("Processing KeepPathSize Empirical P-Value + Std Beta data from", basename(kp_emppstdbeta_file[1]), "\n")
      kp_emppstdbeta_data <- read.csv(kp_emppstdbeta_file[1])
      
      # Check if in_p_beta_top exists (common column name across beta methods)
      if("in_p_beta_top" %in% colnames(kp_emppstdbeta_data)) {
        # Add N and Group columns
        kp_emppstdbeta_data$N <- paste0("Top ", n)
        kp_emppstdbeta_data$Group <- factor(kp_emppstdbeta_data$in_p_beta_top, 
                                      labels = c("Control", method_labels$KeepPathSize_EmpPvalStdBeta))
        kp_emppstdbeta_data$Method_Short <- "KP_EmpPval_StdBeta"
        combined_kp_emppstdbeta_data <- rbind(combined_kp_emppstdbeta_data, kp_emppstdbeta_data)
        
        # Calculate stats for debugging
        random_scores <- kp_emppstdbeta_data$mean_score[kp_emppstdbeta_data$in_p_beta_top]
        control_scores <- kp_emppstdbeta_data$mean_score[!kp_emppstdbeta_data$in_p_beta_top]
        t_test <- t.test(random_scores, control_scores)
        cat(method_labels$KeepPathSize_EmpPvalStdBeta, "advantage:", 
            mean(random_scores, na.rm=TRUE) - mean(control_scores, na.rm=TRUE), 
            "(p =", t_test$p.value, ")\n")
      } else {
        cat("Warning: in_p_beta_top column not found in KeepPathSize Empirical P-Value + Std Beta data\n")
      }
    }
  }
  
  # Order N factor correctly for all datasets
  for(df_name in c("combined_bw_empp_data", "combined_kp_empp_data", "combined_rawp_data", 
                  "combined_pvalbeta_data", "combined_bw_emppstdbeta_data", "combined_kp_emppstdbeta_data")) {
    df <- get(df_name)
    if(!is.null(df) && nrow(df) > 0) {
      df$N <- factor(df$N, levels = paste0("Top ", n_values))
      assign(df_name, df)
    }
  }
  
  # Combine all datasets for direct comparison plots
  # Extract only top pathways (not controls) from each method
  top_pathways_list <- list()

  if(nrow(combined_bw_empp_data) > 0) {
    top_bw_empp <- combined_bw_empp_data[combined_bw_empp_data$Group == method_labels$BireWire_empP, ]
    top_pathways_list[["BireWire_empP"]] <- top_bw_empp
  }

  if(nrow(combined_kp_empp_data) > 0) {
    top_kp_empp <- combined_kp_empp_data[combined_kp_empp_data$Group == method_labels$KeepPathSize_empP, ]
    top_pathways_list[["KeepPathSize_empP"]] <- top_kp_empp
  }

  if(nrow(combined_rawp_data) > 0) {
    top_rawp <- combined_rawp_data[combined_rawp_data$Group == method_labels$RawP, ]
    top_pathways_list[["RawP"]] <- top_rawp
  }

  if(nrow(combined_pvalbeta_data) > 0) {
    top_pvalbeta <- combined_pvalbeta_data[combined_pvalbeta_data$Group == method_labels$PvalueBeta, ]
    top_pathways_list[["PvalueBeta"]] <- top_pvalbeta
  }

  if(nrow(combined_bw_emppstdbeta_data) > 0) {
    top_bw_emppstdbeta <- combined_bw_emppstdbeta_data[combined_bw_emppstdbeta_data$Group == method_labels$BireWire_EmpPvalStdBeta, ]
    top_pathways_list[["BW_EmpPval_StdBeta"]] <- top_bw_emppstdbeta
  }

  if(nrow(combined_kp_emppstdbeta_data) > 0) {
    top_kp_emppstdbeta <- combined_kp_emppstdbeta_data[combined_kp_emppstdbeta_data$Group == method_labels$KeepPathSize_EmpPvalStdBeta, ]
    top_pathways_list[["KP_EmpPval_StdBeta"]] <- top_kp_emppstdbeta
  }
  
  # Combine all top pathways into one dataframe
  top_combined <- do.call(rbind, top_pathways_list)
  
  # Create combined boxplots and save to PDF
  pdf(paste0(trait, "_combined_boxplots.pdf"), width=12, height=10)
  
  # Increase text sizes for all plots
  title_size <- 18
  axis_title_size <- 16
  axis_text_size <- 14
  legend_title_size <- 14
  legend_text_size <- 14
  
  # 1. Create individual method plots
  for(method_data in list(
    list(data=combined_bw_empp_data, name="BireWire_empP", label=method_labels$BireWire_empP),
    list(data=combined_kp_empp_data, name="KeepPathSize_empP", label=method_labels$KeepPathSize_empP),
    list(data=combined_rawp_data, name="RawP", label=method_labels$RawP),
    list(data=combined_pvalbeta_data, name="PvalueBeta", label=method_labels$PvalueBeta),
    list(data=combined_bw_emppstdbeta_data, name="BW_EmpPval_StdBeta", label=method_labels$BireWire_EmpPvalStdBeta),
    list(data=combined_kp_emppstdbeta_data, name="KP_EmpPval_StdBeta", label=method_labels$KeepPathSize_EmpPvalStdBeta)
  )) {
    
    if(!is.null(method_data$data) && nrow(method_data$data) > 0) {
      # Mean score plot
      p <- ggplot(method_data$data, aes(x=N, y=mean_score, fill=Group)) +
        geom_boxplot(position=position_dodge(), width=0.7, alpha=0.8, outlier.shape=NA) +
        geom_point(position=position_jitterdodge(jitter.width=0.15, dodge.width=0.7), 
                  alpha=0.4, size=1) +
        # Use default ggplot colors
        labs(title=paste(method_data$label, ": OpenTargets Mean Score Across Different N Values"),
             x="", y="Mean OpenTargets Score") +
        theme_minimal() +
        theme(
          legend.position = "bottom",
          axis.text.x = element_text(angle = 0, hjust = 0.5, size=axis_text_size),
          axis.text.y = element_text(size=axis_text_size),
          legend.title = element_blank(),
          axis.title = element_text(size=axis_title_size, face="bold"),
          plot.title = element_text(size=title_size, face="bold"),
          legend.text = element_text(size=legend_text_size)
        )
      print(p)
      
      # Evidence density plot
      p <- ggplot(method_data$data, aes(x=N, y=evidence_density, fill=Group)) +
        geom_boxplot(position=position_dodge(), width=0.7, alpha=0.8, outlier.shape=NA) +
        geom_point(position=position_jitterdodge(jitter.width=0.15, dodge.width=0.7), 
                  alpha=0.4, size=1) +
        # Use default ggplot colors
        labs(title=paste(method_data$label, ": Evidence Density Across Different N Values"),
             x="", y="Evidence Density") +
        theme_minimal() +
        theme(
          legend.position = "bottom",
          axis.text.x = element_text(angle = 0, hjust = 0.5, size=axis_text_size),
          axis.text.y = element_text(size=axis_text_size),
          legend.title = element_blank(),
          axis.title = element_text(size=axis_title_size, face="bold"),
          plot.title = element_text(size=title_size, face="bold"),
          legend.text = element_text(size=legend_text_size)
        )
      print(p)
    }
  }
  
  # 2. Direct comparison of all methods (if we have more than one method with data)
  if(length(top_pathways_list) > 1) {
    cat("Creating direct method comparison plots with", length(top_pathways_list), "methods\n")
    
    # Direct comparison of methods - Mean Score
    p <- ggplot(top_combined, aes(x=N, y=mean_score, fill=Method_Short)) +
      geom_boxplot(position=position_dodge(), width=0.7, alpha=0.8, outlier.shape=NA) +
      geom_point(position=position_jitterdodge(jitter.width=0.15, dodge.width=0.7), 
                alpha=0.4, size=1) +
      # Use default ggplot colors
      labs(title="Direct Comparison: OpenTargets Mean Score by Method",
           x="", y="Mean OpenTargets Score") +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 0, hjust = 0.5, size=axis_text_size),
        axis.text.y = element_text(size=axis_text_size),
        legend.title = element_text(face="bold", size=legend_title_size),
        axis.title = element_text(size=axis_title_size, face="bold"),
        plot.title = element_text(size=title_size, face="bold"),
        legend.text = element_text(size=legend_text_size)
      ) +
      guides(fill=guide_legend(title="Null Model"))
    print(p)
    
    # Direct comparison of methods - Evidence Density
    p <- ggplot(top_combined, aes(x=N, y=evidence_density, fill=Method_Short)) +
      geom_boxplot(position=position_dodge(), width=0.7, alpha=0.8, outlier.shape=NA) +
      geom_point(position=position_jitterdodge(jitter.width=0.15, dodge.width=0.7), 
                alpha=0.4, size=1) +
      # Use default ggplot colors
      labs(title="Direct Comparison: Evidence Density by Method",
           x="", y="Evidence Density") +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 0, hjust = 0.5, size=axis_text_size),
        axis.text.y = element_text(size=axis_text_size),
        legend.title = element_text(face="bold", size=legend_title_size),
        axis.title = element_text(size=axis_title_size, face="bold"),
        plot.title = element_text(size=title_size, face="bold"),
        legend.text = element_text(size=legend_text_size)
      ) +
      guides(fill=guide_legend(title="Null Model"))
    print(p)
  }
  
  dev.off()
  cat("Created combined boxplots:", trait, "_combined_boxplots.pdf\n")
  }

cat("\n===== All visualizations complete =====\n")