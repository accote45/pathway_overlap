# ===== Top Pathway Prioritization Analysis - Visualizations =====
# This script creates visualizations for the results of pathway analysis

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

# === Plotting Functions ===

create_boxplot <- function(data, x_col, y_col, fill_col, title, x_lab, y_lab, fill_values) {
  p <- ggplot(data, aes_string(x=x_col, y=y_col, fill=fill_col)) +
    geom_boxplot(position=position_dodge(width=0.8), width=0.7, alpha=0.8, outlier.shape=NA) +
    geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8), alpha=0.6, size=1) +
    scale_fill_manual(values=fill_values) +
    labs(title=title, x=x_lab, y=y_lab) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.title = element_text(face="bold"),
      axis.title = element_text(size=14, face="bold"),
      axis.text = element_text(size=12),
      plot.title = element_text(size=16, face="bold"),
      legend.text = element_text(size=12)
    )
  
  return(p)
}

create_advantage_plot <- function(advantage_data, y_col, title, y_lab) {
  p <- ggplot(advantage_data, aes_string(x="N", y=y_col, fill="method")) +
    geom_bar(stat="identity", position=position_dodge(), width=0.7) +
    geom_text(aes(label=sprintf("%.3f", .data[[y_col]])), 
            position=position_dodge(width=0.9), vjust=-0.5, size=4) +
    scale_fill_manual(values=c("Fix gene frequency and pathway size"="blue", "Fix pathway size"="red")) +
    labs(title=title, x="Top N Pathways", y=y_lab) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face="bold"),
      axis.title = element_text(size=14, face="bold"),
      axis.text = element_text(size=12),
      plot.title = element_text(size=16, face="bold"),
      legend.text = element_text(size=12)
    ) +
    guides(fill=guide_legend(title="Null Model"))
  
  return(p)
}

create_combined_boxplots <- function(n_values, trait, tool_base, data_dir, method_labels) {
  # Load matched data files
  all_bw_files <- list.files(path = data_dir, pattern = paste0(trait, "_", tool_base, "_n\\d+_birewire_matched.csv"), full.names = TRUE)
  all_kp_files <- list.files(path = data_dir, pattern = paste0(trait, "_", tool_base, "_n\\d+_keeppath_matched.csv"), full.names = TRUE)
  
  if(length(all_bw_files) == 0 || length(all_kp_files) == 0) {
    cat("Warning: No matched data files found\n")
    return(NULL)
  }
  
  # Extract and combine data from all N values
  combined_bw_data <- data.frame()
  combined_kp_data <- data.frame()
  
  for(n in n_values) {
    # Find files for this n value
    bw_file <- grep(paste0("_n", n, "_birewire"), all_bw_files, value = TRUE)
    kp_file <- grep(paste0("_n", n, "_keeppath"), all_kp_files, value = TRUE)
    
    if(length(bw_file) > 0 && length(kp_file) > 0) {
      # Process BireWire data
      bw_data <- read.csv(bw_file[1])
      bw_data$N <- paste0("Top ", n)
      bw_data$Group <- factor(bw_data$in_birewire_top, labels = c("Control", "BireWire"))
      combined_bw_data <- rbind(combined_bw_data, bw_data)
      
      # Process KeepPathSize data
      kp_data <- read.csv(kp_file[1])
      kp_data$N <- paste0("Top ", n)
      kp_data$Group <- factor(kp_data$in_keeppath_top, labels = c("Control", "KeepPathSize"))
      combined_kp_data <- rbind(combined_kp_data, kp_data)
    }
  }
  
  # Order N factor correctly
  combined_bw_data$N <- factor(combined_bw_data$N, levels = paste0("Top ", n_values))
  combined_kp_data$N <- factor(combined_kp_data$N, levels = paste0("Top ", n_values))
  
  # Increase text sizes for all plots
  title_size <- 18
  axis_title_size <- 16
  axis_text_size <- 14
  legend_title_size <- 14
  legend_text_size <- 14
  
  # Create boxplots and save to PDF
  pdf(paste0(trait, "_combined_boxplots.pdf"), width=12, height=10)
  
  # BireWire Mean Score
  p1 <- ggplot(combined_bw_data, aes(x=N, y=mean_score, fill=Group)) +
    geom_boxplot(position=position_dodge(), width=0.7, alpha=0.8, outlier.shape=NA) +
    geom_point(position=position_jitterdodge(jitter.width=0.15, dodge.width=0.7), 
              alpha=0.4, size=1) +
    scale_fill_manual(values=c("Control"="gray80", "BireWire"="lightblue")) +
    labs(title=paste(method_labels$birewire, ": OpenTargets Mean Score Across Different N Values"),
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
  print(p1)
  
  # BireWire Evidence Density
  p2 <- ggplot(combined_bw_data, aes(x=N, y=evidence_density, fill=Group)) +
    geom_boxplot(position=position_dodge(), width=0.7, alpha=0.8, outlier.shape=NA) +
    geom_point(position=position_jitterdodge(jitter.width=0.15, dodge.width=0.7), 
              alpha=0.4, size=1) +
    scale_fill_manual(values=c("Control"="gray80", "BireWire"="lightblue")) +
    labs(title=paste(method_labels$birewire, ": Evidence Density Across Different N Values"),
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
  print(p2)
  
  # KeepPathSize Mean Score
  p3 <- ggplot(combined_kp_data, aes(x=N, y=mean_score, fill=Group)) +
    geom_boxplot(position=position_dodge(), width=0.7, alpha=0.8, outlier.shape=NA) +
    geom_point(position=position_jitterdodge(jitter.width=0.15, dodge.width=0.7), 
              alpha=0.4, size=1) +
    scale_fill_manual(values=c("Control"="gray80", "KeepPathSize"="lightpink")) +
    labs(title=paste(method_labels$keeppath, ": OpenTargets Mean Score Across Different N Values"),
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
  print(p3)
  
  # KeepPathSize Evidence Density
  p4 <- ggplot(combined_kp_data, aes(x=N, y=evidence_density, fill=Group)) +
    geom_boxplot(position=position_dodge(), width=0.7, alpha=0.8, outlier.shape=NA) +
    geom_point(position=position_jitterdodge(jitter.width=0.15, dodge.width=0.7), 
              alpha=0.4, size=1) +
    scale_fill_manual(values=c("Control"="gray80", "KeepPathSize"="lightpink")) +
    labs(title=paste(method_labels$keeppath, ": Evidence Density Across Different N Values"),
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
  print(p4)
  
  # Now create a combined plot with both methods for direct comparison
  # First, rename the groups to include method name
  combined_bw_data$Method_Group <- "BireWire"
  combined_kp_data$Method_Group <- "KeepPathSize"
  
  # Filter only the top pathways (not controls)
  top_bw <- combined_bw_data[combined_bw_data$Group == "BireWire", ]
  top_kp <- combined_kp_data[combined_kp_data$Group == "KeepPathSize", ]
  
  # Combine the datasets
  top_combined <- rbind(top_bw, top_kp)
  
  # Direct comparison of methods
  p5 <- ggplot(top_combined, aes(x=N, y=mean_score, fill=Method_Group)) +
    geom_boxplot(position=position_dodge(), width=0.7, alpha=0.8, outlier.shape=NA) +
    geom_point(position=position_jitterdodge(jitter.width=0.15, dodge.width=0.7), 
              alpha=0.4, size=1) +
    scale_fill_manual(values=c("BireWire"="lightblue", "KeepPathSize"="lightpink")) +
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
  print(p5)
  
  p6 <- ggplot(top_combined, aes(x=N, y=evidence_density, fill=Method_Group)) +
    geom_boxplot(position=position_dodge(), width=0.7, alpha=0.8, outlier.shape=NA) +
    geom_point(position=position_jitterdodge(jitter.width=0.15, dodge.width=0.7), 
              alpha=0.4, size=1) +
    scale_fill_manual(values=c("BireWire"="lightblue", "KeepPathSize"="lightpink")) +
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
  print(p6)
  
  dev.off()
}

# Function to find and highlight significant results
find_significant_comparisons <- function(advantage_data, sig_threshold = 0.00625) {
  # Find significant results (p < threshold for either score or density)
  sig_rows <- advantage_data %>%
    filter(p_value_score < sig_threshold | p_value_density < sig_threshold)
  
  return(sig_rows)
}

# Function to create boxplots for significant results
create_significant_boxplots <- function(sig_results, data_dir, trait) {
  if (nrow(sig_results) == 0) {
    cat("No significant results found.\n")
    return(NULL)
  }
  
  cat("Found", nrow(sig_results), "significant results\n")
  
  # Create output directory
  sig_dir <- paste0(trait, "_significant_visualizations")
  dir.create(sig_dir, showWarnings = FALSE)
  
  # Create a plot for each significant finding
  plots <- list()
  
  for (i in 1:nrow(sig_results)) {
    row <- sig_results[i, ]
    n_value <- as.numeric(as.character(row$N))
    method <- as.character(row$method)
    
    # Determine which method is significant
    if (method == "Fix gene frequency and pathway size") {
      method_short <- "birewire"
    } else {
      method_short <- "keeppath"
    }
    
    # Try to find the matched data file
    matched_file <- list.files(path = data_dir, 
                              pattern = paste0(trait, ".*_n", n_value, "_", method_short, "_matched.csv"), 
                              full.names = TRUE)
    
    if (length(matched_file) == 0) {
      cat("Warning: Matched data file not found for", n_value, method_short, "\n")
      next
    }
    
    # Load matched data
    matched_data <- read.csv(matched_file[1])
    
    # Prepare plot data
    if (method_short == "birewire") {
      plot_data <- data.frame(
        group = ifelse(matched_data$in_birewire_top, "Random", "Control"),
        ot_score = matched_data$mean_score,
        method = "BireWire"
      )
      col_name <- "in_birewire_top"
    } else {
      plot_data <- data.frame(
        group = ifelse(matched_data$in_keeppath_top, "Random", "Control"),
        ot_score = matched_data$mean_score,
        method = "KeepPathSize"
      )
      col_name <- "in_keeppath_top"
    }
    
    # Calculate statistics for annotation
    if (method_short == "birewire") {
      random_scores <- matched_data$mean_score[matched_data$in_birewire_top]
      control_scores <- matched_data$mean_score[!matched_data$in_birewire_top]
    } else {
      random_scores <- matched_data$mean_score[matched_data$in_keeppath_top]
      control_scores <- matched_data$mean_score[!matched_data$in_keeppath_top]
    }
    
    t_test <- t.test(random_scores, control_scores)
    p_val <- t_test$p.value
    mean_random <- mean(random_scores, na.rm = TRUE)
    mean_control <- mean(control_scores, na.rm = TRUE)
    percent_diff <- ((mean_random - mean_control) / mean_control) * 100
    
    # Create plot
    p <- ggplot(plot_data, aes(x = group, y = ot_score, fill = group)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.4, size = 2) +
      scale_fill_manual(values = c("Control" = "grey80", "Random" = "#377EB8")) +
      labs(
        title = paste(trait, "-", tool_base, "- Top", n_value, "Pathways"),
        subtitle = paste0(
          method, "\n",
          "p = ", signif(p_val, 3), 
          " | Diff: ", signif(percent_diff, 2), "%"
        ),
        x = "",
        y = "OpenTargets Association Score"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9),
        axis.title.y = element_text(face = "bold")
      )
    
    # Add significance annotation
    sig_pval <- if(method_short == "birewire") row$p_value_score else row$p_value_density
    if (sig_pval < 0.001) {
      p <- p + annotate(
        "text", x = 1.5, y = max(plot_data$ot_score, na.rm = TRUE) * 1.05,
        label = "***", size = 8, fontface = "bold"
      )
    } else if (sig_pval < 0.01) {
      p <- p + annotate(
        "text", x = 1.5, y = max(plot_data$ot_score, na.rm = TRUE) * 1.05,
        label = "**", size = 8, fontface = "bold"
      )
    } else {
      p <- p + annotate(
        "text", x = 1.5, y = max(plot_data$ot_score, na.rm = TRUE) * 1.05,
        label = "*", size = 8, fontface = "bold"
      )
    }
    
    # Store the plot
    plots[[i]] <- p
    
    # Save individual plot
    file_name <- paste0(sig_dir, "/", trait, "_", tool_base, "_N", n_value, "_", method_short, "_boxplot.pdf")
    ggsave(file_name, p, width = 7, height = 6)
  }
  
  # Create combined figure with the most significant results (up to 6)
  if (length(plots) > 1) {
    # Sort by p-value and take up to 6
    sorted_results <- sig_results %>%
      arrange(pmin(p_value_score, p_value_density)) %>%
      head(min(6, nrow(sig_results)))
    
    sorted_plots <- list()
    for (i in 1:nrow(sorted_results)) {
      row_match <- sig_results$N == sorted_results$N[i] & 
                   sig_results$method == sorted_results$method[i]
      
      if (any(row_match)) {
        sorted_plots[[i]] <- plots[[which(row_match)[1]]]
      }
    }
    
    # Determine grid layout
    n_plots <- length(sorted_plots)
    if (n_plots <= 3) {
      ncols <- n_plots
    } else {
      ncols <- 3
    }
    
    # Create combined figure
    combined_file <- paste0(sig_dir, "/top_significant_results.pdf")
    pdf(combined_file, width = min(ncols * 5, 15), height = ceiling(n_plots/ncols) * 5)
    grid.arrange(grobs = sorted_plots, ncol = ncols)
    dev.off()
    cat("Created combined figure with top significant results:", combined_file, "\n")
  }
  
  # Create heatmap of significant results
  if (nrow(sig_results) > 0) {
    heat_data <- sig_results %>%
      mutate(
        sig_type = case_when(
          p_value_score < sig_threshold & p_value_density < sig_threshold ~ "Both",
          p_value_score < sig_threshold ~ "Score",
          p_value_density < sig_threshold ~ "Density",
          TRUE ~ "None"
        ),
        min_p = pmin(p_value_score, p_value_density)
      )
    
    p <- ggplot(heat_data, aes(
      x = N, 
      y = method,
      fill = -log10(min_p)
    )) +
      geom_tile(color = "white", size = 0.5) +
      scale_fill_gradient2(
        low = "white", mid = "orange", high = "red",
        midpoint = -log10(sig_threshold/2),
        name = "-log10(p)"
      ) +
      geom_text(
        aes(label = sig_type),
        size = 3
      ) +
      labs(
        title = paste(trait, tool_base, "- Significant Results Summary"),
        subtitle = paste("Threshold p <", sig_threshold),
        x = "Top N Pathways",
        y = "Method"
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(hjust = 1),
        panel.grid = element_blank()
      )
    
    # Save the heatmap
    ggsave(paste0(sig_dir, "/significant_results_heatmap.pdf"), p, width = 8, height = 4)
  }
  
  return(sig_dir)
}

# === Main Function ===
main <- function() {
  # Load advantage data
  cat("Loading advantage data from:", detailed_advantage_file, "\n")
  advantage_data <- read.csv(detailed_advantage_file)
  
  # Extract n_values from advantage data
  n_values <- unique(as.numeric(as.character(advantage_data$N)))
  cat("Found N values:", paste(n_values, collapse=", "), "\n")
  
  # Set method labels
  method_labels <- list(
    birewire = "Fix gene frequency and pathway size",
    keeppath = "Fix pathway size"
  )
  
  # 1. Create advantage plots
  pdf(paste0(trait, "_combined_advantage_plots.pdf"), width=10, height=8)
  
  # Mean score advantage
  p1 <- create_advantage_plot(
    advantage_data,
    y_col = "mean_score_advantage",
    title = "Mean Score Advantage (Top Pathways vs. Size-Matched Controls)",
    y_lab = "Score Difference"
  )
  print(p1)
  
  # Evidence density advantage
  p2 <- create_advantage_plot(
    advantage_data,
    y_col = "evidence_density_advantage",
    title = "Evidence Density Advantage (Top Pathways vs. Size-Matched Controls)",
    y_lab = "Density Difference"
  )
  print(p2)
  
  dev.off()
  
  # 2. Create combined boxplots from matched data files
  create_combined_boxplots(n_values, trait, tool_base, data_dir, method_labels)
  
  # 3. Find significant results and create special visualizations
  sig_results <- find_significant_comparisons(advantage_data, sig_threshold)
  if (!is.null(sig_results) && nrow(sig_results) > 0) {
    create_significant_boxplots(sig_results, data_dir, trait)
  }
  
  # 4. Create N-specific boxplots 
  for (n in n_values) {
    matched_files <- list.files(path = data_dir, 
                               pattern = paste0(trait, "_", tool_base, "_n", n, "_.*_matched.csv"),
                               full.names = TRUE)
    
    if (length(matched_files) >= 2) {
      pdf(paste0(trait, "_matching_plots_n", n, ".pdf"), width=10, height=8)
      
      # BireWire file
      bw_file <- grep("birewire", matched_files, value=TRUE)
      if (length(bw_file) > 0) {
        bw_data <- read.csv(bw_file[1])
        bw_plot_data <- bw_data %>%
          mutate(Group = factor(in_birewire_top, labels = c("Control", "BireWire")))
        
        # Mean score comparison
        p1 <- create_boxplot(
          bw_plot_data, 
          x_col = "Group", 
          y_col = "mean_score",
          fill_col = "Group",
          title = paste("OpenTargets Mean Score - BireWire Top", n, "vs Size-Matched Controls"),
          x_lab = "Group", 
          y_lab = "Mean OpenTargets Score",
          fill_values = c("Control" = "gray80", "BireWire" = "lightblue")
        )
        print(p1)
        
        # Evidence density comparison
        p2 <- create_boxplot(
          bw_plot_data, 
          x_col = "Group", 
          y_col = "evidence_density",
          fill_col = "Group",
          title = paste("Evidence Density - BireWire Top", n, "vs Size-Matched Controls"),
          x_lab = "Group", 
          y_lab = "Evidence Density",
          fill_values = c("Control" = "gray80", "BireWire" = "lightblue")
        )
        print(p2)
      }
      
      # KeepPathSize file
      kp_file <- grep("keeppath", matched_files, value=TRUE)
      if (length(kp_file) > 0) {
        kp_data <- read.csv(kp_file[1])
        kp_plot_data <- kp_data %>%
          mutate(Group = factor(in_keeppath_top, labels = c("Control", "KeepPathSize")))
        
        # Mean score comparison
        p3 <- create_boxplot(
          kp_plot_data, 
          x_col = "Group", 
          y_col = "mean_score",
          fill_col = "Group",
          title = paste("OpenTargets Mean Score - KeepPathSize Top", n, "vs Size-Matched Controls"),
          x_lab = "Group", 
          y_lab = "Mean OpenTargets Score",
          fill_values = c("Control" = "gray80", "KeepPathSize" = "lightpink")
        )
        print(p3)
        
        # Evidence density comparison
        p4 <- create_boxplot(
          kp_plot_data, 
          x_col = "Group", 
          y_col = "evidence_density",
          fill_col = "Group",
          title = paste("Evidence Density - KeepPathSize Top", n, "vs Size-Matched Controls"),
          x_lab = "Group", 
          y_lab = "Evidence Density",
          fill_values = c("Control" = "gray80", "KeepPathSize" = "lightpink")
        )
        print(p4)
      }
      
      dev.off()
    }
  }
  
  cat("All visualizations complete.\n")
}

# Execute the main function
main()