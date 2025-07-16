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

# Process command line arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 5) {
  stop("Usage: Rscript size_matched_OT_analysis.R <trait> <tool_base> <birewire_results_file> <keeppathsize_results_file> <gmt_file> [top_n_values]")
}

trait <- args[1]
tool_base <- args[2]
birewire_results_file <- args[3]
keeppathsize_results_file <- args[4]
gmt_file <- args[5]

# Optional: parse comma-separated n_values or use default
n_values <- c(10, 20, 50, 100)  # Default
if(length(args) >= 6) {
  n_values <- as.numeric(unlist(strsplit(args[6], ",")))
}

# Define trait mapping
trait_mapping <- list(
  "t2d" = "MONDO_0005148",
  "cad" = "EFO_0001645", 
  "ad" = "MONDO_0004975",
  "mdd" = "MONDO_0002050",
  "scz" = "MONDO_0005090",
  "ibd" = "EFO_0000555",
  "breast" = "MONDO_0007254",
  "HDL_cholesterol" = "EFO_0004612",
  "Lymphocyte_count" = "EFO_0004587",
  "Platelet_crit" = "EFO_0007985",
  "Alkaline_phosphatase" = "EFO_0004533"
)

# === 1. Data Loading Functions ===

load_pathway_data <- function(gmt_file) {
  cat("Loading pathway data from", gmt_file, "...\n")
  dat <- GSA.read.gmt(gmt_file)
  path.list <- dat$genesets
  names(path.list) <- dat$geneset.names
  
  # Convert the path list to a data table
  genes_long <- rbindlist(lapply(names(path.list), function(name) {
    data.table(value = path.list[[name]], name = name)
  }))
  genes_long$targetId <- genes_long$value
  
  return(genes_long)
}

process_disease_data <- function(files, disease_id) {
  cat("Reading JSON files for disease", disease_id, "...\n")
  
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

load_magma_results <- function(real_results, emp_p_file) {
  cat("Reading MAGMA results...\n")
  
  # Read empirical p-values file
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
  
  return(list(
    birewire = birewire_results,
    keeppath = keeppath_results
  ))
}

calculate_pathway_scores <- function(genes_long, disease_targets) {
  cat("Calculating pathway-level scores...\n")
  
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
  
  return(pathway_scores)
}

# === 2. Matching and Analysis Functions ===

prepare_matching_data <- function(birewire_results, keeppath_results, pathway_scores, n) {
  cat(paste("\nPreparing data for top", n, "pathways...\n"))
  
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
  
  # Add OpenTargets evidence information
  matching_data <- merge(matching_data, pathway_scores, by="name", all.x=TRUE)
  
  return(matching_data)
}

perform_matching <- function(matching_data, target_col) {
  # Perform size matching for target column (either in_birewire_top or in_keeppath_top)
  tryCatch({
    formula <- as.formula(paste(target_col, "~ size"))
    m.out <- matchit(formula, data = matching_data, method = "nearest", ratio = 1)
    matched_data <- match.data(m.out)
    return(matched_data)
  }, error = function(e) {
    cat("Error in matching:", e$message, "\n")
    return(NULL)
  })
}

calculate_comparison_metrics <- function(matched_data, target_col, metrics=c("mean_score", "evidence_density", "max_score", "median_score")) {
  # Extract group of interest and control group based on target column
  target_group <- matched_data[matched_data[[target_col]] == TRUE, ]
  control_group <- matched_data[matched_data[[target_col]] == FALSE, ]
  
  # Compare metrics
  comparison <- data.frame()
  
  for(metric in metrics) {
    if(metric %in% colnames(matched_data)) {
      target_value <- mean(target_group[[metric]], na.rm=TRUE)
      control_value <- mean(control_group[[metric]], na.rm=TRUE)
      
      # Run t-test and extract all statistics
      t_test_result <- t.test(target_group[[metric]], control_group[[metric]])
      
      comparison <- rbind(comparison, data.frame(
        metric = metric,
        target_value = target_value,
        control_value = control_value,
        difference = target_value - control_value,
        percent_diff = 100 * (target_value - control_value) / ifelse(control_value == 0, 1, control_value),
        t_statistic = t_test_result$statistic,
        df = t_test_result$parameter,  # degrees of freedom
        p_value = t_test_result$p.value,
        significant = t_test_result$p.value < 0.05
      ))
    }
  }
  
  # Size verification
  target_size_avg <- mean(target_group$size)
  control_size_avg <- mean(control_group$size)
  size_diff_pct <- 100 * abs(target_size_avg - control_size_avg) / control_size_avg
  
  attr(comparison, "size_balance") <- list(
    target_size = target_size_avg,
    control_size = control_size_avg,
    difference_pct = size_diff_pct
  )
  
  return(comparison)
}

# === 3. Plotting Functions ===

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

# === 4. Individual Analysis Functions ===

analyze_method_at_n <- function(birewire_results, keeppath_results, pathway_scores, n, disease_name) {
  cat(paste("\nAnalyzing top", n, "pathways...\n"))
  
  # 1. Prepare data for matching
  matching_data <- prepare_matching_data(birewire_results, keeppath_results, pathway_scores, n)
  
  # 2. Perform matching for BireWire
  bw_matched <- perform_matching(matching_data, "in_birewire_top")
  if(is.null(bw_matched)) {
    cat("Matching failed for BireWire top", n, "pathways\n")
    return(NULL)
  }
  bw_matched$N <- n
  bw_matched$method <- "BireWire"
  
  # 3. Calculate metrics for BireWire
  bw_comparison <- calculate_comparison_metrics(bw_matched, "in_birewire_top")
  size_balance_bw <- attr(bw_comparison, "size_balance")
  
  # 4. Perform matching for KeepPathSize
  kp_matched <- perform_matching(matching_data, "in_keeppath_top")
  if(is.null(kp_matched)) {
    cat("Matching failed for KeepPathSize top", n, "pathways\n")
    return(NULL)
  }
  kp_matched$N <- n
  kp_matched$method <- "KeepPathSize"
  
  # 5. Calculate metrics for KeepPathSize
  kp_comparison <- calculate_comparison_metrics(kp_matched, "in_keeppath_top")
  size_balance_kp <- attr(kp_comparison, "size_balance")
  
  # 6. Display results
  cat("\nSize-matched comparison (BireWire top", n, "vs. size-matched pathways):\n")
  cat("- Average size in BireWire top:", round(size_balance_bw$target_size, 1), 
      "vs. matched pathways:", round(size_balance_bw$control_size, 1), 
      "(", round(size_balance_bw$difference_pct, 1), "% diff)\n")
  
  for(i in 1:nrow(bw_comparison)) {
    metric_name <- bw_comparison$metric[i]
    target_val <- bw_comparison$target_value[i]
    control_val <- bw_comparison$control_value[i]
    diff <- bw_comparison$difference[i]
    p_val <- bw_comparison$p_value[i]
    
    cat("- Average", metric_name, "in BireWire top", n, ":", round(target_val, 4), 
        "vs. matched pathways:", round(control_val, 4), 
        "(diff:", round(diff, 4), ",", 
        ifelse(diff > 0, "BireWire better", "Control better"),
        ", p =", format.pval(p_val, digits=3), ")\n")
  }
  
  cat("\nSize-matched comparison (KeepPathSize top", n, "vs. size-matched pathways):\n")
  cat("- Average size in KeepPathSize top:", round(size_balance_kp$target_size, 1), 
      "vs. matched pathways:", round(size_balance_kp$control_size, 1),
      "(", round(size_balance_kp$difference_pct, 1), "% diff)\n")
  
  for(i in 1:nrow(kp_comparison)) {
    metric_name <- kp_comparison$metric[i]
    target_val <- kp_comparison$target_value[i]
    control_val <- kp_comparison$control_value[i]
    diff <- kp_comparison$difference[i]
    p_val <- kp_comparison$p_value[i]
    
    cat("- Average", metric_name, "in KeepPathSize top", n, ":", round(target_val, 4), 
        "vs. matched pathways:", round(control_val, 4), 
        "(diff:", round(diff, 4), ",", 
        ifelse(diff > 0, "KeepPathSize better", "Control better"),
        ", p =", format.pval(p_val, digits=3), ")\n")
  }
  
  # 7. Create individual plots for this N value
  pdf(paste0(disease_name, "_matching_plots_n", n, ".pdf"), width=10, height=8)
  
  # Create BireWire plots
  bw_plot_data <- bw_matched %>%
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
  
  # Create KeepPathSize plots
  kp_plot_data <- kp_matched %>%
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
  
  dev.off()
  
  # 8. Return results for this N value
  return(list(
    birewire_matched = bw_matched,
    keeppath_matched = kp_matched,
    birewire_metrics = bw_comparison,
    keeppath_metrics = kp_comparison
  ))
}

calculate_method_advantages <- function(all_results, n_values, new_labels) {
  # Convert method names
  new_bw_label <- new_labels$birewire
  new_kp_label <- new_labels$keeppath
  
  advantage_data <- data.frame()
  
  for(n in n_values) {
    result_name <- paste0("n_", n)
    
    if(!is.null(all_results[[result_name]])) {
      bw_data <- all_results[[result_name]]$birewire_matched
      kp_data <- all_results[[result_name]]$keeppath_matched
      
      # Calculate BireWire advantage
      bw_top <- bw_data[bw_data$in_birewire_top == TRUE, ]
      bw_control <- bw_data[bw_data$in_birewire_top == FALSE, ]
      bw_advantage <- data.frame(
        N = factor(n, levels=n_values),
        method = new_bw_label,
        mean_score_top = mean(bw_top$mean_score, na.rm=TRUE),
        mean_score_control = mean(bw_control$mean_score, na.rm=TRUE),
        mean_score_advantage = mean(bw_top$mean_score, na.rm=TRUE) - mean(bw_control$mean_score, na.rm=TRUE),
        evidence_density_top = mean(bw_top$evidence_density, na.rm=TRUE),
        evidence_density_control = mean(bw_control$evidence_density, na.rm=TRUE),
        evidence_density_advantage = mean(bw_top$evidence_density, na.rm=TRUE) - mean(bw_control$evidence_density, na.rm=TRUE),
        p_value_score = all_results[[result_name]]$birewire_metrics$p_value[all_results[[result_name]]$birewire_metrics$metric == "mean_score"],
        t_statistic_score = all_results[[result_name]]$birewire_metrics$t_statistic[all_results[[result_name]]$birewire_metrics$metric == "mean_score"],
        df_score = all_results[[result_name]]$birewire_metrics$df[all_results[[result_name]]$birewire_metrics$metric == "mean_score"],
        p_value_density = all_results[[result_name]]$birewire_metrics$p_value[all_results[[result_name]]$birewire_metrics$metric == "evidence_density"],
        t_statistic_density = all_results[[result_name]]$birewire_metrics$t_statistic[all_results[[result_name]]$birewire_metrics$metric == "evidence_density"],
        df_density = all_results[[result_name]]$birewire_metrics$df[all_results[[result_name]]$birewire_metrics$metric == "evidence_density"]
      )
      
      # Calculate KeepPathSize advantage
      kp_top <- kp_data[kp_data$in_keeppath_top == TRUE, ]
      kp_control <- kp_data[kp_data$in_keeppath_top == FALSE, ]
      kp_advantage <- data.frame(
        N = factor(n, levels=n_values),
        method = new_kp_label,
        mean_score_top = mean(kp_top$mean_score, na.rm=TRUE),
        mean_score_control = mean(kp_control$mean_score, na.rm=TRUE),
        mean_score_advantage = mean(kp_top$mean_score, na.rm=TRUE) - mean(kp_control$mean_score, na.rm=TRUE),
        evidence_density_top = mean(kp_top$evidence_density, na.rm=TRUE),
        evidence_density_control = mean(kp_control$evidence_density, na.rm=TRUE),
        evidence_density_advantage = mean(kp_top$evidence_density, na.rm=TRUE) - mean(kp_control$evidence_density, na.rm=TRUE),
        p_value_score = all_results[[result_name]]$keeppath_metrics$p_value[all_results[[result_name]]$keeppath_metrics$metric == "mean_score"],
        t_statistic_score = all_results[[result_name]]$keeppath_metrics$t_statistic[all_results[[result_name]]$keeppath_metrics$metric == "mean_score"],
        df_score = all_results[[result_name]]$keeppath_metrics$df[all_results[[result_name]]$keeppath_metrics$metric == "mean_score"],
        p_value_density = all_results[[result_name]]$keeppath_metrics$p_value[all_results[[result_name]]$keeppath_metrics$metric == "evidence_density"],
        t_statistic_density = all_results[[result_name]]$keeppath_metrics$t_statistic[all_results[[result_name]]$keeppath_metrics$metric == "evidence_density"],
        df_density = all_results[[result_name]]$keeppath_metrics$df[all_results[[result_name]]$keeppath_metrics$metric == "evidence_density"]
      )
      
      advantage_data <- rbind(advantage_data, bw_advantage, kp_advantage)
    }
  }
  
  return(advantage_data)
}

create_combined_visualizations <- function(advantage_data, disease_name) {
  # Add significance indicators
  advantage_data$mean_score_sig <- ifelse(advantage_data$p_value_score < 0.05, "*", "")
  advantage_data$evidence_density_sig <- ifelse(advantage_data$p_value_density < 0.05, "*", "")
  
  # Create combined PDF with all advantage plots
  pdf(paste0(disease_name, "_combined_advantage_plots.pdf"), width=10, height=8)
  
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
  
  # Create summary table
  summary_data <- advantage_data %>%
    select(N, method, 
           mean_score_advantage, t_statistic_score, df_score, mean_score_sig, 
           evidence_density_advantage, t_statistic_density, df_density, evidence_density_sig) %>%
    pivot_wider(
      names_from = method,
      values_from = c(mean_score_advantage, t_statistic_score, df_score, mean_score_sig, 
                     evidence_density_advantage, t_statistic_density, df_density, evidence_density_sig)
    )
  
  # Write summary data to CSV
  write.csv(summary_data, paste0(disease_name, "_advantage_summary.csv"), row.names=FALSE)
  
  # Write detailed advantage data to CSV
  write.csv(advantage_data, paste0(disease_name, "_detailed_advantage.csv"), row.names=FALSE)
  
  return(list(summary = summary_data, detailed = advantage_data))
}

create_combined_boxplots <- function(all_results, disease_name, method_labels) {
  # Extract and combine data from all N values
  combined_bw_data <- data.frame()
  combined_kp_data <- data.frame()
  
  for(result_name in names(all_results)) {
    if(!is.null(all_results[[result_name]])) {
      n <- as.numeric(gsub("n_", "", result_name))
      
      # Process BireWire data
      bw_data <- all_results[[result_name]]$birewire_matched
      bw_data$N <- paste0("Top ", n)
      bw_data$Group <- factor(bw_data$in_birewire_top, 
                             labels = c("Control", "BireWire"))
      combined_bw_data <- rbind(combined_bw_data, bw_data)
      
      # Process KeepPathSize data
      kp_data <- all_results[[result_name]]$keeppath_matched
      kp_data$N <- paste0("Top ", n)
      kp_data$Group <- factor(kp_data$in_keeppath_top, 
                             labels = c("Control", "KeepPathSize"))
      combined_kp_data <- rbind(combined_kp_data, kp_data)
    }
  }
  
  # Order N factor correctly
  combined_bw_data$N <- factor(combined_bw_data$N, 
                              levels = paste0("Top ", c(10, 20, 50, 100)))
  combined_kp_data$N <- factor(combined_kp_data$N, 
                              levels = paste0("Top ", c(10, 20, 50, 100)))
  
  # Increase text sizes for all plots
  title_size <- 18
  axis_title_size <- 16
  axis_text_size <- 14
  legend_title_size <- 14
  legend_text_size <- 14
  
  # Create BireWire plots
  pdf(paste0(disease_name, "_combined_boxplots.pdf"), width=12, height=10)
  
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
  
  return(list(
    birewire_plots = list(mean_score=p1, evidence_density=p2),
    keeppath_plots = list(mean_score=p3, evidence_density=p4),
    comparison_plots = list(mean_score=p5, evidence_density=p6)
  ))
}


# === Main Execution ===
main <- function() {
  # Get disease ID from mapping
  disease_id <- trait_mapping[[trait]]
  if(is.null(disease_id)) {
    stop(paste("No mapping found for trait:", trait))
  }
  
  cat("======= Starting Size-Matched Analysis for", trait, "=======\n")
  cat("Using disease ID:", disease_id, "\n")
  
  # 1. Load data
  genes_long <- load_pathway_data(gmt_file)
  
  # Process OpenTargets JSON files for disease
  files.temp <- list.files(path = "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/results/drugtarget_test/associationByDatatypeDirect", pattern = "*.json",full.names=TRUE)
  files <- c(files.temp,list.files(path = "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/results/drugtarget_test/associationByDatatypeDirect", pattern = "*.json",full.names = TRUE))
  if(length(files) == 0) {
    stop("No JSON files found in current directory")
  }
  
  disease_data <- process_disease_data(files, disease_id)
  
  # Select max evidence count for each target
  disease_targets <- disease_data %>%
    group_by(targetId) %>%
    slice_max(evidenceCount, with_ties = FALSE) %>%
    ungroup()
  
  # Load results files based on tool type
  birewire_data <- read.table(paste0('./results/empirical_pvalues/magma_birewire/cad/',birewire_results_file))
  keeppath_data <- read.table(paste0('./results/empirical_pvalues/magma_keeppath/cad/',keeppathsize_results_file))

  cat("Loaded", nrow(birewire_data), "pathways from BireWire results\n")
  cat("Loaded", nrow(keeppath_data), "pathways from KeepPathSize results\n")
  
  # Calculate pathway scores
  pathway_scores <- calculate_pathway_scores(genes_long, disease_targets)
  
  # 2. Run analysis for each N value
  all_results <- list()
  
  for(n in n_values) {
    result_key <- paste0("n_", n)
    all_results[[result_key]] <- analyze_method_at_n(
      birewire_data, 
      keeppath_data, 
      pathway_scores, 
      n, 
      trait
    )
  }
  
  # 3. Create advantage calculations and visualizations
  method_labels <- list(
    birewire = "Fix gene frequency and pathway size",
    keeppath = "Fix pathway size"
  )
  
  advantage_data <- calculate_method_advantages(all_results, n_values, method_labels)
  visualization_results <- create_combined_visualizations(advantage_data, trait)
  combined_boxplot_results <- create_combined_boxplots(all_results, trait, method_labels)
  
  # 4. Create a summary file for Nextflow tracking
  summary_file <- paste0(tolower(trait), "_", tolower(tool_base), "_size_matched_analysis_summary.tsv")
  write.table(
    data.frame(
      trait = trait,
      disease_id = disease_id,
      tool_base = tool_base,
      n_values = paste(n_values, collapse = ","),
      num_birewire_pathways = nrow(birewire_data),
      num_keeppath_pathways = nrow(keeppath_data)
    ),
    file = summary_file,
    sep = "\t", 
    row.names = FALSE, 
    quote = FALSE
  )
  
  cat("Summary written to", summary_file, "\n")
  cat("======= Analysis Complete =======\n")
}

# Execute the main function
main()