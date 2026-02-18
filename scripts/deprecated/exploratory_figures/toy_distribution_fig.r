# Load required libraries
library(ggplot2)

# Create sample data for the null distribution
set.seed(123)
null_scores <- rnorm(5000, mean = 0, sd = 0.8)

# Set the observed enrichment score (the vertical line)
observed_score <- 1.8

# Create dataframe for plotting
null_df <- data.frame(score = null_scores)

pdf('test.pdf')
# Create the plot with white background and teal accents
p <- ggplot() +
  # Add the null distribution histogram
  geom_histogram(data = null_df, aes(x = score, y = after_stat(density)), 
                 bins = 40, fill = "gray80", color = "gray50", alpha = 0.8) +
  # Add the observed score line
  geom_vline(xintercept = observed_score, color = "#008B8B", linewidth = 2) +
  # Add labels and title
  labs(
    title = "Enrichment Statistic Calibration",
    x = "Enrichment score",
    y = "Density under swap-null"
  ) +
  # Set the x axis limits
  scale_x_continuous(limits = c(-3.5, 3.5)) +
  # Add a clean theme with white background
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90"),
    text = element_text(color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    # Add a teal border around the plot like in the flowchart
    plot.margin = margin(20, 20, 20, 20),
    panel.border = element_rect(color = "#70C8C8", fill = NA, linewidth = 2)
  )

# Instead of adding a legend with segment and text, directly add an annotation
# for the observed score at the top of the plot
p <- p + 
  annotate("text", x = observed_score, y = 0.45, 
           label = "Observed pathway\nenrichment", 
           color = "black", hjust = 0.5, vjust = 0, fontface = "bold",
           size = 5)

# Make plot area smaller within the margins by adjusting the aspect ratio
p <- p + coord_cartesian(expand = FALSE) + 
         theme(aspect.ratio = 0.75)  # Make plot more compact vertically

print(p)
dev.off()