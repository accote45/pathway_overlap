library(data.table)
library(GSA)
library(tidyverse)
library(ggplot2)


# Read the gene set data
dat <- GSA.read.gmt('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/pathway_db/msigdb/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt')
path.list <- dat$genesets
names(path.list) <- dat$geneset.names

# Convert the path list to a data table
genes_long <- rbindlist(lapply(names(path.list), function(name) {
  data.table(value = path.list[[name]], name = name)
}))

# Calculate pathway sizes
pathsize <- as.data.frame(table(genes_long$name))
colnames(pathsize) <- c("Pathway", "Size")

# Create histogram of pathway sizes with improved x-axis
pdf('pathway_size_histogram.pdf', width = 8, height = 6)
ggplot(pathsize, aes(x = Freq)) +
  geom_histogram(binwidth = 10, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(x = 'Number of genes in pathway', 
       y = 'Count', 
       title = 'Distribution of pathway sizes in MSigDB') +
  theme_minimal() +
  # Use fewer breaks with better spacing
  scale_x_continuous(breaks = seq(0, max(pathsize$Freq), by = 200)) +
  geom_vline(xintercept = median(pathsize$Freq), linetype = "dashed", color = "red") +
  geom_text(aes(x = median(pathsize$Freq) + 20, 
                y = max(table(cut(pathsize$Freq, breaks = seq(0, max(pathsize$Freq) + 10, by = 10))))/2),
            label = paste("Median =", median(pathsize$Freq)), color = "red") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))  # Keep x-axis labels horizontal
dev.off()

