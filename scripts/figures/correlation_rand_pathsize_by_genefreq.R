
## determine if genes in larger random pathways have higher GWAS association
library(data.table)
library(ggplot2)

# average pathway size per gene for random GMTs
gene_pathway_summary <- fread("/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/explore/gene_pathway_size_summary.csv")

# magma gene level results
magma_results <- fread("/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_real/height/height_real_set.gsa.genes.out")

gene_pathway_summary$GENE <- gene_pathway_summary$gene

# Merge by gene
master <- merge(gene_pathway_summary, magma_results, by = "GENE")

# Plot mean pathway size vs. z-statistic
pdf('test.pdf')
ggplot(master, aes(x = mean_size, y = ZSTAT)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  xlab("Mean pathway size (random GMTs)") +
  ylab("MAGMA gene z-statistic")
dev.off()

master$freq <- master$n/1000
pdf('test.pdf',width=7,height=7)
ggplot(master, aes(x = mean_size, y = freq)) +
  geom_point(alpha = 0.5) +
  xlab("Mean pathway size (random GMTs)") +
  ylab("Gene frequency") + theme_classic() + ylim(-1, NA)
dev.off()

# Correlation test
cor_test <- cor.test(master$freq, master$mean_size, method = "spearman")
print(cor_test)




