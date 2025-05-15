library(tidyverse)
library(data.table)
library(ggplot2)
library(patchwork)

read_files <- function(path) {
  setwd(path)
  files <- list.files(pattern="*gsa.out")
  ldf <- lapply(files, function(f) tryCatch(read.table(f, header=T), error=function(e) NULL))
  ldf <- ldf[!sapply(ldf, is.null)]  # Remove any NULL entries
  return(ldf)
}

calc_fpr <- function(background, random, output_name=NULL) {
  # Base path for all traits
  base_path <- paste0('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_random/',random,'/',background)
  
  # Get trait directories directly from specified path
  trait_dirs <- list.dirs(base_path, full.names=FALSE, recursive=FALSE)
  trait_dirs <- trait_dirs[trait_dirs != ""] # Filter out empty names
  
  # Create paths list from available directories
  paths <- setNames(
    lapply(trait_dirs, function(trait_dir) {
      paste0(base_path, '/', trait_dir)
    }),
    trait_dirs
  )
  
  message(paste("Processing", length(paths), "traits"))
  
  # Process each path
  results <- lapply(names(paths), function(name) {
    tryCatch({
      message(paste("Processing", name, "at", paths[[name]]))
      if(dir.exists(paths[[name]])) {
        ldf <- read_files(paths[[name]])
        if(length(ldf) > 0) {
          all.ldf <- as.data.frame(do.call(rbind, ldf))
          nom <- as.data.frame(table(all.ldf[all.ldf$P < 0.05, ]$FULL_NAME))
          nom$FPR <- nom$Freq/length(ldf)
          nom$dataset <- name
          return(nom)
        } else {
          message(paste("No .gsa.out files found in:", paths[[name]]))
          return(NULL)
        }
      } else {
        message(paste("Directory does not exist:", paths[[name]]))
        return(NULL)
      }
    }, error = function(e) {
      message(paste("Error processing", name, ":", e$message))
      return(NULL)
    })
  })
  
  # Filter out NULL results
  results <- results[!sapply(results, is.null)]
  
  if(length(results) == 0) {
    message("No valid results found. Check paths and file availability.")
    return(NULL)
  }
  
  master <- bind_rows(results)

  if (!is.null(output_name)) {
    assign(output_name, master, envir = .GlobalEnv)
  }
  
  return(master)
}

# Run calculations with all traits from the specified directory
calc_fpr("msigdbgenes", "birewire", output_name="df_msigdbgenes_birewire")

##### master boxplot

# Create a list of data frames with associated labels
dfs <- list(
  df_msigdbgenes_birewire = list(rand = "Keep gene frequency and pathway size", background = "Pathway database genes")
)

# Apply labels programmatically
for (dfname in names(dfs)) {
  if (exists(dfname)) {
    assign(
      dfname,
      within(get(dfname), {
        rand <- dfs[[dfname]]$rand
        background <- dfs[[dfname]]$background
      }))
  }
}

# Combine all data frames
datlist <- mget(ls(pattern = "^df_"))
options(stringsAsFactors=FALSE)
master <- rbindlist(datlist)

# Optimize plot layout based on number of traits
n_traits <- length(unique(master$dataset))
height_factor <- max(10, n_traits * 0.5)  # Adjust height based on trait count

# First, let's recreate the column distribution in a different way
# Get unique datasets and distribute them into columns
unique_datasets <- unique(master$dataset)
n_datasets <- length(unique_datasets)
n_columns <- 3
n_per_column <- ceiling(n_datasets / n_columns)

# Create a mapping - which column each dataset should appear in
dataset_column <- data.frame(
  dataset = unique_datasets,
  column = rep(1:n_columns, each=n_per_column)[1:n_datasets]
)

# Merge this mapping back into master data
master <- merge(master, dataset_column, by="dataset")

# Now create the plot with facet_wrap
pdf('trait_comparison_fixed.pdf', width=18, height=12)
ggplot(master, aes(x=background, y=FPR, fill=rand)) + 
  geom_boxplot() + 
  theme_bw() +
  stat_summary(fun=mean, geom="point", shape=18, size=2.5, color="darkred", position=position_dodge(0.75)) +
  geom_hline(yintercept=0.05, linetype="dashed") +
  ylab("FPR") + 
  xlab("Background") +
  # Use facet_wrap with nested vars
  facet_wrap(~ column + dataset, ncol=3, scales="free_y") +  
  theme(strip.text = element_text(size=10),
        strip.background = element_rect(fill="lightgray"),
        panel.spacing = unit(1, "lines"))
dev.off()

###### correlation btw pathway size and FPR
dat <- read.table("/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/results/pathwaydb_enrichment_OLD/msigdball_pathway_size.txt")
dat <- dat %>% rownames_to_column('Var1')
master <- merge(master,dat,by="Var1")



library(ggpubr)
library(patchwork)

# Create correlation plot with stat_cor
pdf('pathway_correlation.pdf', height=20, width=20)
ggplot(master, aes(x=pathway_size, y=FPR)) + 
  geom_point(shape=1) + 
  geom_smooth(method="lm", se=FALSE, linetype="dashed", alpha=0.3) +
  theme_classic() +
  geom_hline(yintercept=0.05, linetype="dashed") + 
  ylab("FPR") + 
  xlab("Pathway size") +
  # Add correlation statistics
  stat_cor(
    method = "pearson",
    label.x.npc = 0.85,  # Position at right side
    label.y.npc = 0.3,    # Position at top
    size = 5,               # Text size
    cor.coef.name = "r",    # Label for coefficient
    p.accuracy = 0.001,     # P-value decimal places
    r.accuracy = 0.01,
    label.sep = "\n"       # Correlation decimal places
  ) +
  facet_wrap(~ dataset, ncol=3, scales="free_y") +  
  theme(strip.text = element_text(size=15),
        strip.background = element_rect(fill="lightgray"),
        panel.spacing = unit(1, "lines"))
dev.off()
