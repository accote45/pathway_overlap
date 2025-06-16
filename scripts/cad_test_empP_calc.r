library(data.table)
library(GSA)
library(tidyverse)
library(ggplot2)

cad <- read.table('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_real/cad/cad_real_set.gsa.out',header=T)

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
  trait_dirs <- trait_dirs[trait_dirs %like% "cad"]

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

# Load data for both randomization methods
calc_fpr("msigdbgenes", "keeppathsize", output_name="df_msigdbgenes_keeppathsize")
calc_fpr("msigdbgenes", "birewire", output_name="df_msigdbgenes_birewire")

# Combine results
df_msigdbgenes_keeppathsize$randomization <- "keeppathsize"
df_msigdbgenes_birewire$randomization <- "birewire"
df_combined <- rbind(df_msigdbgenes_keeppathsize, df_msigdbgenes_birewire)

calc_empirical_pvalues_fast <- function(real_results, background, random_method) {
  # Convert real results to data.table
  real_dt <- as.data.table(real_results)

  # Base path for random results
  base_path <- paste0('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/magma_random/',
                     random_method,'/',background,'/cad')
  
  if(!dir.exists(base_path)) {
    stop("Directory not found: ", base_path)
  }
  
  # Step 1: Read all random results at once
  message("Reading all random results files...")
  setwd(base_path)
  random_files <- list.files(pattern="*gsa.out")
  total_files <- length(random_files)
  message("Found ", total_files, " random results files")
  
  # Create a master table of all random results
  message("Creating master random results table...")
  all_random_results <- rbindlist(
    lapply(random_files, function(file) {
      tryCatch({
        dt <- fread(file, select = c("FULL_NAME", "P"))
        dt[, file_id := file]  # Add file identifier
        return(dt)
      }, error = function(e) {
        message("Error reading file ", file, ": ", e$message)
        return(NULL)
      })
    }),
    fill = TRUE,
    use.names = TRUE
  )
  
  # Keep only necessary columns from real data
  real_pathways_dt <- real_dt[, .(FULL_NAME, P)]
  
  # Step 2: Calculate empirical p-values directly
  message("Calculating empirical p-values...")
  
  # Count occurrences per pathway where random P ≤ real P
  counts <- real_pathways_dt[, {
    pathway <- FULL_NAME
    real_p <- P
    
    # Extract relevant random results for this pathway
    random_subset <- all_random_results[FULL_NAME == pathway]
    
    # Count how many are more significant than real
    more_significant_count <- sum(random_subset$P <= real_p, na.rm = TRUE)
    
    # Count unique files this pathway appears in
    unique_files <- uniqueN(random_subset$file_id)
    
    # Calculate empirical p-value
    list(
      real_pvalue = real_p,
      more_significant_count = more_significant_count,
      total_random = unique_files,
      empirical_pvalue = (more_significant_count + 1) / (unique_files + 1)
    )
  }, by = FULL_NAME]
  
  message("Done! Processed ", uniqueN(real_pathways_dt$FULL_NAME), " pathways.")
  
  return(as.data.frame(counts))
}


# Calculate empirical p-values for different randomization methods
empirical_keeppathsize <- calc_empirical_pvalues_fast(cad, "msigdbgenes", "keeppathsize")
empirical_birewire <- calc_empirical_pvalues_fast(cad, "msigdbgenes", "birewire")

# Add method information
empirical_keeppathsize$method <- "keeppathsize"
empirical_birewire$method <- "birewire"

# Combine results
empirical_combined <- rbind(empirical_keeppathsize, empirical_birewire)

# Save results
write.csv(empirical_combined, "cad_empirical_pvalues.csv", row.names = FALSE)






