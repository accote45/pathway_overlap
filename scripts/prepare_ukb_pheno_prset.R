


# organize existing phenotype files
library(data.table)
library(tidyverse)

setwd('/sc/arion/projects/paul_oreilly/lab/shared/pheno')

cad <- read.csv('CAD.csv') %>% select(IID, CAD)
bmi <- read.csv('BMI.csv') %>% select(IID, BMI)
breast <- read.csv('BreastCancer.csv') %>% select(IID, BreastCancer)
ibd_firstdx <- read.csv('IBD.csv') %>% select(IID, IBD)
mdd <- read.csv('MDD.csv') %>% select(IID, MDD)
prostate <- read.csv('ProstateCancer.csv') %>% select(IID, ProstateCancer)
t2d <- read.csv('T2D.csv') %>% select(IID, T2D)

ad_cases <- read.table('./clive/Alzheimer_disease_case.txt',header=T)
all_iids <- unique(cad$IID) 
# Create AD dataframe with 0/1 encoding
ad <- data.frame(
  IID = all_iids,
  AD = ifelse(all_iids %in% ad_cases$sample_id, 1, 0)  
)

ibd <- ibd_firstdx %>%
  mutate(IBD = ifelse(is.na(IBD), 0, ifelse(IBD > 0, 1, 0))) %>%
  select(IID, IBD)

# List of all phenotype data frames
pheno_list <- list(cad, bmi, breast, ibd, mdd, prostate, t2d, ad)

# Merge all at once
combined_pheno <- pheno_list %>%
  reduce(full_join, by = "IID")

# read in file with lab values and basic covariates
pheno <- fread('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/ukb/ukb_phenofile_forgwas.txt')
pheno$BMI <- NULL

final_pheno <- pheno %>%
  left_join(combined_pheno, by = "IID")

# match column names to trait identifiers
colnames(final_pheno) <- c(
  "FID", "IID", "sample_id", "Total_bilirubin", "Alkaline_phosphatase", 
  "HDL_cholesterol", "IGF_1", "Lipoprotein_A", "Urea", "Vitamin_D", 
  "Eosinophill_percentage", "Lymphocyte_count", "Mean_platelet_thrombocyte_volume", 
  "Monocyte_percentage", "Platelet_crit", "FID.x", "Age", "Age2", "Sex", 
  "Centre", "SES", "FID.y", "Batch", "PC1", "PC2", "PC3", "PC4", "PC5", 
  "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", 
  "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "PC21", "PC22", "PC23", 
  "PC24", "PC25", "PC26", "PC27", "PC28", "PC29", "PC30", "PC31", "PC32", 
  "PC33", "PC34", "PC35", "PC36", "PC37", "PC38", "PC39", "PC40", "cad", 
  "bmi", "breast", "ibd", "mdd", "prostate", "t2d", "ad"
)

# Calculate sample sizes for continuous traits
trait_columns <- c("bmi")

sample_sizes <- sapply(trait_columns, function(trait) {
  sum(!is.na(final_pheno[[trait]]))
})

# For binary traits, also get case/control counts
binary_traits <- c("cad", "breast", "ibd", "mdd", "prostate", "t2d", "ad")

case_control_counts <- sapply(binary_traits, function(trait) {
  if(trait %in% colnames(final_pheno)) {
    cases <- sum(final_pheno[[trait]] == 1, na.rm = TRUE)
    controls <- sum(final_pheno[[trait]] == 0, na.rm = TRUE)
    total <- sum(!is.na(final_pheno[[trait]]))
    return(c(cases = cases, controls = controls, total = total))
  }
})

# Create summary table
trait_summary <- data.frame(
  trait = names(sample_sizes),
  sample_size = sample_sizes,
  cases = case_control_counts["cases", ],
  controls = case_control_counts["controls", ],
  case_percentage = round(case_control_counts["cases", ] / sample_sizes * 100, 2)
)

# Print summary
print(trait_summary)

write.table(final_pheno, file = '/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/ukb/ukb_phenofile_forprset.txt', sep = "\t", row.names = FALSE, quote = FALSE)







