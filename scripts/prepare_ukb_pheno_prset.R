


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


# Adjust outcomes for covariates, obtain residuals

# Define outcome variables (lab values + diseases)
continuous_outcomes <- c("Total_bilirubin", "Alkaline_phosphatase", "HDL_cholesterol", 
                  "IGF_1", "Lipoprotein_A", "Urea", "Vitamin_D", 
                  "Eosinophill_percentage", "Lymphocyte_count", 
                  "Mean_platelet_thrombocyte_volume", "Monocyte_percentage", 
                  "Platelet_crit", "bmi")

binary_outcomes <- c("cad", "breast", "ibd", "mdd", "prostate", "t2d", "ad")

# Define covariates
covariates <- c("Sex", "Age", "Batch", "Centre", 
                paste0("PC", 1:15))

# Create formula for covariates
covariate_formula <- paste(covariates, collapse = " + ")

# Function to get residuals from regression
get_residuals <- function(outcome, data) {
  # Create complete cases for this outcome
  complete_data <- data[!is.na(data[[outcome]]), ]
  
  # Create formula
  formula_str <- paste(outcome, "~", covariate_formula)
  
  # Fit model
  model <- lm(as.formula(formula_str), data = complete_data)
  
  # Get residuals and put them back in original data structure
  residuals_vector <- rep(NA, nrow(data))
  residuals_vector[!is.na(data[[outcome]])] <- residuals(model)
  
  return(residuals_vector)
}

# Apply regression adjustment to continuous outcomes
for(outcome in continuous_outcomes) {
  if(outcome %in% colnames(final_pheno)) {
    cat("Adjusting", outcome, "for covariates...\n")
    final_pheno[[paste0(outcome, "_resid")]] <- get_residuals(outcome, final_pheno)
  }
}

# For binary outcomes, use logistic regression
get_logistic_residuals <- function(outcome, data) {
  # Create complete cases for this outcome
  complete_data <- data[!is.na(data[[outcome]]), ]
  
  # Create formula
  formula_str <- paste(outcome, "~", covariate_formula)
  
  # Fit logistic model
  model <- glm(as.formula(formula_str), data = complete_data, family = binomial)
  
  # Get Pearson residuals
  residuals_vector <- rep(NA, nrow(data))
  residuals_vector[!is.na(data[[outcome]])] <- residuals(model, type = "pearson")
  
  return(residuals_vector)
}

# Apply logistic regression adjustment to binary outcomes
for(outcome in binary_outcomes) {
  if(outcome %in% colnames(final_pheno)) {
    cat("Adjusting", outcome, "for covariates using logistic regression...\n")
    final_pheno[[paste0(outcome, "_resid")]] <- get_logistic_residuals(outcome, final_pheno)
  }
}

# Save the adjusted phenotype file
write.table(final_pheno, 
            file = '/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/ukb/ukb_phenofile_forprset.txt', 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Print summary of adjustments
cat("\nAdjusted outcomes available:\n")
resid_cols <- colnames(final_pheno)[grepl("_resid$", colnames(final_pheno))]
cat(paste(resid_cols, collapse = "\n"))