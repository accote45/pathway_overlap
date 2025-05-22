#!/bin/bash

# Set paths
BASE_PATH="/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/ukb"
OUTPUT_PATH="${BASE_PATH}/ukb_labvalue_gwas"
PHENO_FILE="${BASE_PATH}/ukb_phenofile_forgwas.txt"
PLINK2_PATH="/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/software/plink2"
GENOTYPE_FILE="${BASE_PATH}/ukb18177_eur_autosomes"
SCRIPT_DIR="/Users/cotea02/Desktop/pathway_overlap"

# Create output directories
mkdir -p ${OUTPUT_PATH}
mkdir -p ${OUTPUT_PATH}/jobs
mkdir -p ${OUTPUT_PATH}/logs
mkdir -p ${OUTPUT_PATH}/by_chr

ml R
# Create 70/30 train/test split
echo "Creating 70/30 train/test split..."
Rscript - <<EOF
library(data.table)
set.seed(42)

# Read phenotype file
pheno <- fread("${PHENO_FILE}")
total_samples <- nrow(pheno)

# Create 70% training set
train_idx <- sample(1:total_samples, size = round(0.7 * total_samples))
train_samples <- pheno[train_idx, c("FID", "IID")]
test_samples <- pheno[-train_idx, c("FID", "IID")]

# Write train/test sets
write.table(train_samples, "${BASE_PATH}/ukb_train_samples.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(test_samples, "${BASE_PATH}/ukb_test_samples.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

cat("Split complete: ", nrow(train_samples), " samples in training set, ", 
    nrow(test_samples), " samples in test set\n")
EOF

# First, submit jobs to split by chromosome
echo "Preparing chromosome-specific datasets..."

for chr in {1..22}; do
  job_file="${OUTPUT_PATH}/jobs/split_chr${chr}.lsf"
  
  cat > $job_file <<EOL
#!/bin/bash
#BSUB -L /bin/bash
#BSUB -J split_chr${chr}
#BSUB -q express
#BSUB -n 4
#BSUB -R "rusage[mem=20000]"
#BSUB -W 1:00
#BSUB -P acc_psychgen
#BSUB -o ${OUTPUT_PATH}/logs/split_chr${chr}.out
#BSUB -e ${OUTPUT_PATH}/logs/split_chr${chr}.err

# Extract chromosome data for training set
${PLINK2_PATH} \\
--bfile ${GENOTYPE_FILE} \\
--keep ${BASE_PATH}/ukb_train_samples.txt \\
--chr ${chr} \\
--make-bed \\
--out ${OUTPUT_PATH}/by_chr/chr${chr}

# Signal completion
touch ${OUTPUT_PATH}/by_chr/.chr${chr}_done
EOL

  # Submit chromosome splitting job
  bsub < $job_file
  echo "Submitted job to extract chromosome ${chr}"
done

# Read traits from CSV, removing quotes
echo "Reading traits from CSV file..."
traits=($(awk -F, 'NR>1 {gsub(/"/,"",$1); print $1}' /sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/heritability/selected_traits_for_gwas.csv))

echo "Found ${#traits[@]} traits to analyze"

# Create job submission script for each trait and chromosome
for trait in "${traits[@]}"; do
  echo "Creating jobs for trait: $trait"
  
  # Clean trait name for file naming
  clean_trait=$(echo $trait | tr ' ' '_')
  
  # Create trait directory
  mkdir -p ${OUTPUT_PATH}/${clean_trait}
  
  for chr in {1..22}; do
    job_file="${OUTPUT_PATH}/jobs/gwas_${clean_trait}_chr${chr}.lsf"
    
    cat > $job_file <<EOL
#!/bin/bash
#BSUB -L /bin/bash
#BSUB -J gwas_${clean_trait}_${chr}
#BSUB -q premium
#BSUB -n 2
#BSUB -R "rusage[mem=65000]"
#BSUB -W 4:00
#BSUB -P acc_psychgen
#BSUB -o ${OUTPUT_PATH}/logs/gwas_${clean_trait}_chr${chr}.out
#BSUB -e ${OUTPUT_PATH}/logs/gwas_${clean_trait}_chr${chr}.err
#BSUB -w "done(wait_split)"

# Run GWAS for this trait and chromosome
${PLINK2_PATH} \\
--bfile ${OUTPUT_PATH}/by_chr/chr${chr} \\
--glm hide-covar \\
--covar-variance-standardize \\
--covar ${PHENO_FILE} \\
--covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,Age,Sex,Centre \\
--threads 2 \\
--memory 50000 \\
--pheno ${PHENO_FILE} \\
--pheno-name "${trait}" \\
--no-input-missing-phenotype \\
--out ${OUTPUT_PATH}/${clean_trait}/${clean_trait}_chr${chr}_gwas

# Signal completion
touch ${OUTPUT_PATH}/${clean_trait}/.chr${chr}_completed
EOL

    # Submit the job with dependency on chromosome splitting
    bsub < $job_file
    
    echo "Job submitted for $trait chromosome ${chr}"
  done
done






#!/bin/bash

# Set paths
BASE_PATH="/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/ukb"
OUTPUT_PATH="${BASE_PATH}/ukb_labvalue_gwas"
PHENO_FILE="${BASE_PATH}/ukb_phenofile_forgwas.txt"

# Read traits from CSV, removing quotes
traits=($(awk -F, 'NR>1 {gsub(/"/,"",$1); print $1}' /sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/heritability/selected_traits_for_gwas.csv))

echo "Found ${#traits[@]} traits to process"

# Function to process a specific trait
process_trait() {
  local trait=$1
  local clean_trait=$(echo $trait | tr ' ' '_')
  local trait_dir="${OUTPUT_PATH}/${clean_trait}"
  local output_file="${trait_dir}/${clean_trait}_gwas.txt"
  
  echo "Processing ${clean_trait}..."
  
  # Create header for the merged file
  echo -e "CHROM\tPOS\tID\tREF\tALT\tPROVISIONAL_REF?\tA1\tOMITTED\tA1_FREQ\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP\tERRCODE\tN" > ${output_file}
  
  # Navigate to trait directory
  cd ${trait_dir}
  
  # Merge all chromosome files for this trait
  for linear_file in $(ls ${clean_trait}_chr*_gwas.${clean_trait}.glm.linear 2>/dev/null); do
    if [ -f "$linear_file" ]; then
      # Skip header and add sample size column
      grep -v "^#" $linear_file | grep -v "CHROM" | awk '{print $0"\t173200"}' >> ${output_file}
    fi
  done
  
  # Check if any data was merged
  if [ $(wc -l < ${output_file}) -le 1 ]; then
    echo "Warning: No data found for ${clean_trait}"
    return
  fi
  
  echo "Merged ${clean_trait} results to ${output_file}"
  
  # Add the A2 column using R
  Rscript - <<EOF
library(data.table)
# Read merged data
dat <- read.table("${output_file}",header=T,sep="\t")
dat <- as.data.table(dat)

# Create A2 column based on relationship between A1, REF, and ALT
dat[, A2 := ifelse(A1 == REF, ALT, 
             ifelse(A1 == ALT, REF, NA_character_))]

# Check for any rows where A2 is NA
na_count <- sum(is.na(dat\$A2))
if (na_count > 0) {
  cat("Warning:", na_count, "rows have NA values for A2 in ${clean_trait}\n")
}

# remove any rows with NA
dat <- dat[complete.cases(dat), ]

# Write updated data back to file
fwrite(dat, "${output_file}", sep = "\t")
cat("Added A2 column to ${clean_trait}\n")
EOF

  # Copy final file to the GWAS directory
  cp ${output_file} /sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/gwas/
  echo "Copied ${clean_trait} to GWAS directory"
  
  # Return to original directory
  cd ${OUTPUT_PATH}
}

# Process each trait in parallel (up to 4 at a time)
for trait in "${traits[@]}"; do
  process_trait "${trait}"
done


