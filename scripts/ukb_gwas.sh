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


# merge gwas results for all chromosomes

for i in $(ls *linear);do
cat ${i} | grep -v CHROM | awk '{print $0"\t""173200"}'>> IGF-1_gwas.txt
done
CHROM	POS	ID	REF	ALT	PROVISIONAL_REF?	A1	OMITTED	A1_FREQ	TEST	OBS_CT	BETA	SE	T_STAT	P	ERRCODE N

cp *txt /sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/gwas