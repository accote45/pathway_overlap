process gwas_remove_dup_snps {
  executor 'lsf'
  tag "${trait}_deduplicate_${rand_method}"
  
  input:
  tuple val(trait),
        path(gwas_file),
        val(binary_target),
        val(effect_allele),
        val(other_allele),
        val(rsid_col),
        val(pval_col),
        val(summary_statistic_name),
        val(summary_statistic_type),
        val(rand_method)
  
  output:
  tuple val(trait),
        path("${trait}_deduplicated.txt"),
        val(binary_target),
        val(effect_allele),
        val(other_allele),
        val(rsid_col),
        val(pval_col),
        val(summary_statistic_name),
        val(summary_statistic_type),
        val(rand_method)

  script:
  """
  #!/bin/bash
  module load R
  
  # Run R script
  Rscript - <<EOF
  # Load required libraries
  library(data.table)
  
  # Read the GWAS file
  cat("Reading GWAS file: ${gwas_file}\\n")
  gwas <- read.table("${gwas_file}", header=TRUE, na.strings=c("", "NA", "N/A", "."),check.names=FALSE, stringsAsFactors=FALSE)
  gwas <- as.data.table(gwas)
  
  # Replace empty cells with NA
  cat("Replacing empty cells with NA\\n")
  for(col in names(gwas)) {
    # Replace empty strings with NA
    if(is.character(gwas[[col]])) {
      empty_mask <- gwas[[col]] == "" | gwas[[col]] == " "
      if(sum(empty_mask, na.rm=TRUE) > 0) {
        cat("  - Found", sum(empty_mask, na.rm=TRUE), "empty cells in column", col, "\\n")
        gwas[empty_mask, (col) := NA]
      }
    }
  }
  
  # Report initial counts
  initial_snps <- nrow(gwas)
  cat("Initial number of SNPs:", initial_snps, "\\n")
  
  # Store column names as R variables
  rsid_col <- "${rsid_col}"
  
  # Remove rows where the SNP ID is NA or empty
  na_rsid <- is.na(gwas[[rsid_col]])
  if(sum(na_rsid) > 0) {
    cat("Removing", sum(na_rsid), "rows with NA values in", rsid_col, "column\\n")
    gwas <- gwas[!na_rsid]
  }
  
  # Identify duplicated SNP IDs
  cat("Checking for duplicated SNPs in column:", rsid_col, "\\n")
  
  # Count occurrences of each SNP ID
  snp_counts <- table(gwas[[rsid_col]])
  
  # Get IDs that appear more than once
  dup_ids <- names(snp_counts[snp_counts > 1])
  dup_count <- length(dup_ids)
  
  cat("Found", dup_count, "unique SNP IDs with duplicates\\n")
  
  if (dup_count > 0) {
    # Strategy: remove ALL instances of duplicate SNPs
    # Keep only SNPs that appear exactly once
    gwas_clean <- gwas[!(gwas[[rsid_col]] %in% dup_ids)]
    
    # Sort the data to maintain the original order where possible
    if ("CHROM" %in% colnames(gwas_clean) && "POS" %in% colnames(gwas_clean)) {
      setorder(gwas_clean, CHROM, POS)
    }
    
    # Report final counts
    final_snps <- nrow(gwas_clean)
    removed_snps <- initial_snps - final_snps
    cat("Final number of SNPs after deduplication:", final_snps, "\\n")
    cat("Removed", removed_snps, "SNPs with duplicate IDs\\n")
    
    # Write the deduplicated file
    fwrite(gwas_clean, "${trait}_deduplicated.txt", sep="\t", na="NA",quote=FALSE)
    cat("Wrote deduplicated file to ${trait}_deduplicated.txt\\n")
  } else {
    # Even if no duplicates, still write the file with NA values properly handled
    cat("No duplicates found. Writing file with properly handled NA values.\\n")
    fwrite(gwas, "${trait}_deduplicated.txt", sep="\t", na="NA",quote=FALSE)
    cat("Wrote file to ${trait}_deduplicated.txt\\n")
  }
  EOF
  """
}

process run_random_sets_prset {
  executor 'lsf'
  tag "${trait}_set_random${perm}_${rand_method}"
  
  input:
  tuple val(trait),
        path(gwas_file),
        val(binary_target),
        val(effect_allele),
        val(other_allele),
        val(rsid_col),
        val(pval_col),
        val(summary_statistic_name),
        val(summary_statistic_type),
        val(rand_method),
        val(perm)
  
  output:
  path("${trait}_set_random${perm}*")

  publishDir "${params.outdir}/prset_random/${rand_method}/${params.background}/${trait}", mode: 'copy', overwrite: true

  script:
  // Determine the correct GMT directory based on randomization method
  def gmt_base_dir = "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/randomized_gene_sets"
  def gmt_dir = rand_method == "birewire" ? "${gmt_base_dir}/random_birewire" : "${gmt_base_dir}/random_keeppathsize"
  
  """
  /sc/arion/projects/psychgen/cotea02_prset/PRSice_linux \\
    --a1 ${effect_allele} \\
    --a2 ${other_allele} \\
    --background /sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/msigdb.genes.txt:gene \\
    --bar-levels 1 \\
    --base ${gwas_file} \\
    --binary-target ${binary_target} \\
    --clump-kb 1000kb \\
    --clump-p 1.000000 \\
    --clump-r2 0.100000 \\
    --extract ${params.ukb_dir}/ukb18177-qc.snplist \\
    --fastscore \\
    --gtf /sc/arion/projects/paul_oreilly/lab/cotea02/project/data/reference/Homo_sapiens.GRCh37.75.gtf.gz \\
    --keep ${params.ukb_dir}/ukb_test_samples.txt \\
    --msigdb ${gmt_dir}/GeneSet.random${perm}.gmt \\
    --num-auto 22 \\
    --out ${trait}_set_random${perm} \\
    --pheno ${params.ukb_dir}/ukb_phenofile_forprset.txt \\
    --pheno-col ${trait}_resid \\
    --print-snp \\
    --pvalue ${pval_col} \\
    --set-perm 10000 \\
    --snp ${rsid_col} \\
    --stat ${summary_statistic_name} \\
    --${summary_statistic_type} \\
    --target ${params.ukb_dir}/ukb18177_chr1.22 \\
    --ultra \\
    --thread ${task.cpus} \\
    --wind-3 35kb \\
    --wind-5 35kb

  rm ${trait}_set_random${perm}.best
  """
}

