process gwas_remove_dup_snps {
  executor 'lsf'
  tag "${trait}_deduplicate"
  
  input:
    tuple val(trait),
          path(gwas_file),
          val(binary_target),
          val(effect_allele),
          val(other_allele),
          val(rsid_col),
          val(pval_col),
          val(summary_statistic_name),
          val(summary_statistic_type)
  
  output:
    tuple val(trait),
          path("${trait}_deduplicated.txt"),
          val(binary_target),
          val(effect_allele),
          val(other_allele),
          val(rsid_col),
          val(pval_col),
          val(summary_statistic_name),
          val(summary_statistic_type)

  script:
  """
  #!/usr/bin/env Rscript
  
  # Load required libraries
  library(data.table)
  
  # Read the GWAS file
  cat("Reading GWAS file: ${gwas_file}\\n")
  gwas <- fread("${gwas_file}", header=TRUE)
  
  # Report initial counts
  initial_snps <- nrow(gwas)
  cat("Initial number of SNPs:", initial_snps, "\\n")
  
  # Identify duplicated SNPs
  cat("Checking for duplicated SNPs in column: ${rsid_col}\\n")
  dup_snps <- duplicated(gwas[[${rsid_col}]]) | duplicated(gwas[[${rsid_col}]], fromLast=TRUE)
  dup_count <- sum(dup_snps)
  cat("Found", dup_count, "SNPs in duplicated groups\\n")
  
  if (dup_count > 0) {
    # Get the list of duplicated SNP IDs
    dup_ids <- unique(gwas[[${rsid_col}]][dup_snps])
    cat("Number of unique SNP IDs with duplicates:", length(dup_ids), "\\n")
    
    # Strategy: remove duplicate SNPs    
    # Create a clean dataset
    gwas_clean <- data.table()
    
    # Process non-duplicated SNPs
    gwas_clean <- rbind(gwas_clean, gwas[!gwas[[${rsid_col}]] %in% dup_ids])
    
    # Sort the data to maintain the original order where possible
    if ("CHROM" %in% colnames(gwas_clean) && "POS" %in% colnames(gwas_clean)) {
      setorder(gwas_clean, CHROM, POS)
    }
    
    # Report final counts
    final_snps <- nrow(gwas_clean)
    cat("Final number of SNPs after deduplication:", final_snps, "\\n")
    cat("Removed", initial_snps - final_snps, "duplicate SNPs\\n")
    
    # Write the deduplicated file
    fwrite(gwas_clean, "${trait}_deduplicated.txt", sep="\t")
    cat("Wrote deduplicated file to ${trait}_deduplicated.txt\\n")
  } else {
    cat("No duplicates found. Creating symlink to original file.\\n")
    # If no duplicates, just create a symlink to the original file
    system("ln -s ${gwas_file} ${trait}_deduplicated.txt")
  }
  """
}



process run_random_sets_prset {
  executor 'lsf'
  tag "${trait}_set_random${perm}"
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
          val(perm)
  
  output:
    path("${trait}_set_random${perm}*")

  publishDir "${params.outdir}/prset_random/${params.randomization_method}/${params.background}/${trait}", mode: 'copy'

  script:
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
    --feature exon,gene,protein_coding,CDS \\
    --gtf /sc/arion/projects/paul_oreilly/lab/cotea02/project/data/reference/Homo_sapiens.GRCh37.75.gtf.gz \\
    --keep ${params.ukb_dir}/ukb_test_samples.txt \\
    --msigdb ${params.gmt_dir}/GeneSet.random${perm}.gmt \\
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
    --target ${params.ukb_dir}/ukb_labvalue_gwas/by_chr/chr# \\
    --thread 48 \\
    --ultra \\
    --wind-3 35kb \\
    --wind-5 35kb
  """
}

