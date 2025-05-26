process run_random_sets_prset {
  executor 'lsf'
  tag "${trait}_set_random${perm}"
    input:
    tuple val(trait),
          path(gwas_file),
          val(binary_trait),
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
    --binary-target ${binary_trait} \\
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
    --pheno-col ${trait} \\
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

