process prepare_input	{
  executor 'local'
  input:
  		tuple	val(trait),
  				path(gwas_file),
  				val(rsid_col),
	  			val(chr_col),
	  			val(pos_col),
	  			val(pval_col),
	  			val(n_col)

  output:
  		tuple	val(trait),
  				path(gwas_file),
  				val(rsid_col),
  				val(pval_col),
  				val(n_col),
  				path("snp_loc.txt"), emit: snp_loc_data

  script:
  """
awk -v rsid="${rsid_col}" -v chr="${chr_col}" -v pos="${pos_col}" 'NR==1 {for (i=1; i<=NF; i++) {f[\$i] = i}} NR>1 {print \$(f[rsid]), \$(f[chr]), \$(f[pos])}' ${gwas_file} > snp_loc.txt
  """
}



process annotate_genes	{
  executor 'local'
  input:
  		tuple	val(trait),
  				path(gwas_file),
  				val(rsid_col),
  				val(pval_col),
  				val(n_col),
  				path(snp_file),
  				path(gene_file)

  output:
  		tuple	val(trait),
  				path(gwas_file),
  				val(rsid_col),
  				val(pval_col),
  				val(n_col),
  				path("*.genes.annot"), emit: gene_annot_data

  script:
  """
  module load magma_gwas/1.10

  magma --annotate window=35,35 \\
        --snp-loc ${snp_file} \\
        --gene-loc ${gene_file} \\
        --out ${trait}_case.control
  """
}



process run_gene_analysis	{
  input:
  		tuple	val(trait),
  				path(gwas_file),
  				val(rsid_col),
  				val(pval_col),
  				val(n_col),
  				path(gene_annot)

  output:
  		tuple	val(trait),
  				path("${trait}_case.control_genebased.genes.raw"), emit: gene_results

  publishDir { "${params.outdir}/magma_real/${trait}" }, mode: 'copy'

  script:
  """
  module load magma_gwas/1.10
  
  magma \\
    --bfile ${params.bfile} \\
    --pval ${gwas_file} use=${rsid_col},${pval_col} ncol=${n_col} \\
    --gene-annot ${gene_annot} \\
    --out ${trait}_case.control_genebased
  """
}


process run_real_geneset	{
  executor 'lsf'
  input:
  		tuple	val(trait),
  				path(gene_results_file)

  output:
    path("${trait}_real_set.*")

  publishDir { "${params.outdir}/magma_real/${trait}" }, mode: 'copy'

  script:
  """
  module load magma_gwas/1.10

  magma \\
    --gene-results ${gene_results_file} \\
    --set-annot ${params.geneset_real} \\
    --out ${trait}_real_set
  """
}


process run_random_sets {
  executor 'lsf'
  tag "${trait}_set_random${perm}"
  input:
    tuple val(trait),
          path(gene_results_file),
	  val(perm)
  output:
    path("${trait}_set_random${perm}*")

  publishDir "${params.outdir}/magma_random/${params.randomization_method}/${params.background}/${trait}", mode: 'copy'

  script:
  """
  module load magma_gwas/1.10

  magma \\
    --gene-results ${gene_results_file} \\
    --set-annot /sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/randomized_gene_sets/random_${params.randomization_method}/GeneSet.random${perm}.gmt \\
    --out ${trait}_set_random${perm}
  """

}

