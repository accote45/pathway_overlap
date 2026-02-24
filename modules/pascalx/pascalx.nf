process prepare_pascalx_gwas {
  executor 'local'
  tag "${trait}"
  
  input:
  tuple val(trait),
        path(gwas_file),
        val(rsid_col),
        val(pval_col),
        val(effect_allele),
        val(other_allele),
        val(beta_col)

  output:
  tuple val(trait),
        path("${trait}_pascalx.txt")

  script:
  """
  # Extract columns in order: rsid, a1, a2, pval, beta
  # This matches the expected column order for PascalX: rscol=0, a1col=1, a2col=2, pcol=3, bcol=4
  awk -v rsid="${rsid_col}" \\
      -v a1="${effect_allele}" \\
      -v a2="${other_allele}" \\
      -v pval="${pval_col}" \\
      -v beta="${beta_col}" \\
      'BEGIN {OFS="\\t"} 
       NR==1 {
          for (i=1; i<=NF; i++) {
              if (\$i == rsid) rsid_idx = i
              if (\$i == a1) a1_idx = i
              if (\$i == a2) a2_idx = i
              if (\$i == pval) pval_idx = i
              if (\$i == beta) beta_idx = i
          }
          print "RSID", "A1", "A2", "PVAL", "BETA"
       } 
       NR>1 {
          print \$(rsid_idx), \$(a1_idx), \$(a2_idx), \$(pval_idx), \$(beta_idx)
       }' ${gwas_file} > ${trait}_pascalx.txt
  """
}

process run_pascalx_genes {
  tag "${trait}"
  container "${params.pascalx_sif}"
  publishDir "${params.outdir}/pascalx_genes/${trait}", mode: 'copy', overwrite: true
  
  input:
  tuple val(trait),
        path(gwas_file)

  output:
  tuple val(trait),
        path("gene_scores.txt"),
        emit: gene_scores

  script:
  """
  python3 /scripts/tool_specific/pascalx/run_pascalx_genes.py \
    ${trait} \
    ${gwas_file} \
    ${params.pascalx_ref_panel} \
    ${params.pascalx_genome_annot}
  """
}

process run_real_pascalx {
  tag "${trait}"
  container "${params.pascalx_sif}"
  publishDir "${params.outdir}/pascalx_real/${trait}", mode: 'copy', overwrite: true
  
  input:
  tuple val(trait),
        path(gene_scores),
        val(rand_method)

  output:
  tuple val(trait),
        path("${trait}_real_pascalx.csv"),
        val(rand_method)

  script:
  """
  python3 /scripts/tool_specific/pascalx/run_pascalx_pathways.py \
    ${trait} \
    ${gene_scores} \
    ${params.gmt_file} \
    ${params.pascalx_genome_annot} \
    real
  """
}

process run_random_sets_pascalx {
  tag "${trait}_random${perm}_${rand_method}"
  container "${params.pascalx_sif}"
  publishDir "${params.outdir}/pascalx_random/${rand_method}/${params.background}/${trait}", mode: 'copy', overwrite: true
  
  input:
  tuple val(trait),
        path(gene_scores),
        val(rand_method),
        val(perm)
        
  output:
  tuple val(trait),
        path("${trait}_random${perm}.${rand_method}.csv"),
        val(rand_method)

  script:
  // Use randomized GMT directory pattern matching other tools
  def random_gmt = "${params.gmt_dirs[rand_method]}/GeneSet.random${perm}.gmt"
  
  """
  python3 /scripts/tool_specific/pascalx/run_pascalx_pathways.py \
    ${trait} \
    ${gene_scores} \
    ${random_gmt} \
    ${params.pascalx_genome_annot} \
    random${perm}.${rand_method}
  """
}
