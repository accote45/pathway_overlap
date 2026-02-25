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

process setup_pascalx_reference {
  executor 'local'
  publishDir "${params.outdir}/pascalx_ref_temp", mode: 'copy', overwrite: true
  
  output:
  path("EUR.1KG.GRCh37*"), emit: ref_files
  
  script:
  """
  # Copy reference panel files to writable location
  cp /sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/pascalx_reference/EUR.1KG.GRCh37* .
  """
}

process run_pascalx_genes {
  executor 'lsf'
  tag "${trait}"
  container "${params.pascalx_sif}"
  publishDir "${params.outdir}/pascalx_genes/${trait}", mode: 'copy', overwrite: true
  
  input:
  tuple val(trait),
        path(gwas_file),
        path('ref_files/*')  // Collect all reference files

  output:
  tuple val(trait),
        path("gene_scores.txt"),
        emit: gene_scores

  script:
  """
  python3 /scripts/tool_specific/pascalx/run_pascalx_genes.py \
    ${trait} \
    ${gwas_file} \
    /pascalx_ref/EUR.1KG.GRCh37 \
    ${params.pascalx_genome_annot}
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
  // Map to container-internal path
  def random_gmt = "/randomized_gene_sets/random_${rand_method}/GeneSet.random${perm}.gmt"
  
  """
  python3 /scripts/tool_specific/pascalx/run_pascalx_pathways.py \
    ${trait} \
    ${gene_scores} \
    ${random_gmt} \
    ${params.pascalx_genome_annot} \
    random${perm}.${rand_method}
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
  // Real GMT is already in /data mount
  def real_gmt = "/data/pathway_db/msigdb/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt"
  
  """
  python3 /scripts/tool_specific/pascalx/run_pascalx_pathways.py \
    ${trait} \
    ${gene_scores} \
    ${real_gmt} \
    ${params.pascalx_genome_annot} \
    real
  """
}
