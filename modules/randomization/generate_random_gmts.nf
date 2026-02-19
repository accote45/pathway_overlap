process generate_birewire_random_gmts {
  executor 'lsf'
  queue 'premium'
  cpus 1
  memory '15G'
  time '12h'
  clusterOptions '-P acc_paul_oreilly'
  
  input:
  path(input_gmt)
  val(num_random_sets)
  
  output:
  path("random_birewire/*.gmt"), emit: gmt_files
  val("${params.outdir}/randomized_gene_sets/random_birewire"), emit: birewire_dir
  
  publishDir "${params.outdir}/randomized_gene_sets", mode: 'copy', overwrite: false, pattern: "random_birewire/*.gmt"
  
  script:
  """
  module load R
  
  echo "========================================="
  echo "BiRewire GMT Generation Process"
  echo "========================================="
  echo "Generating ${num_random_sets} randomized GMT files"
  echo "This will take approximately 6-12 hours"
  echo "========================================="
  
  Rscript ${params.scripts_dir}/core/generate_birewire_gmts.R \\
    ${input_gmt} \\
    random_birewire \\
    ${num_random_sets}
  
  echo ""
  echo "BiRewire randomization complete!"
  echo "Generated ${num_random_sets} GMT files in random_birewire/"
  """
}

process generate_keeppathsize_random_gmts {
  executor 'lsf'
  queue 'premium'
  cpus 1
  memory '8G'
  time '4h'
  clusterOptions '-P acc_paul_oreilly'
  
  input:
  path(input_gmt)
  val(num_random_sets)
  
  output:
  path("random_keeppathsize/*.gmt"), emit: gmt_files
  val("${params.outdir}/randomized_gene_sets/random_keeppathsize"), emit: keeppathsize_dir
  
  publishDir "${params.outdir}/randomized_gene_sets", mode: 'copy', overwrite: false, pattern: "random_keeppathsize/*.gmt"
  
  script:
  """
  module load R
  
  echo "========================================="
  echo "KeepPathSize GMT Generation Process"
  echo "========================================="
  echo "Generating ${num_random_sets} randomized GMT files"
  echo "This will take approximately 2-4 hours"
  echo "========================================="
  
  Rscript ${params.scripts_dir}/core/generate_keeppathsize_gmts.R \\
    ${input_gmt} \\
    random_keeppathsize \\
    ${num_random_sets}
  
  echo ""
  echo "KeepPathSize randomization complete!"
  echo "Generated ${num_random_sets} GMT files in random_keeppathsize/"
  """
}