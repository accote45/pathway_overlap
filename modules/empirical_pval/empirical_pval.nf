process calc_empirical_pvalues {
    executor 'lsf'
    tag "${trait}_${tool}_empP"
    
    input:
    tuple val(trait),
          val(tool),
          path(real_results),
          val(random_dir)
    
    output:
    tuple val(trait),
          val(tool), 
          path("${trait}_${tool}_empirical_pvalues.txt")
    
    publishDir "${params.outdir}/empirical_pvalues/${tool}/${trait}", mode: 'copy', overwrite: true
    
    script:
    
    """
    module load R
    
    Rscript ${params.scripts_dir}/calc_empirical.R "${trait}" "${tool}" "${real_results}" "${random_dir}"
    """
}