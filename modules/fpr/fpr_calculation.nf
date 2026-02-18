nextflow.enable.dsl=2

process calculate_fpr {
    executor 'lsf'
    tag "${trait}_${tool_base}_${rand_method}_fpr"
    
    publishDir "${params.outdir}/fpr_analysis/${tool_base}/${rand_method}/${trait}", mode: 'copy'

    input:
    tuple val(trait),
          val(tool_base),
          val(rand_method),
          path(random_files),
          val(random_dir)
    
    output:
    tuple val(trait),
          val(tool_base), 
          val(rand_method),
          path("${trait}_${tool_base}_${rand_method}_fpr_results.csv")

    script:
    """
    module load R
    
    Rscript ${params.scripts_dir}/core/calc_fpr.r "${trait}" "${tool_base}" "${rand_method}" "${random_dir}"
    """
}