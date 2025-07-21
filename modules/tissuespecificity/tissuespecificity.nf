process tissue_specificity_analysis {
    executor 'lsf'
    tag "${trait}_${tool_base}_tissue_specificity"
    
    input:
    tuple val(trait),
          val(tool_base),  // 'magma' or 'prset' without randomization method
          path(birewire_results),
          path(keeppathsize_results)
    
    output:
    tuple val(trait), 
          val(tool_base),
          path("${trait}_${tool_base}_size_matched_analysis_summary.tsv"),
          path("${trait}_${tool_base}_all_detailed_metrics.csv")
    
    publishDir "${params.outdir}/tissue_specificity/${tool_base}/${trait}", mode: 'copy', overwrite: true
    
    script:
    """
    module load R
    
    Rscript ${params.scripts_dir}/size_matched_tissuespecificity.r "${trait}" "${tool_base}" "${birewire_results}" "${keeppathsize_results}" "${params.geneset_real}" "${params.opentargets_n_values}" "${params.tissue_expression_data}"
    """
}