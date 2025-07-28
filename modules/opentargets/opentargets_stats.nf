process opentargets_statistics {
    executor 'lsf'
    tag "${trait}_${tool_base}_opentargets_stats"
    
    input:
    tuple val(trait),
          val(tool_base),  // 'magma' or 'prset' without randomization method
          path(birewire_results),
          path(keeppathsize_results)
    
    output:
    tuple val(trait), 
          val(tool_base),
          path("${trait}_detailed_advantage.csv"),
          path("*_n*_birewire*matched.csv"),
          path("*_n*_keeppath*matched.csv"),
          path("${trait}_${tool_base}_gene_disease_associations.csv"),
          path("${trait}_advantage_summary*.csv")  // Include all ranking method outputs

    publishDir "${params.outdir}/size_matched_analysis/${tool_base}/${trait}/data", mode: 'copy', overwrite: true
    
    script:
    """
    module load R
    
    # Run the Size-Matched OT statistics R script
    Rscript ${params.scripts_dir}/size_matched_OT_stats.R "${trait}" "${tool_base}" "${birewire_results}" "${keeppathsize_results}" "${params.geneset_real}" "${params.opentargets_n_values}"
    """
}