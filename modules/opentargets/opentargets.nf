process opentargets_comparison {
    executor 'lsf'
    tag "${trait}_${tool_base}_size_matched"
    
    input:
    tuple val(trait),
          val(tool_base),  // 'magma' or 'prset' without randomization method
          path(birewire_results),
          path(keeppathsize_results)
    
    output:
    tuple val(trait), 
          val(tool_base),
          path("${trait}_${tool_base}_size_matched_analysis_summary.tsv"),
          path("${trait}_combined_advantage_plots.pdf"),
          path("${trait}_combined_boxplots.pdf"),
          path("${trait}_advantage_summary.csv"),
          path("${trait}_detailed_advantage.csv"),
          path("${trait}_matching_plots_n*.pdf")

    
    publishDir "${params.outdir}/size_matched_analysis/${tool_base}/${trait}", mode: 'copy', overwrite: true
    
    script:
    """
    module load R
    
    # Run the Size-Matched OT analysis R script
    Rscript ${params.scripts_dir}/size_matched_OT_analysis.R "${trait}" "${tool_base}" "${birewire_results}" "${keeppathsize_results}" "${params.geneset_real}" "10,20,50,100"
    """
}