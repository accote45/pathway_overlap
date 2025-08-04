process opentargets_statistics {
    executor 'lsf'
    tag "${trait}_${tool_base}_opentargets_stats"
    
    publishDir "${params.outdir}/opentargets_stats/${trait}", mode: 'copy'

    input:
    tuple val(trait),
          val(tool_base),  // 'magma' or 'prset' without randomization method
          path(birewire_results),
          path(keeppathsize_results)
    
    output:
    tuple val(trait), 
          val(tool_base),
          path("${trait}_${tool_base}_detailed_advantage.csv"),
          path("*_n*_birewire_empp*.csv", optional: true),       // Changed from just birewire.csv
          path("*_n*_keeppathsize_empp*.csv", optional: true),       // Changed from just keeppathsize.csv
          path("*_n*_rawp*.csv", optional: true),
          path("*_n*_pvaluebeta*.csv", optional: true),          // New method
          path("*_n*_birewire_emppvalstdbeta*.csv", optional: true), // New method
          path("*_n*_keeppathsize_emppvalstdbeta*.csv", optional: true), // New method
          path("${trait}_${tool_base}_gene_disease_associations.csv"),
          path("${trait}_${tool_base}_advantage_summary*.csv")

    script:
    """
    # Run the modified Size-Matched OT statistics R script that supports all four ranking methods

ml R   

 Rscript ${params.scripts_dir}/size_matched_OT_stats_optimized.R "${trait}" "${tool_base}" "${birewire_results}" "${keeppathsize_results}" "${params.geneset_real}" "${params.opentargets_n_values}"
    """
}
