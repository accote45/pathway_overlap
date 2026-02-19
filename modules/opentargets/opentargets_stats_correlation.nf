process opentargets_stats_correlation {
    executor 'lsf'
    tag "${trait}_${tool_base}_correlation_stats"
    
    publishDir "${params.outdir}/opentargets_correlation/${trait}", mode: 'copy'

    input:
    tuple val(trait),
          val(tool_base),  // 'magma' or 'prset' without randomization method
          path(birewire_results),
          path(keeppathsize_results)
    
      output:
      tuple val(trait), 
              val(tool_base),
              path("${trait}_${tool_base}_rank_correlation_summary.csv"),
              path("${trait}_${tool_base}_gene_disease_associations.csv")

    script:
    """
    module load R
    Rscript ${params.scripts_dir}/validation/opentargets/OT_correlation_stats.R \\
      "${trait}" \\
      "${tool_base}" \\
      "${birewire_results}" \\
      "${keeppathsize_results}" \\
      "${params.geneset_real}" \\
      "${params.opentargets_n_values}"
    """
}