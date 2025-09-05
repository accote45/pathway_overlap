process malacards_correlation {
    executor 'lsf'
    tag "${trait}_${tool_base}_malacards_correlation"

    publishDir "${params.outdir}/malacards_correlation/${trait}", mode: 'copy', overwrite: true

    input:
    tuple val(trait),
          val(tool_base),          // 'magma' or 'prset' (base tool name)
          path(birewire_results),  // empirical p-values for birewire
          path(keeppathsize_results)

    output:
    tuple val(trait),
          val(tool_base),
          path("${trait}_${tool_base}_malacards_rank_correlation_summary.csv")

    script:
    """
    module load R
    Rscript ${params.scripts_dir}/malacards_correlation.R \
      "${trait}" \
      "${tool_base}" \
      "${params.malacards_path}" \
      "${birewire_results}" \
      "${keeppathsize_results}" \
      "${params.geneset_real}"
    """
}