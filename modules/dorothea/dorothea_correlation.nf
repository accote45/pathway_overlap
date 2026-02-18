process dorothea_correlation {
    executor 'lsf'
    tag "${trait}_${tool_base}_dorothea_correlation"

    publishDir "${params.outdir}/dorothea_correlation/${trait}", mode: 'copy', overwrite: true

    input:
    tuple val(trait),
          val(tool_base),
          path(birewire_results),
          path(keeppathsize_results)

    output:
    tuple val(trait),
          val(tool_base),
          path("${trait}_${tool_base}_dorothea_rank_correlation_summary.csv")

    script:
    """
    module load R
    Rscript ${params.scripts_dir}/validation/dorothea/dorothea_correlation.R \\
      "${trait}" \\
      "${tool_base}" \\
      "${params.dorothea_scores_path}" \\
      "${birewire_results}" \\
      "${keeppathsize_results}"
    """
}