process delta_rank_dorothea_correlation {
    executor 'lsf'
    tag "${trait}_${tool_base}_delta_rank_dorothea"

    publishDir "${params.outdir}/delta_rank_dorothea/${trait}", mode: 'copy', overwrite: true

    input:
    tuple val(trait),
          val(tool_base),
          path(birewire_results)

    output:
    tuple val(trait),
          val(tool_base),
          path("${trait}_${tool_base}_delta_rank_with_dorothea_scores.csv"),
          path("${trait}_${tool_base}_delta_rank_dorothea_correlation_summary.csv")

    script:
    """
    module load R
    Rscript ${params.scripts_dir}/delta_rank_dorothea_correlation.R \\
      "${trait}" \\
      "${tool_base}" \\
      "${birewire_results}" \\
      "${params.dorothea_scores_path}" \\
      "${params.delta_rank_top_ns}"
    """
}