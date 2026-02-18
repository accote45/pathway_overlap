process delta_rank_malacards_correlation {
    executor 'lsf'
    tag "${trait}_${tool_base}_delta_rank_malacards"

    publishDir "${params.outdir}/delta_rank_malacards/${trait}", mode: 'copy', overwrite: true

    input:
    tuple val(trait),
          val(tool_base),
          path(birewire_results)

    output:
    tuple val(trait),
          val(tool_base),
          path("${trait}_${tool_base}_delta_rank_with_MC_scores.csv"),
          path("${trait}_${tool_base}_delta_rank_malacards_correlation_summary.csv"),
          path("figs_delta_rank_malacards/*.pdf", optional: true)

    script:
    """
    module load R
    Rscript ${params.scripts_dir}/validation/delta_rank_malacards_correlation.r \\
      "${trait}" \\
      "${tool_base}" \\
      "${birewire_results}" \\
      "${params.geneset_real}" \\
      "${params.malacards_path}" \\
      "${params.delta_rank_top_ns}"
    """
}