process delta_rank_ot_correlation {
    executor 'lsf'
    tag "${trait}_${tool_base}_delta_rank_ot"

    publishDir "${params.outdir}/delta_rank_ot/${trait}", mode: 'copy', overwrite: true

    input:
    tuple val(trait),
          val(tool_base),
          path(birewire_results)

    output:
    tuple val(trait),
          val(tool_base),
          path("${trait}_${tool_base}_delta_rank_with_OT_scores.csv"),
          path("${trait}_${tool_base}_delta_rank_ot_correlation_summary.csv"),
          path("figs_delta_rank_ot/*.pdf", optional: true)

    script:
    """
    module load R
    Rscript ${params.scripts_dir}/validation/delta_rank_ot_correlation.R \\
      "${trait}" \\
      "${tool_base}" \\
      "${birewire_results}" \\
      "${params.geneset_real}" \\
      "${params.opentargets_json_dir}" \\
      "${params.delta_rank_top_ns}"
    """
}