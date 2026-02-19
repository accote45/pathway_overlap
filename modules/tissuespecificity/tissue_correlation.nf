process tissue_correlation_analysis {
    executor 'lsf'
    tag "${trait}_${tool_base}_tissue_correlation"
    
    publishDir "${params.outdir}/tissue_correlation/${trait}", mode: 'copy'

    input:
    tuple val(trait),
          val(tool_base),  // 'magma' or 'prset' without randomization method
          path(birewire_results),
          path(keeppathsize_results)
    
    output:
    tuple val(trait), 
          val(tool_base),
          path("${trait}_${tool_base}_tissue_correlation_summary.csv"),
          path("${trait}_${tool_base}_*_all_pathways_tissue_correlation_*.pdf"),
          path("${trait}_${tool_base}_*_top_500_pathways_tissue_correlation_*.pdf"),
          path("${trait}_${tool_base}_*_top_250_pathways_tissue_correlation_*.pdf"),
          path("${trait}_${tool_base}_*_top_100_pathways_tissue_correlation_*.pdf"),
          path("${trait}_${tool_base}_*_top_50_pathways_tissue_correlation_*.pdf"),
          path("${trait}_${tool_base}_best_method_by_tissue.csv")

    script:
    """
    module load R
    Rscript ${params.scripts_dir}/validation/tissue/tissue_correlation_stats.R \\
      "${trait}" \\
      "${tool_base}" \\
      "${params.tissue_expression_data}" \\
      "${birewire_results}" \\
      "${keeppathsize_results}" \\
      "${params.geneset_real}"
    """
}