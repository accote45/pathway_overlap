process opentargets_visualization {
    executor 'lsf'
    tag "${trait}_${tool_base}_opentargets_viz"
    
    publishDir "${params.outdir}/opentargets_viz/${trait}", mode: 'copy'

    input:
    tuple val(trait), 
          val(tool_base),
          path(detailed_advantage),
          path("birewire_files/*"),
          path("keeppath_files/*"),
          path(gene_disease_associations),
          path("summary_files/*")
    
    output:
    tuple val(trait), 
          val(tool_base),
          path("${trait}_combined_advantage_plots*.pdf"),
          path("${trait}_combined_boxplots*.pdf"),
          path("${trait}_matching_plots_*.pdf"),
          path("${trait}_significant_visualizations/*"),
          path("${trait}_size_groups/*"),
          path("${trait}_ranking_comparison_*.pdf")  // New output for ranking comparisons
          path("*.csv"), optional: true

    script:
    """
    module load R
    
    # Create directories
    mkdir -p birewire_files keeppath_files rawp_files sigbeta_files summary_files ${trait}_significant_visualizations ${trait}_size_groups
    
    # Run visualization script with support for all ranking methods
    Rscript ${params.scripts_dir}/size_matched_OT_viz.R "${trait}" "${tool_base}" "${detailed_advantage}" "${params.opentargets_sig_threshold}" "."
    
    # Move any specific files to proper directories
    if [ -f ${trait}_*_bin_*.csv ]; then
      mv ${trait}_*_bin_*.csv ${trait}_size_groups/
    fi
    
    if [ -f ${trait}_all_top_pathways.csv ]; then
      mv ${trait}_all_top_pathways.csv ${trait}_size_groups/
    fi
    """
}