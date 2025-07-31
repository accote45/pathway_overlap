process opentargets_visualization {
    executor 'lsf'
    tag "${trait}_${tool_base}_opentargets_viz"
    
    publishDir "${params.outdir}/opentargets_viz/${trait}", mode: 'copy'

    input:
    tuple val(trait), 
          val(tool_base),
          path(detailed_advantage),
          path(gene_disease_associations)
    
    output:
    tuple val(trait), 
          val(tool_base),
          path("${trait}_combined_advantage_plots.pdf"),
          path("${trait}_combined_boxplots.pdf")

    script:
    """
    module load R
    
    # Create directories
    mkdir -p birewire_files keeppath_files rawp_files sigbeta_files summary_files ${trait}_significant_visualizations ${trait}_size_groups

    # get directory where OT stats output is stored
    DATA_DIR="${params.outdir}/size_matched_analysis/${tool_base}/${trait}/data/"
    
    # Run visualization script with support for all ranking methods
    Rscript ${params.scripts_dir}/size_matched_OT_viz.R "${trait}" "${tool_base}" "${detailed_advantage}" "${params.opentargets_sig_threshold}" "\$DATA_DIR"
    """
}