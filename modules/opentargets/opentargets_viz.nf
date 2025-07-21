process opentargets_visualization {
    executor 'lsf'
    tag "${trait}_${tool_base}_opentargets_viz"
    
    input:
    tuple val(trait),
          val(tool_base),
          path(detailed_advantage_file),
          path(birewire_matched),
          path(keeppath_matched)
    
    output:
    tuple val(trait), 
          val(tool_base),
          path("${trait}_combined_advantage_plots.pdf"),
          path("${trait}_combined_boxplots.pdf"),
          path("${trait}_matching_plots_n*.pdf"),
          path("${trait}_significant_visualizations/*")
    
    publishDir "${params.outdir}/opentargets_visualizations/${tool_base}/${trait}", mode: 'copy', overwrite: true
    
    script:
    """
    module load R
    
    # Create a directory for the data files
    mkdir -p data_files
    cp ${birewire_matched} ${keeppath_matched} data_files/
    
    # Run the visualization script
    Rscript ${params.scripts_dir}/size_matched_OT_viz.R \
        "${trait}" \
        "${tool_base}" \
        "${detailed_advantage_file}" \
        "${params.opentargets_sig_threshold}" \
        "./data_files"
    """
}