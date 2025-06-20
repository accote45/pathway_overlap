process opentargets_comparison {
    executor 'lsf'
    tag "${trait}_${tool_base}_targets"
    
    input:
    tuple val(trait),
          val(tool_base),  // 'magma' or 'prset' without randomization method
          path(birewire_results),
          path(keeppathsize_results)
    
    output:
    tuple val(trait), 
          val(tool_base),
          path("${trait}_${tool_base}_opentargets_comparison.tsv")
    
    publishDir "${params.outdir}/opentargets_comparison/${tool_base}/${trait}", mode: 'copy', overwrite: true
    
    script:
    """
    module load R
    
    # Run the OpenTargets comparison R script
    Rscript ${params.scripts_dir}/opentargets.R "${trait}" "${tool_base}" "${birewire_results}" "${keeppathsize_results}"
    """
}