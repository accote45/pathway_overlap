process run_real_pascalx {
    executor 'lsf'
    clusterOptions '-P acc_paul_oreilly'
    tag "${trait}_pascalx_real"
    
    publishDir "${params.outdir}/pascalx_real/${trait}", mode: 'copy', overwrite: true
    
    input:
    tuple val(trait),
          path(standardized_gwas_file)
    
    output:
    tuple val(trait),
          path("${trait}_pascalx_real.txt")
    
    script:
    """
    python ${params.scripts_dir}/run_pascalx.py \\
        --gwas-file ${standardized_gwas_file} \\
        --trait ${trait} \\
        --geneset-file ${params.geneset_real} \\
        --output ${trait}_pascalx_real.txt
    """
}

process run_random_pascalx {
    executor 'lsf'
    clusterOptions '-P acc_paul_oreilly'
    tag "${trait}_pascalx_random${perm}_${rand_method}"
    
    publishDir "${params.outdir}/pascalx_random/${rand_method}/${params.background}/${trait}", mode: 'copy', overwrite: true
    
    input:
    tuple val(trait),
          path(standardized_gwas_file),
          val(rand_method),
          val(perm)
    
    output:
    tuple val(trait),
          path("${trait}_pascalx_random${perm}.${rand_method}.txt"),
          val(rand_method)
    
    script:
    def gmt_dir = params.gmt_dirs[rand_method]
    """
    python ${params.scripts_dir}/run_pascalx.py \\
        --gwas-file ${standardized_gwas_file} \\
        --trait ${trait} \\
        --geneset-file ${gmt_dir}/GeneSet.random${perm}.gmt \\
        --output ${trait}_pascalx_random${perm}.${rand_method}.txt
    """
}