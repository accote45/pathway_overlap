nextflow.enable.dsl=2

process prepare_gsamixer_sumstats {
  executor 'lsf'
  tag "${trait}_gsamixer_prepare"

  input:
  tuple val(trait),
        path(gwas_file)

  output:
  tuple val(trait),
        path("${trait}.sumstats.gz")

  publishDir "${params.outdir}/gsamixer/${trait}", mode: 'copy', overwrite: true

  script:
  """
  module load R
  Rscript ${params.scripts_dir}/gsamixer/prepare_sumstats.R \\
    --trait "${trait}" \\
    --traits-config "${params.traits_config}" \\
    --input "${gwas_file}" \\
    --output "${trait}.sumstats.gz"
  """
}

process split_gsamixer_sumstats {
  executor 'lsf'
  tag "${trait}_gsamixer_split"

  input:
  tuple val(trait),
        path(sumstats_gz)

  output:
  tuple val(trait),
        path("${trait}.chr*.sumstats.gz")

  publishDir "${params.outdir}/gsamixer/${trait}", mode: 'copy', overwrite: true

  script:
  """
  module load singularity
  ${params.mixer_py} split_sumstats \\
    --trait1-file ${sumstats_gz} \\
    --out ${trait}.chr@.sumstats.gz
  """
}

process gsamixer_plsa_base {
  executor 'lsf'
  tag "${trait}_gsamixer_base"

  input:
  tuple val(trait),
        path(chrom_sumstats)

  output:
  tuple val(trait),
        path("${trait}_base.json"),
        path("${trait}_base.log")

  publishDir "${params.outdir}/gsamixer/${trait}", mode: 'copy', overwrite: true

  script:
  """
  module load singularity
  ${params.mixer_py} plsa --gsa-base \\
    --trait1-file ${trait}.chr@.sumstats.gz \\
    --out ${trait}_base \\
    --bim-file ${params.mixer_ref_bim} \\
    --use-complete-tag-indices \\
    --loadlib-file ${params.mixer_ref_loadlib} \\
    --go-file ${params.mixer_go_base} \\
    --annot-file ${params.mixer_ref_annot} \\
    ${params.mixer_extra_flags:-}
  """
}

process gsamixer_plsa_full {
  executor 'lsf'
  tag "${trait}_gsamixer_full"

  input:
  tuple val(trait),
        path(base_json),
        path(base_log)

  output:
  tuple val(trait),
        path("${trait}_full.json"),
        path("${trait}_full.log")

  publishDir "${params.outdir}/gsamixer/${trait}", mode: 'copy', overwrite: true

  script:
  """
  module load singularity
  ${params.mixer_py} plsa --gsa-full \\
    --trait1-file ${trait}.chr@.sumstats.gz \\
    --out ${trait}_full \\
    --bim-file ${params.mixer_ref_bim} \\
    --use-complete-tag-indices \\
    --loadlib-file ${params.mixer_ref_loadlib} \\
    --go-file ${params.mixer_go_full} \\
    --go-file-test ${params.mixer_go_full_test} \\
    --annot-file ${params.mixer_ref_annot} \\
    --load-params-file ${base_json} \\
    ${params.mixer_extra_flags:-}
  """
}