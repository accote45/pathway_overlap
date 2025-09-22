process convert_gmt_for_gsamixer {
  executor 'lsf'
  tag "convert_gmt_for_gsamixer"

  input:
  path gmt_file
  path gtf_file

  output:
  path "baseline.txt"
  path "full_gene.txt"
  path "full_gene_set.txt"

  publishDir "${params.outdir}/gsamixer_reference", mode: 'copy', overwrite: true

  script:
  """
  module load R
  Rscript ${params.scripts_dir}/convert_geneset_gsamixer.R ${gmt_file} ${gtf_file}
  """
}

process prepare_gsamixer_sumstats {
  executor 'lsf'
  tag "${trait}_gsamixer_prepare"

  input:
  tuple val(trait),
        path(gwas_file),
        val(rsid_col),
        val(chr_col),
        val(pos_col),
        val(pval_col),
        val(n_col),
        val(binary_target),
        val(effect_allele),
        val(other_allele),
        val(summary_statistic_name),
        val(summary_statistic_type),
        val(se_col),
        val(neff_col)

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
  ml python
  ${params.mixer_py} split_sumstats \\
    --trait1-file ${sumstats_gz} \\
    --out ${trait}.chr\\@.sumstats.gz
  """
}

process gsamixer_plsa_base {
  executor 'lsf'
  tag "${trait}_gsamixer_base"

  input:
  tuple val(trait),
        path(chrom_sumstats),
        path(baseline_txt)

  output:
  tuple val(trait),
        path("${trait}_base.json"),
        path("${trait}_base.log")

  publishDir "${params.outdir}/gsamixer/${trait}", mode: 'copy', overwrite: true

  script:
  """
  module load singularity
  ml python
  ${params.mixer_py} plsa --gsa-base \
    --trait1-file ${trait}.chr\@.sumstats.gz \
    --out ${trait}_base \
    --bim-file ${params.mixer_ref_bim} \
    --loadlib-file ${params.mixer_ref_loadlib} \
    --go-file ${baseline_txt} \
    --annot-file ${params.mixer_ref_annot} \
    --go-extend-bp 35000 \
    ${params.mixer_extra_flags ?: ''}
  """
}

process gsamixer_plsa_full {
  executor 'lsf'
  tag "${trait}_gsamixer_full"

  input:
  tuple val(trait),
        path(base_json),
        path(base_log),
        path(full_gene_txt),
        path(full_gene_set_txt)

  output:
  tuple val(trait),
        path("${trait}_full.json"),
        path("${trait}_full.log")

  publishDir "${params.outdir}/gsamixer/${trait}", mode: 'copy', overwrite: true

  script:
  """
  module load singularity
  ml python
  ${params.mixer_py} plsa --gsa-full \\
    --trait1-file ${trait}.chr\\@.sumstats.gz \\
    --out ${trait}_full \\
    --bim-file ${params.mixer_ref_bim} \\
    --use-complete-tag-indices \\
    --loadlib-file ${params.mixer_ref_loadlib} \\
    --go-file ${full_gene_txt} \\
    --go-file-test ${full_gene_set_txt} \\
    --annot-file ${params.mixer_ref_annot} \\
    --go-extend-bp 35000 \\
    --load-params-file ${base_json} \\
    ${params.mixer_extra_flags ?: ''}
  """
}