nextflow.enable.dsl=2

// Process to convert random GMT files to GSA-MiXeR format
process convert_random_gmt_for_gsamixer {
  executor 'lsf'
  tag "${rand_method}_random${perm}"

  input:
  tuple val(trait),
        val(rand_method),
        val(perm),
        path(gmt_file),
        path(gtf_file)

  output:
  tuple val(trait),
        val(rand_method),
        val(perm),
        path("${rand_method}_random${perm}_baseline.txt"),
        path("${rand_method}_random${perm}_full_gene.txt"),
        path("${rand_method}_random${perm}_full_gene_set.txt")

  publishDir "${params.outdir}/gsamixer_random/${rand_method}/${trait}/reference", mode: 'copy', overwrite: true

  script:
  """
  module load R
  Rscript ${params.scripts_dir}/gsamixer/convert_random_gmt_for_gsamixer.R \
    ${gmt_file} \
    ${gtf_file} \
    ${rand_method}_random${perm}
  """
}

// Process to run GSA-MiXeR base model for random GMTs
process gsamixer_plsa_base_random {
  executor 'lsf'
  tag "${trait}_${rand_method}_random${perm}_base"

  input:
  tuple val(trait),
        val(rand_method),
        val(perm),
        path(baseline_txt),
        path(full_gene_txt),
        path(full_gene_set_txt),
        path(chrom_sumstats)

  output:
  tuple val(trait),
        val(rand_method),
        val(perm),
        path("${trait}_${rand_method}_random${perm}_base.json"),
        path("${trait}_${rand_method}_random${perm}_base.log")

  publishDir "${params.outdir}/gsamixer_random/${rand_method}/${trait}/random${perm}", mode: 'copy', overwrite: true

  script:
  """
  module load singularity
  ${params.mixer_py} plsa --gsa-base \
    --trait1-file ${trait}.chr@.sumstats.gz \
    --out ${trait}_${rand_method}_random${perm}_base \
    --bim-file ${params.mixer_ref_bim} \
    --use-complete-tag-indices \
    --loadlib-file ${params.mixer_ref_loadlib} \
    --go-file ${baseline_txt} \
    --annot-file ${params.mixer_ref_annot} \
    ${params.mixer_extra_flags:-}
  """
}

// Process to run GSA-MiXeR full model for random GMTs
process gsamixer_plsa_full_random {
  executor 'lsf'
  tag "${trait}_${rand_method}_random${perm}_full"

  input:
  tuple val(trait),
        val(rand_method),
        val(perm),
        path(baseline_txt),
        path(full_gene_txt),
        path(full_gene_set_txt),
        path(base_json),
        path(base_log)

  output:
  tuple val(trait),
        val(rand_method),
        val(perm),
        path("${trait}_${rand_method}_random${perm}_full.json"),
        path("${trait}_${rand_method}_random${perm}_full.log")

  publishDir "${params.outdir}/gsamixer_random/${rand_method}/${trait}/random${perm}", mode: 'copy', overwrite: true

  script:
  """
  module load singularity
  ${params.mixer_py} plsa --gsa-full \
    --trait1-file ${trait}.chr@.sumstats.gz \
    --out ${trait}_${rand_method}_random${perm}_full \
    --bim-file ${params.mixer_ref_bim} \
    --use-complete-tag-indices \
    --loadlib-file ${params.mixer_ref_loadlib} \
    --go-file ${full_gene_txt} \
    --go-file-test ${full_gene_set_txt} \
    --annot-file ${params.mixer_ref_annot} \
    --load-params-file ${base_json} \
    ${params.mixer_extra_flags:-}
  """
}