// ---------------------------------------------------------------------------
// Setup stages: reference preparation, the one-time MAGMA annotation, and the
// per-replicate architecture + summary-statistic simulation.
// ---------------------------------------------------------------------------

// Reference prep and annotation depend only on the gene-location file and the
// LD panel, so they run ONCE and are shared by every condition and replicate.
process prepare_reference {
  publishDir "${params.outdir}/reference", mode: 'copy', overwrite: true,
             pattern: "ref_{gene_pool.txt,reference_diag.txt,snp_loc.txt}"

  output:
  path "ref_reference.rds", emit: reference
  path "ref_gene_pool.txt", emit: gene_pool
  path "ref_snp_loc.txt",   emit: snp_loc
  path "ref_reference_diag.txt"

  stub:
  """
  touch ref_reference.rds ref_gene_pool.txt ref_snp_loc.txt ref_reference_diag.txt
  """

  script:
  """
  ${params.load_r}
  export SIM_SCRIPTS=${params.sim_scripts}

  for f in "${params.gene_file}" "${params.bfile}.bim" "${params.bfile}.bed" "${params.bfile}.fam"; do
    [ -s "\$f" ] || { echo "MISSING REQUIRED ASSET: \$f" >&2; exit 1; }
  done

  Rscript ${params.sim_scripts}/02_prepare_reference.R \\
    gene_loc=${params.gene_file} \\
    bim=${params.bfile}.bim \\
    out_prefix=ref \\
    min_snps=${params.min_snps_per_gene}
  """
}

process magma_annotate {
  publishDir "${params.outdir}/reference", mode: 'copy', overwrite: true

  input:
  path snp_loc

  output:
  path "sim_annot.genes.annot", emit: annot

  stub:
  """
  touch sim_annot.genes.annot
  """

  script:
  """
  ${params.load_magma}

  magma --annotate window=35,35 \\
        --snp-loc ${snp_loc} \\
        --gene-loc ${params.gene_file} \\
        --out sim_annot
  """
}

process build_architecture {
  tag "${cond}_rep${rep}"
  publishDir "${params.outdir}/architectures/${cond}/${rep}", mode: 'copy', overwrite: true,
             pattern: "a_{ground_truth.tsv,architecture_diag.txt}"

  input:
  tuple val(cond), val(rep), val(overlap_frac), val(hub_signal),
        val(mu_low), val(mu_med), val(mu_high), path(gene_pool)

  output:
  tuple val(cond), val(rep), path("a_geneset.gmt"),      emit: gmt
  tuple val(cond), val(rep), path("a_gene_signal.tsv"),  emit: signal
  tuple val(cond), val(rep), path("a_ground_truth.tsv"), emit: truth
  path "a_architecture_diag.txt"

  stub:
  """
  touch a_geneset.gmt a_gene_signal.tsv a_ground_truth.tsv a_architecture_diag.txt
  """

  script:
  // Seed is deterministic in (condition, replicate) so any replicate can be
  // reproduced in isolation.
  def seed = (cond.hashCode() & 0xffff) * 100000 + (rep as int)
  """
  ${params.load_r}
  export SIM_SCRIPTS=${params.sim_scripts}

  Rscript ${params.sim_scripts}/01_build_architecture.R \\
    gene_pool_file=${gene_pool} \\
    out_prefix=a \\
    seed=${seed} \\
    n_pathways=${params.n_pathways} \\
    pathway_size=${params.pathway_size} \\
    variable_sizes=${params.variable_sizes} \\
    n_enriched_per_tier=${params.n_enriched_per_tier} \\
    n_contaminable_nulls=${params.n_contaminable_nulls} \\
    overlap_frac=${overlap_frac} \\
    hub_degree=${params.hub_degree} \\
    hub_slots_in_high=${params.hub_slots_in_high} \\
    hub_signal=${hub_signal} \\
    mu_low=${mu_low} \\
    mu_med=${mu_med} \\
    mu_high=${mu_high} \\
    mu_hub=${params.mu_hub} \\
    effect_scale=${params.effect_scale} \\
    frac_causal=${params.frac_causal}
  """
}

process simulate_sumstats {
  tag "${cond}_rep${rep}"
  publishDir "${params.outdir}/sumstats/${cond}/${rep}", mode: 'copy', overwrite: true,
             pattern: "s_{sumstats_diag.txt,causal_snps.tsv}"

  input:
  tuple val(cond), val(rep), path(gene_signal), path(reference)

  output:
  tuple val(cond), val(rep), path("s_sumstats.txt"), emit: sumstats
  tuple val(cond), val(rep), path("s_causal_snps.tsv"), emit: causal
  path "s_sumstats_diag.txt"

  stub:
  """
  touch s_sumstats.txt s_causal_snps.tsv s_sumstats_diag.txt
  """

  script:
  def seed = (cond.hashCode() & 0xffff) * 100000 + (rep as int) + 7
  """
  ${params.load_r}
  export SIM_SCRIPTS=${params.sim_scripts}

  Rscript ${params.sim_scripts}/03_simulate_sumstats.R \\
    reference_rds=${reference} \\
    gene_signal=${gene_signal} \\
    bfile=${params.bfile} \\
    out_prefix=s \\
    seed=${seed} \\
    max_block_snps=${params.max_block_snps} \\
    ld_lambda=${params.ld_lambda} \\
    independent=${params.independent_z}
  """
}
