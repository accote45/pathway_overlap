// ---------------------------------------------------------------------------
// Enrichment stages. Every statistical step here is the PARENT pipeline's own
// code: generate_birewire_gmts.R, generate_keeppathsize_gmts.R and
// calc_empirical.r are called from ${params.core_scripts} unmodified. The only
// simulation-specific machinery is the batching of MAGMA gene-set calls, which
// is verified to be result-identical by `verify_batching`.
// ---------------------------------------------------------------------------

process magma_gene_analysis {
  tag "${cond}_rep${rep}"
  publishDir "${params.outdir}/magma_genes/${cond}/${rep}", mode: 'copy', overwrite: true,
             pattern: "*.genes.out"

  input:
  tuple val(cond), val(rep), path(sumstats), path(annot)

  output:
  tuple val(cond), val(rep), path("g.genes.raw"), emit: raw
  path "g.genes.out"

  stub:
  """
  touch g.genes.raw g.genes.out
  """

  script:
  """
  ${params.load_magma}

  magma --bfile ${params.bfile} \\
        --pval ${sumstats} use=SNP,P N=${params.gwas_n} \\
        --gene-annot ${annot} \\
        --out g
  """
}

// One process for both null models; only the parent script differs. The
// simulated matrix is ~100 x 50, so this takes minutes rather than the
// 6-12 h the parent pipeline needs on the 5878-pathway MSigDB database.
process randomize_gmts {
  tag "${cond}_rep${rep}_${method}"

  input:
  tuple val(cond), val(rep), path(gmt), val(method)

  output:
  tuple val(cond), val(rep), val(method), path("random_gmts"), emit: gmt_dir

  stub:
  """
  mkdir -p random_gmts
  for i in \$(seq 1 ${params.num_random_sets}); do touch random_gmts/GeneSet.random\$i.gmt; done
  """

  script:
  if (method == 'birewire')
    """
    ${params.load_r}
    Rscript ${params.core_scripts}/generate_birewire_gmts.R \\
      ${gmt} random_gmts ${params.num_random_sets}
    """
  else if (method == 'keeppathsize')
    """
    ${params.load_r}
    Rscript ${params.core_scripts}/generate_keeppathsize_gmts.R \\
      ${gmt} random_gmts ${params.num_random_sets} ${task.cpus}
    """
  else
    error "Unknown randomization method: ${method}"
}

process magma_geneset_real {
  tag "${cond}_rep${rep}"
  publishDir "${params.outdir}/magma_real/${cond}/${rep}", mode: 'copy', overwrite: true

  input:
  tuple val(cond), val(rep), path(gene_raw), path(gmt)

  output:
  tuple val(cond), val(rep), path("real_set.gsa.out"), emit: real_gsa

  stub:
  """
  touch real_set.gsa.out
  """

  script:
  """
  ${params.load_magma}
  ${params.load_r}
  export SIM_SCRIPTS=${params.sim_scripts}

  magma --gene-results ${gene_raw} --set-annot ${gmt} --out raw_real_set

  # Guarantee a FULL_NAME column: MAGMA only emits one when set names are long
  # enough to be truncated, and the simulated names are short. calc_empirical.r
  # reads FULL_NAME for MAGMA, so this keeps that script usable unchanged.
  Rscript ${params.sim_scripts}/05_split_batch_gsa.R \\
    gsa_out=raw_real_set.gsa.out out_dir=. strip_prefix=false
  mv raw_real_set.gsa.out real_set.gsa.out
  """
}

process magma_geneset_random_batch {
  tag "${cond}_rep${rep}_${method}_b${batch}"

  input:
  tuple val(cond), val(rep), val(method), path(gmt_dir), path(gene_raw),
        val(batch), val(first), val(last)

  output:
  tuple val(cond), val(rep), val(method), path("split/*.gsa.out"), emit: split

  stub:
  """
  mkdir -p split
  for i in \$(seq ${first} ${last}); do touch split/${cond}rep${rep}_set_random\$i.${method}.gsa.out; done
  """

  script:
  """
  ${params.load_magma}
  ${params.load_r}
  export SIM_SCRIPTS=${params.sim_scripts}

  Rscript ${params.sim_scripts}/04_batch_random_gmts.R \\
    gmt_dir=${gmt_dir} first=${first} last=${last} \\
    out_file=batch.setannot index_file=batch_index.tsv

  magma --gene-results ${gene_raw} --set-annot batch.setannot --out batch

  mkdir -p split
  Rscript ${params.sim_scripts}/05_split_batch_gsa.R \\
    gsa_out=batch.gsa.out out_dir=split \\
    tag=${cond}rep${rep} method=${method}
  """
}

process calc_empirical {
  tag "${cond}_rep${rep}_${method}"
  publishDir "${params.outdir}/empirical/${cond}/${rep}", mode: 'copy', overwrite: true

  input:
  tuple val(cond), val(rep), val(method), path(real_gsa), path(gmt),
        path(random_files, stageAs: 'random_sets/*')

  output:
  tuple val(cond), val(rep), val(method),
        path("${cond}rep${rep}_${method}_magma_empirical_pvalues.txt"), emit: empirical

  stub:
  """
  touch ${cond}rep${rep}_${method}_magma_empirical_pvalues.txt
  """

  script:
  // The parent pipeline's empirical-p / SES code, called verbatim.
  """
  ${params.load_r}

  n=\$(ls random_sets/*.gsa.out | wc -l)
  echo "Random databases staged: \$n (expected ${params.num_random_sets})"
  [ "\$n" -eq "${params.num_random_sets}" ] || \\
    { echo "ERROR: wrong number of random result files" >&2; exit 1; }

  Rscript ${params.core_scripts}/calc_empirical.r \\
    "${cond}rep${rep}_${method}" "magma" "${real_gsa}" "random_sets" "${gmt}"
  """
}

// Asserts the core feasibility assumption: scoring K null databases in one
// --set-annot file gives the same per-set result as scoring them separately.
process verify_batching {
  tag "${cond}_rep${rep}_${method}"
  publishDir "${params.outdir}/verification", mode: 'copy', overwrite: true

  input:
  tuple val(cond), val(rep), val(method), path(gmt_dir), path(gene_raw)

  output:
  path "batching_equivalence_${cond}_rep${rep}_${method}.txt"

  stub:
  """
  touch batching_equivalence_${cond}_rep${rep}_${method}.txt
  """

  script:
  """
  ${params.load_magma}
  ${params.load_r}
  export SIM_SCRIPTS=${params.sim_scripts}

  N=5
  magma --gene-results ${gene_raw} \\
        --set-annot ${gmt_dir}/GeneSet.random1.gmt --out standalone

  Rscript ${params.sim_scripts}/04_batch_random_gmts.R \\
    gmt_dir=${gmt_dir} first=1 last=\$N out_file=batch.setannot index_file=idx.tsv
  magma --gene-results ${gene_raw} --set-annot batch.setannot --out batched
  mkdir -p split
  Rscript ${params.sim_scripts}/05_split_batch_gsa.R \\
    gsa_out=batched.gsa.out out_dir=split tag=v method=${method}

  Rscript -e '
    suppressPackageStartupMessages(library(data.table))
    rd <- function(f){ x <- readLines(f); h <- which(!startsWith(x,"#") & nzchar(trimws(x)))[1]
      d <- fread(text=paste(x[h:length(x)],collapse="\\n"))
      if (!"FULL_NAME" %in% names(d)) d[, FULL_NAME := VARIABLE]
      setkey(d, FULL_NAME); d[] }
    a <- rd("standalone.gsa.out"); b <- rd("split/v_set_random1.${method}.gsa.out")
    j <- merge(a, b, by="FULL_NAME", suffixes=c(".s",".b"))
    dP <- max(abs(j\$P.s - j\$P.b)); dB <- max(abs(j\$BETA.s - j\$BETA.b))
    dN <- max(abs(j\$NGENES.s - j\$NGENES.b))
    ok <- nrow(j)==nrow(a) && dP < 1e-12 && dB < 1e-12 && dN == 0
    out <- c(sprintf("sets compared        : %d (standalone %d)", nrow(j), nrow(a)),
             sprintf("max |dP|             : %.3g", dP),
             sprintf("max |dBETA|          : %.3g", dB),
             sprintf("max |dNGENES|        : %g", dN),
             sprintf("BATCHING EQUIVALENT  : %s", ok))
    writeLines(out, "batching_equivalence_${cond}_rep${rep}_${method}.txt")
    cat(paste(out, collapse="\\n"), "\\n")
    if (!ok) { cat("\\nERROR: batched MAGMA results differ from standalone.\\n"); quit(status=1) }
  '
  """
}
