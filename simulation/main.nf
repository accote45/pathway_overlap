#!/usr/bin/env nextflow
// ============================================================================
//  GSR simulation sub-pipeline (see simulation_instructions.md)
//
//  Ground-truth simulation showing that pathway-overlap-driven false positives
//  are removed by Gene Swap Randomization while true enrichment is retained.
//
//  Self-contained: reads the parent repository's scripts/core/*.R by path but
//  changes nothing outside this directory.
//
//  Run order:
//    1) nextflow run main.nf -entry calibrate      # solve for mu (Sec. 3)
//    2) paste results/calibration/calibration_mu.env into nextflow.config
//    3) nextflow run main.nf -resume               # the sweep
// ============================================================================

nextflow.enable.dsl = 2

include { prepare_reference; magma_annotate;
          build_architecture; simulate_sumstats } from './modules/sim_setup.nf'
include { magma_gene_analysis; randomize_gmts; magma_geneset_real;
          magma_geneset_random_batch; calc_empirical;
          verify_batching } from './modules/sim_enrichment.nf'
include { metrics; metrics_original; aggregate_results; figures;
          birewire_diagnostics; fit_calibration } from './modules/sim_analysis.nf'

// ---------------------------------------------------------------------------
// Shared front half: reference -> architecture -> sumstats -> gene analysis
// ---------------------------------------------------------------------------
workflow SIMULATE {
  take:
  cells          // tuple(cond, rep, overlap_frac, hub_signal, mu_low, mu_med, mu_high)

  main:
  prepare_reference()
  magma_annotate(prepare_reference.out.snp_loc)

  build_architecture(cells.combine(prepare_reference.out.gene_pool))

  simulate_sumstats(
    build_architecture.out.signal.combine(prepare_reference.out.reference))

  magma_gene_analysis(
    simulate_sumstats.out.sumstats.combine(magma_annotate.out.annot))

  magma_geneset_real(
    magma_gene_analysis.out.raw.join(build_architecture.out.gmt, by: [0, 1]))

  emit:
  gmt      = build_architecture.out.gmt
  truth    = build_architecture.out.truth
  gene_raw = magma_gene_analysis.out.raw
  real_gsa = magma_geneset_real.out.real_gsa
}

// ---------------------------------------------------------------------------
// Main sweep
// ---------------------------------------------------------------------------
workflow {
  if (params.mu_low == null || params.mu_med == null || params.mu_high == null) {
    error """
    mu_low / mu_med / mu_high are unset.

    Instructions Sec. 3 requires these to be CALIBRATED, not guessed. Run:
        nextflow run main.nf -entry calibrate
    then copy results/calibration/calibration_mu.env into nextflow.config.
    """.stripIndent()
  }

  log.info "GSR simulation sweep"
  log.info "  conditions      : ${params.conditions.collect{ it.name }.join(', ')}"
  log.info "  replicates      : ${params.n_replicates}"
  log.info "  null databases  : ${params.num_random_sets} per method"
  log.info "  mu low/med/high : ${params.mu_low} / ${params.mu_med} / ${params.mu_high}"
  log.info "  effect scale s  : ${params.effect_scale}"
  if (params.independent_z) {
    log.warn "independent_z = true: Z is drawn WITHOUT LD. Smoke-test only - " +
             "these results must not appear in final figures (Sec. 4)."
  }

  cells = Channel.fromList(params.conditions)
    .combine(Channel.of(1..params.n_replicates))
    .map { c, r -> tuple(c.name, r, c.overlap_frac, c.hub_signal,
                         params.mu_low, params.mu_med, params.mu_high) }

  SIMULATE(cells)

  // ---- null databases, per condition x replicate x method ----------------
  // The randomized databases depend only on the architecture, so they are
  // generated once per (condition, replicate) and reused across all batches.
  randomize_gmts(
    SIMULATE.out.gmt.combine(Channel.fromList(params.randomization_methods)))

  n_batches = (int) Math.ceil(params.num_random_sets / (double) params.batch_size)
  batches = Channel.of(0..<n_batches).map { b ->
    tuple(b + 1,
          b * params.batch_size + 1,
          Math.min((b + 1) * params.batch_size, params.num_random_sets))
  }

  magma_geneset_random_batch(
    randomize_gmts.out.gmt_dir
      .combine(SIMULATE.out.gene_raw, by: [0, 1])
      .combine(batches))

  random_grouped = magma_geneset_random_batch.out.split
    .groupTuple(by: [0, 1, 2], size: n_batches)
    .map { c, r, m, fl -> tuple(c, r, m, fl.flatten()) }

  calc_empirical(
    random_grouped
      .combine(SIMULATE.out.real_gsa, by: [0, 1])
      .combine(SIMULATE.out.gmt,      by: [0, 1])
      .map { c, r, m, files, gsa, gmt -> tuple(c, r, m, gsa, gmt, files) })

  // ---- metrics -----------------------------------------------------------
  emp = calc_empirical.out.empirical
  emp_bw = emp.filter { it[2] == 'birewire'     }.map { c, r, m, f -> tuple(c, r, f) }
  emp_kp = emp.filter { it[2] == 'keeppathsize' }.map { c, r, m, f -> tuple(c, r, f) }

  metrics(
    SIMULATE.out.truth
      .join(SIMULATE.out.real_gsa, by: [0, 1])
      .join(emp_bw, by: [0, 1])
      .join(emp_kp, by: [0, 1]))

  aggregate_results(metrics.out.metrics.collect(),
                    metrics.out.pathway_stats.collect())

  figures(aggregate_results.out.summary
            .combine(aggregate_results.out.replicates)
            .combine(aggregate_results.out.pathways)
            .combine(aggregate_results.out.hubreg))

  // ---- controls and diagnostics ------------------------------------------
  if (params.run_birewire_diagnostics) {
    // one representative architecture per condition
    birewire_diagnostics(SIMULATE.out.gmt.filter { it[1] == 1 })
  }

  if (params.run_batch_verification) {
    verify_batching(
      randomize_gmts.out.gmt_dir
        .combine(SIMULATE.out.gene_raw, by: [0, 1])
        .filter { it[1] == 1 && it[2] == 'birewire' }
        .first())
  }
}

// ---------------------------------------------------------------------------
// Calibration pilot (instructions Sec. 3)
//
// 0% overlap, all three tiers sharing one mu, swept over the grid. Calibrates
// on Original (raw MAGMA p), so no null databases are needed.
// ---------------------------------------------------------------------------
workflow calibrate {
  log.info "GSR simulation - mu calibration pilot"
  log.info "  grid       : ${params.calibration_mu_grid.join(', ')}"
  log.info "  replicates : ${params.n_replicates}"

  grid = Channel.fromList(params.calibration_mu_grid.withIndex()
           .collect { mu, i -> tuple("cal${i + 1}", mu) })

  cells = grid
    .combine(Channel.of(1..params.n_replicates))
    .map { name, mu, r -> tuple(name, r, 0.0, true, mu, mu, mu) }

  SIMULATE(cells)

  metrics_original(SIMULATE.out.truth.join(SIMULATE.out.real_gsa, by: [0, 1]))

  grid_file = grid
    .map { name, mu -> "${name}\t${mu}" }
    .collectFile(name: 'pilot_grid.tsv', newLine: true, sort: true,
                 seed: "condition\tmu")

  fit_calibration(metrics_original.out.metrics.collect(), grid_file)
}
