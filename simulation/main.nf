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
//    1) nextflow run main.nf --stage calibrate     # solve for mu (Sec. 3)
//    2) paste results/calibration/calibration_mu.env into nextflow.config
//    3) nextflow run main.nf -resume               # the sweep
//
//  Stage selection is a PARAM, not `-entry`: Nextflow's strict parser (25.x+)
//  requires the entry workflow to be the anonymous one.
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
workflow SWEEP {
  if (params.mu_low == null || params.mu_med == null || params.mu_high == null) {
    error """
    mu_low / mu_med / mu_high are unset.

    Instructions Sec. 3 requires these to be CALIBRATED, not guessed. Run:
        nextflow run main.nf --stage calibrate
    then copy results/calibration/calibration_mu.env into nextflow.config.
    """.stripIndent()
  }

  log.info "GSR simulation sweep"
  log.info "  conditions      : ${params.conditions.collect{ c -> c.name }.join(', ')}"
  log.info "  replicates      : ${params.n_replicates}"
  log.info "  null databases  : ${params.num_random_sets} per method"
  log.info "  mu low/med/high : ${params.mu_low} / ${params.mu_med} / ${params.mu_high}"
  log.info "  effect scale s  : ${params.effect_scale}"
  if (params.independent_z) {
    log.warn "independent_z = true: Z is drawn WITHOUT LD. Smoke-test only - " +
             "these results must not appear in final figures (Sec. 4)."
  }

  def n_reps   = params.n_replicates as Integer
  def n_null   = params.num_random_sets as Integer
  def batch_sz = params.batch_size as Integer

  cells = channel.fromList(params.conditions)
    .combine(channel.of(1..n_reps))
    .map { c, r -> tuple(c.name, r, c.overlap_frac, c.hub_signal,
                         params.mu_low, params.mu_med, params.mu_high) }

  SIMULATE(cells)

  // ---- null databases, per condition x replicate x method ----------------
  // The randomized databases depend only on the architecture, so they are
  // generated once per (condition, replicate) and reused across all batches.
  randomize_gmts(
    SIMULATE.out.gmt.combine(channel.fromList(params.randomization_methods)))

  // Batch boundaries come straight from `collate` rather than any arithmetic:
  // Nextflow 26's runtime does not expose intdiv/Math to the DSL.
  batches = channel.of(1..n_null)
    .collate(batch_sz)
    .map { grp -> tuple(grp[0], grp[-1]) }

  magma_geneset_random_batch(
    randomize_gmts.out.gmt_dir
      .combine(SIMULATE.out.gene_raw, by: [0, 1])
      .combine(batches))

  // No `size:` - the number of batches is now implicit in the collate above,
  // and groupTuple closes each group when the upstream channel completes.
  random_grouped = magma_geneset_random_batch.out.split
    .groupTuple(by: [0, 1, 2])
    .map { c, r, m, fl -> tuple(c, r, m, fl.flatten()) }

  calc_empirical(
    random_grouped
      .combine(SIMULATE.out.real_gsa, by: [0, 1])
      .combine(SIMULATE.out.gmt,      by: [0, 1])
      .map { c, r, m, files, gsa, gmt -> tuple(c, r, m, gsa, gmt, files) })

  // ---- metrics -----------------------------------------------------------
  emp = calc_empirical.out.empirical
  emp_bw = emp.filter { t -> t[2] == 'birewire'     }.map { c, r, _m, f -> tuple(c, r, f) }
  emp_kp = emp.filter { t -> t[2] == 'keeppathsize' }.map { c, r, _m, f -> tuple(c, r, f) }

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
    birewire_diagnostics(SIMULATE.out.gmt.filter { t -> t[1] == 1 })
  }

  if (params.run_batch_verification) {
    verify_batching(
      randomize_gmts.out.gmt_dir
        .combine(SIMULATE.out.gene_raw, by: [0, 1])
        .filter { t -> t[1] == 1 && t[2] == 'birewire' }
        .first())
  }
}

// ---------------------------------------------------------------------------
// Calibration pilot (instructions Sec. 3)
//
// 0% overlap, all three tiers sharing one mu, swept over the grid. Calibrates
// on Original (raw MAGMA p), so no null databases are needed.
// ---------------------------------------------------------------------------
workflow CALIBRATE {
  log.info "GSR simulation - mu calibration pilot"
  log.info "  grid       : ${params.calibration_mu_grid.join(', ')}"
  log.info "  replicates : ${params.n_replicates}"

  grid = channel.fromList(params.calibration_mu_grid.withIndex()
           .collect { mu, i -> tuple("cal${i + 1}", mu) })

  def n_reps = params.n_replicates as Integer

  cells = grid
    .combine(channel.of(1..n_reps))
    .map { name, mu, r -> tuple(name, r, 0.0, true, mu, mu, mu) }

  SIMULATE(cells)

  metrics_original(SIMULATE.out.truth.join(SIMULATE.out.real_gsa, by: [0, 1]))

  grid_file = grid
    .map { name, mu -> "${name}\t${mu}" }
    .collectFile(name: 'pilot_grid.tsv', newLine: true, sort: true,
                 seed: "condition\tmu")

  fit_calibration(metrics_original.out.metrics.collect(), grid_file)
}

// ---------------------------------------------------------------------------
// Entry point. The strict parser requires exactly one anonymous workflow, so
// the stage is selected by `--stage` rather than `-entry`.
// ---------------------------------------------------------------------------
workflow {
  if (params.stage == 'calibrate') {
    CALIBRATE()
  }
  else if (params.stage == 'sweep') {
    SWEEP()
  }
  else {
    error "Unknown --stage '${params.stage}' (expected 'sweep' or 'calibrate')"
  }
}
