// ---------------------------------------------------------------------------
// Analysis stages: per-replicate metrics, sweep aggregation, figures,
// BiRewire diagnostics, and the mu calibration fit.
// ---------------------------------------------------------------------------

process metrics {
  tag "${cond}_rep${rep}"
  publishDir path: { "${params.outdir}/metrics/${cond}/${rep}" }, mode: 'copy', overwrite: true

  input:
  tuple val(cond), val(rep), path(truth), path(real_gsa),
        path(emp_birewire), path(emp_keeppath)

  output:
  path "${cond}_rep${rep}_metrics.tsv",       emit: metrics
  path "${cond}_rep${rep}_pathway_stats.tsv", emit: pathway_stats


  script:
  """
  ${params.load_r}
  export SIM_SCRIPTS=${params.sim_scripts}

  Rscript ${params.sim_scripts}/06_metrics.R \\
    ground_truth=${truth} \\
    real_gsa=${real_gsa} \\
    emp_birewire=${emp_birewire} \\
    emp_keeppath=${emp_keeppath} \\
    condition=${cond} replicate=${rep} \\
    alpha=${params.alpha} top_k=${params.top_k} \\
    out_prefix=${cond}_rep${rep}
  """

  stub:
  """
  touch ${cond}_rep${rep}_metrics.tsv ${cond}_rep${rep}_pathway_stats.tsv
  """
}

// Calibration pilot needs Original only (raw MAGMA p), so it skips
// randomization and the empirical step entirely.
process metrics_original {
  tag "${cond}_rep${rep}"

  input:
  tuple val(cond), val(rep), path(truth), path(real_gsa)

  output:
  path "${cond}_rep${rep}_metrics.tsv", emit: metrics


  script:
  """
  ${params.load_r}
  export SIM_SCRIPTS=${params.sim_scripts}

  Rscript ${params.sim_scripts}/06_metrics.R \\
    ground_truth=${truth} \\
    real_gsa=${real_gsa} \\
    condition=${cond} replicate=${rep} \\
    alpha=${params.alpha} top_k=${params.top_k} \\
    out_prefix=${cond}_rep${rep}
  """

  stub:
  """
  touch ${cond}_rep${rep}_metrics.tsv
  """
}

process aggregate_results {
  publishDir "${params.outdir}/summary", mode: 'copy', overwrite: true

  input:
  path metrics_files, stageAs: 'metrics/*'
  path pathway_files, stageAs: 'pathways/*'

  output:
  path "summary_metrics.csv",       emit: summary
  path "all_replicate_metrics.csv", emit: replicates
  path "all_pathway_stats.csv",     emit: pathways
  path "hub_regression.csv",        emit: hubreg
  path "sanity_checks.txt"


  script:
  """
  ${params.load_r}
  export SIM_SCRIPTS=${params.sim_scripts}

  Rscript ${params.sim_scripts}/07_aggregate.R in_dir=. out_dir=.
  """

  stub:
  """
  touch summary_metrics.csv all_replicate_metrics.csv all_pathway_stats.csv hub_regression.csv sanity_checks.txt
  """
}

process figures {
  publishDir "${params.outdir}/summary", mode: 'copy', overwrite: true

  input:
  tuple path(summary), path(replicates), path(pathways), path(hubreg)

  output:
  path "figures/*"


  script:
  """
  ${params.load_r}
  export SIM_SCRIPTS=${params.sim_scripts}

  Rscript ${params.sim_scripts}/08_figures.R \\
    summary_csv=${summary} replicate_csv=${replicates} \\
    pathway_csv=${pathways} hubreg_csv=${hubreg} \\
    out_dir=figures \\
    overlap_map="${params.overlap_map}" \\
    main_cond=${params.main_cond} nullhub_cond=${params.nullhub_cond}
  """

  stub:
  """
  mkdir -p figures && touch figures/fig1.png
  """
}

process birewire_diagnostics {
  tag "${cond}_rep${rep}"
  publishDir "${params.outdir}/birewire_diagnostics", mode: 'copy', overwrite: true

  input:
  tuple val(cond), val(rep), path(gmt)

  output:
  path "bw_${cond}_*"


  script:
  """
  ${params.load_r}
  export SIM_SCRIPTS=${params.sim_scripts}

  Rscript ${params.sim_scripts}/09_birewire_diagnostics.R \\
    gmt=${gmt} out_prefix=bw_${cond} seed=${rep}
  """

  stub:
  """
  touch bw_${cond}_report.txt bw_${cond}_convergence.csv
  """
}

process fit_calibration {
  publishDir "${params.outdir}/calibration", mode: 'copy', overwrite: true

  input:
  path metrics_files, stageAs: 'metrics/*'
  path grid

  output:
  path "calibration_mu.tsv"
  path "calibration_mu.env", emit: mu_env
  path "calibration_report.txt"
  path "calibration_power_curve.csv"
  path "calibration_clean_null_fpr.csv"


  script:
  """
  ${params.load_r}
  export SIM_SCRIPTS=${params.sim_scripts}

  Rscript ${params.sim_scripts}/10_calibrate.R \\
    pilot_metrics=metrics pilot_grid=${grid} \\
    out_prefix=calibration
  """

  stub:
  """
  touch calibration_mu.tsv calibration_mu.env calibration_report.txt calibration_power_curve.csv calibration_clean_null_fpr.csv
  """
}
