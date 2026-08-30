#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# GSR simulation - end-to-end runner (instructions Sec. 10).
#
# Run from this directory on the cluster. Each stage is resumable; re-running is
# safe. Stages are deliberately separate because stage 2 requires a human to
# paste the calibrated mu values into nextflow.config.
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "$0")"

STAGE="${1:-help}"
NF="${NF:-nextflow}"

case "$STAGE" in

  smoke)
    # Sec. 4 fast fallback: iid Z, tiny sweep. Proves the plumbing only -
    # these results must NOT appear in any figure.
    echo ">>> SMOKE TEST (independent Z, no LD - plumbing check only)"
    $NF run main.nf -resume \
      --n_replicates 2 --num_random_sets 20 --batch_size 10 \
      --independent_z true \
      --mu_low 0.10 --mu_med 0.20 --mu_high 0.35 \
      --outdir results_smoke
    ;;

  calibrate)
    # Sec. 3: solve for mu at 0% overlap against Original's power curve.
    echo ">>> CALIBRATION PILOT"
    $NF run main.nf --stage calibrate -resume --n_replicates 10
    echo
    echo "Calibrated values:"
    cat results/calibration/calibration_mu.env
    echo
    echo "NEXT: paste those three lines into nextflow.config as"
    echo "      params.mu_low / params.mu_med / params.mu_high,"
    echo "      then run:  ./run_all.sh sweep"
    ;;

  sweep)
    # Sec. 10 step 1 (the MVP): overlap {0%, 5%} + null-hub control.
    echo ">>> MAIN SWEEP"
    $NF run main.nf -resume
    echo
    echo "Figures + tidy CSVs: results/summary/figures/"
    echo "Sanity checks      : results/summary/sanity_checks.txt"
    echo "Batching check     : results/verification/"
    ;;

  extended)
    # Sec. 10 step 3: add the 10% / 20% overlap levels. Figure 1 needs at least
    # TWO non-zero overlap levels to be a line rather than a single point.
    echo ">>> EXTENDED SWEEP (adds 10% and 20% overlap)"
    $NF run main.nf -resume \
      --overlap_map "ov00=0,ov05=5,ov10=10,ov20=20,ov05_nullhub=5" \
      -c conf/extended.config
    ;;

  *)
    cat <<'EOF'
Usage: ./run_all.sh <stage>

  smoke      tiny end-to-end run with independent Z (plumbing check only)
  calibrate  solve for mu_low / mu_med / mu_high  (RUN THIS FIRST)
  sweep      the MVP sweep: overlap {0%, 5%} + null-hub control
  extended   adds the 10% and 20% overlap conditions

Typical order:  ./run_all.sh calibrate  ->  edit nextflow.config  ->  ./run_all.sh sweep
EOF
    exit 1
    ;;
esac
