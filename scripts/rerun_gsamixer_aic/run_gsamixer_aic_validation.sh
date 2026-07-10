#!/usr/bin/env bash
# ============================================================================
# run_gsamixer_aic_validation.sh
#
# Re-run ONLY the GSA-MiXeR validation analyses on an OLD results directory,
# switching the top-500 pathway pool from "union of top pathways across
# PRSet/MAGMA/PascalX" to "top 500 by GSA-MiXeR's native AIC".
#
# Nothing upstream is recomputed: empirical p-values and standardized-effect-
# size / GSR corrections are taken verbatim from the old run. We only:
#   1. append a `mixer_aic` column to a COPY of each empirical table
#      (joined from the raw real *_full.go_test_enrich.csv), then
#   2. re-run the 4 validation scripts (which already select top-500 by AIC).
#
# Run this ON THE CLUSTER, interactively (R on PATH). No scheduler needed.
# ============================================================================
set -euo pipefail

# ----------------------------------------------------------------------------
# CONFIG  --  edit these to match the cluster
# ----------------------------------------------------------------------------
OLD_DIR="/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf"

# This repo's scripts dir -- MUST be the current (AIC-aware) validation scripts,
# NOT the old project's union-based copies. geneoverlap_nf is NOT a git clone, so
# its scripts/ holds the OLD union-based versions -- do NOT point here. Instead
# scp this repo's current scripts/ to the location below (see walkthrough).
SCRIPTS_DIR="${OLD_DIR}/aic_rerun_code/scripts"   # <-- current (finalize-for-release) scripts, copied here

# Where patched tables + new validation outputs go (non-destructive).
OUT_DIR="${OLD_DIR}/aic_topselection_rerun"

# Reference data -- recovered from this repo's git history (commit ce9c909^,
# before paths were replaced with placeholders for release).
GENESET_GMT="/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt"
TISSUE_DATA="/sc/arion/projects/psychgen/cotea02_prset/judit_revisions/software/1kg_test/GeneExpressionLandscape/data/Exp_Spe_DataTables/specificity"
OT_JSON_DIR="/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/results/drugtarget_test/associationByDatatypeDirect"
MALACARDS_PATH="/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/Malacards"
DOROTHEA_PATH="/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/dorothea_pairwise_scores.csv"

RAND_METHODS=(birewire keeppathsize)

# Traits with GSA-MiXeR results to re-validate (all traits get tissue + dorothea).
# The 4 quantitative traits are only used by the tissue + dorothea panels (the OT
# and MalaCards aggregators don't include them), but they must be here or gsamixer
# is blank for them in those two panels. Traits with no gsamixer inputs skip safely.
TRAITS=(ibd scz bmi cad t2d mdd breast ad \
        Alkaline_phosphatase Eosinophill_percentage HDL_cholesterol Mean_platelet_thrombocyte_volume)

# Per-dataset trait filters (mirror the Nextflow workflow).
OT_TRAITS=(t2d cad ad mdd scz ibd breast bmi)
MALACARDS_TRAITS=(bmi cad t2d mdd ad scz ibd breast)

# Which validation datasets to run (set to 0 to skip one).
RUN_OPENTARGETS=1
RUN_TISSUE=1
RUN_MALACARDS=1
RUN_DOROTHEA=1

# ----------------------------------------------------------------------------
# Path patterns in the OLD dir  --  edit if your layout differs
# ----------------------------------------------------------------------------
# Raw real result (confirmed from your example): results/gsamixer/<trait>/<trait>_full.go_test_enrich.csv
real_file()  { echo "${OLD_DIR}/results/gsamixer/$1/$1_full.go_test_enrich.csv"; }

# Old GSA-MiXeR standardized-effects table for a trait + randomization method.
# Confirmed filename: <trait>_gsamixer_<rm>_standardized_effects.txt
# The directory below is a best guess; if it misses, resolve_emp() falls back to
# a `find` under OLD_DIR, so a wrong directory won't stop the run.
EMP_FILENAME() { echo "$1_gsamixer_$2_standardized_effects.txt"; }
emp_file()     { echo "${OLD_DIR}/results/empirical_pvalues/gsamixer_$2/$1/$(EMP_FILENAME "$1" "$2")"; }

# Resolve the empirical table: use the pattern if present, else search OLD_DIR.
resolve_emp() {
  local direct; direct="$(emp_file "$1" "$2")"
  if [[ -f "$direct" ]]; then echo "$direct"; return 0; fi
  find "${OLD_DIR}" -type f -name "$(EMP_FILENAME "$1" "$2")" 2>/dev/null | head -1
}

# ----------------------------------------------------------------------------
PATCH="${SCRIPTS_DIR}/rerun_gsamixer_aic/patch_add_mixer_aic.R"
V_OT="${SCRIPTS_DIR}/validation/opentargets/OT_correlation_stats.R"
V_TIS="${SCRIPTS_DIR}/validation/tissue/tissue_correlation_stats.R"
V_MAL="${SCRIPTS_DIR}/validation/malacards/malacards_correlation.R"
V_DOR="${SCRIPTS_DIR}/validation/dorothea/dorothea_correlation.R"

in_list() { local x="$1"; shift; for e in "$@"; do [[ "$e" == "$x" ]] && return 0; done; return 1; }

mkdir -p "${OUT_DIR}" "${OUT_DIR}/patched_empirical"

# ----------------------------------------------------------------------------
# PREFLIGHT: fail fast if the first trait's inputs are not where we expect.
# ----------------------------------------------------------------------------
echo "== Preflight =="
t0="${TRAITS[0]}"
pf_real="$(real_file "$t0")"
pf_bw="$(resolve_emp "$t0" birewire)"
pf_kp="$(resolve_emp "$t0" keeppathsize)"
for f in "$pf_real" "$pf_bw" "$pf_kp"; do
  if [[ -n "$f" && -f "$f" ]]; then echo "  OK   $f"; else echo "  MISSING  (trait ${t0}) -- fix path patterns at the top."; exit 1; fi
done
for s in "$PATCH" "$V_OT" "$V_TIS" "$V_MAL" "$V_DOR"; do
  [[ -f "$s" ]] || { echo "  MISSING script: $s"; exit 1; }
done
echo "  Preflight passed."
echo

# ----------------------------------------------------------------------------
# MAIN LOOP
# ----------------------------------------------------------------------------
for trait in "${TRAITS[@]}"; do
  echo "==================== ${trait} ===================="
  real="$(real_file "$trait")"
  if [[ ! -f "$real" ]]; then echo "  skip: no real file ($real)"; continue; fi

  # --- Step 1: patch each rand-method empirical table with mixer_aic ---
  declare -A patched=()
  ok=1
  for rm in "${RAND_METHODS[@]}"; do
    emp_in="$(resolve_emp "$trait" "$rm")"
    emp_out="${OUT_DIR}/patched_empirical/${trait}_gsamixer_${rm}_standardized_effects_with_aic.txt"
    if [[ -z "$emp_in" || ! -f "$emp_in" ]]; then echo "  skip ${rm}: no standardized_effects table found"; ok=0; break; fi
    Rscript "$PATCH" "$real" "$emp_in" "$emp_out"
    patched[$rm]="$emp_out"
  done
  [[ $ok -eq 1 ]] || { echo "  skip ${trait}: incomplete empirical tables"; continue; }

  bw="${patched[birewire]}"
  kp="${patched[keeppathsize]}"

  # --- Step 2: run validation scripts in a per-trait output dir ---
  tdir="${OUT_DIR}/${trait}"
  mkdir -p "$tdir"
  pushd "$tdir" >/dev/null

  if [[ "$RUN_TISSUE" == 1 ]]; then
    echo "  -> tissue"
    Rscript "$V_TIS" "$trait" gsamixer "$bw" "$kp" "$GENESET_GMT" "$TISSUE_DATA"
  fi
  if [[ "$RUN_DOROTHEA" == 1 ]]; then
    echo "  -> dorothea"
    Rscript "$V_DOR" "$trait" gsamixer "$DOROTHEA_PATH" "$bw" "$kp"
  fi
  if [[ "$RUN_OPENTARGETS" == 1 ]] && in_list "$trait" "${OT_TRAITS[@]}"; then
    echo "  -> opentargets"
    Rscript "$V_OT" "$trait" gsamixer "$bw" "$kp" "$GENESET_GMT" "$OT_JSON_DIR"
  fi
  if [[ "$RUN_MALACARDS" == 1 ]] && in_list "$trait" "${MALACARDS_TRAITS[@]}"; then
    echo "  -> malacards"
    Rscript "$V_MAL" "$trait" gsamixer "$MALACARDS_PATH" "$bw" "$kp" "$GENESET_GMT"
  fi

  popd >/dev/null
done

echo
echo "Done. Patched tables:   ${OUT_DIR}/patched_empirical/"
echo "      Validation output: ${OUT_DIR}/<trait>/"
