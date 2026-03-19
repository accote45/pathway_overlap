#!/bin/bash
# ============================================================
# Partitioned heritability LDSC — per-trait script
# Munges sumstats inline, then submits partitioned h2 LSF job.
# Assumes setup_sldsc.sh has already been run and LD scores exist.
#
# Usage:
#   bash run_partitioned_h2.sh \
#       --trait    ibd \
#       --gwas     /path/to/gwas.txt \
#       --snp      Rsid \
#       --a1       Allele1 \
#       --a2       Allele2 \
#       --p        P.value \
#       --n-col    N \
#       --beta-col Effect
# ============================================================

set -euo pipefail

# --- Shared configuration ------------------------------------
PROJ=acc_paul_oreilly
QUEUE=premium

WORKDIR=/sc/arion/projects/psychgen/cotea02_prset/pathway_overlap/exploratory/sldsc

SNPLIST=${WORKDIR}/reference/w_hm3.snplist
BASELINE_LDSCORES=${WORKDIR}/reference/baseline/baselineLD.
WLDCHR=${WORKDIR}/reference/weights/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
FRQFILE_CHR=${WORKDIR}/reference/frq/1000G_Phase3_frq/1000G.EUR.QC.

# --- Parse arguments -----------------------------------------
TRAIT=""
GWAS_FILE=""
SNP_COL="SNP"
A1_COL="A1"
A2_COL="A2"
P_COL="P"
N_COL="N"
BETA_COL="BETA"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --trait)    TRAIT="$2";     shift 2 ;;
        --gwas)     GWAS_FILE="$2"; shift 2 ;;
        --snp)      SNP_COL="$2";   shift 2 ;;
        --a1)       A1_COL="$2";    shift 2 ;;
        --a2)       A2_COL="$2";    shift 2 ;;
        --p)        P_COL="$2";     shift 2 ;;
        --n-col)    N_COL="$2";     shift 2 ;;
        --beta-col) BETA_COL="$2";  shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [[ -z "${TRAIT}" || -z "${GWAS_FILE}" ]]; then
    echo "Error: --trait and --gwas are required."
    exit 1
fi

cd "${WORKDIR}"
mkdir -p logs results scripts_temp

submit() { bsub < "$1" | awk '{print $2}' | tr -d '<'; }


# ============================================================
# INLINE: Munge sumstats
# ============================================================
echo "[1/2] Munging sumstats for ${TRAIT}..."
ml ldsc/1.0.1

munge_sumstats.py \
    --sumstats "${GWAS_FILE}" \
    --merge-alleles "${SNPLIST}" \
    --snp "${SNP_COL}" \
    --a1  "${A1_COL}" \
    --a2  "${A2_COL}" \
    --p   "${P_COL}" \
    --N-col "${N_COL}" \
    --signed-sumstats "${BETA_COL}",0 \
    --chunksize 500000 \
    --out "${TRAIT}"

echo "[2/2] Submitting partitioned h2 job for ${TRAIT}..."


# ============================================================
# LSF JOB: Partitioned heritability
# No LSF dependency — assumes LD scores already exist on disk
# ============================================================
cat > scripts_temp/stage_part_h2_${TRAIT}.sh << EOF
#!/bin/bash
#BSUB -J sldsc_part_h2_${TRAIT}
#BSUB -P ${PROJ}
#BSUB -q ${QUEUE}
#BSUB -n 1
#BSUB -R "rusage[mem=16000]"
#BSUB -W 2:00
#BSUB -o ${WORKDIR}/logs/o.sldsc_part_h2_${TRAIT}
#BSUB -e ${WORKDIR}/logs/e.sldsc_part_h2_${TRAIT}

cd ${WORKDIR}
ml ldsc/1.0.1

ldsc.py \\
    --h2 ${TRAIT}.sumstats.gz \\
    --ref-ld-chr ${BASELINE_LDSCORES},annot/multipathway.,annot/othergenes. \\
    --frqfile-chr ${FRQFILE_CHR} \\
    --w-ld-chr ${WLDCHR} \\
    --overlap-annot \\
    --print-coefficients \\
    --print-delete-vals \\
    --out results/${TRAIT}.multipathway_partitioned_h2
EOF
JOB_H2=$(submit scripts_temp/stage_part_h2_${TRAIT}.sh)

echo ""
echo "============================================================"
echo "Submitted sldsc_part_h2_${TRAIT}: job ID ${JOB_H2}"
echo "Monitor with: bjobs -w"
echo "Output: results/${TRAIT}.multipathway_partitioned_h2"
echo "============================================================"

# IMPORTANT — Symlink requirement for LDSC prefixes:
# LDSC uses glob.glob() internally to validate --ref-ld-chr and --w-ld-chr prefixes.
# Bare filename prefixes (e.g. "baselineLD.") don't match as filesystem entries,
# causing: ValueError: No objects to concatenate
# Fix: create a self-referential symlink in each LD score directory:
#   cd ldscores/baseline/ && ln -s . "baselineLD."
#   cd annot/            && ln -s . "multipathway." && ln -s . "othergenes."
#   cd weights_dir/      && ln -s . "weights.hm3_noMHC."
# This only needs to be done once per directory