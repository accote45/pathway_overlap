#!/bin/bash
# ============================================================
# Partitioned heritability LDSC — one-time reference setup
# Builds gene lists, allele freqs, SNP annotations, and LD scores.
# Run this once before any per-trait analyses.
#
# Dependency graph:
#   [inline] extract_refs, hm3_snpids, gene_lists, gene_coords, annot[1-22]
#            → LSF: sldsc_l2[1-22]
# ============================================================

set -euo pipefail

# --- Shared configuration ------------------------------------
PROJ=acc_paul_oreilly
QUEUE=premium

WORKDIR=/sc/arion/projects/psychgen/cotea02_prset/pathway_overlap/exploratory/sldsc
LDSCORE_DIR=/sc/arion/projects/data-ark/Public_Unrestricted/LDSCORE
GMT_FILE=/sc/arion/projects/psychgen/cotea02_prset/pathway_overlap/data/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt

# Reference paths (populated by extraction step below)
BIMDIR=${WORKDIR}/reference/plink/1000G_EUR_Phase3_plink
SNPLIST=${WORKDIR}/reference/w_hm3.snplist

cd "${WORKDIR}"
mkdir -p logs annot data/gene_lists data/gene_coords results scripts_temp reference/plink reference/baseline reference/weights reference/frq

submit() { bsub < "$1" | awk '{print $2}' | tr -d '<'; }


# ============================================================
# INLINE: Extract LDSCORE references (skips if already done)
# ============================================================
echo "[0/4] Extracting LDSCORE references..."

# HapMap3 SNP list
if [[ ! -f "${SNPLIST}" ]]; then
    bunzip2 -c "${LDSCORE_DIR}/w_hm3.snplist.bz2" > "${WORKDIR}/reference/w_hm3.snplist"
fi

# 1000G Phase3 PLINK bfiles (needed for make_annot.py and ldsc.py --l2)
# NOTE: verify subdir name after extraction: ls reference/plink/
if [[ ! -d "${BIMDIR}" ]]; then
    tar -xzf "${LDSCORE_DIR}/1000G_Phase3_plinkfiles.tgz" -C "${WORKDIR}/reference/plink/"
fi

# Pre-computed allele frequencies (replaces inline PLINK step)
# NOTE: verify prefix pattern after extraction matches 1000G.EUR.QC.{chr}.frq
if [[ ! -d "${WORKDIR}/reference/frq/1000G_Phase3_frq" ]] || [[ -z "$(ls -A ${WORKDIR}/reference/frq/1000G_Phase3_frq)" ]]; then
    tar -xzf "${LDSCORE_DIR}/1000G_Phase3_frq.tgz" -C "${WORKDIR}/reference/frq/"
fi

# Baseline-LD v2.2 scores (recommended for h2 enrichment per README)
if [[ ! -d "${WORKDIR}/reference/baseline" ]] || [[ -z "$(ls -A ${WORKDIR}/reference/baseline)" ]]; then
    tar -xzf "${LDSCORE_DIR}/1000G_Phase3_baselineLD_v2.2_ldscores.tgz" -C "${WORKDIR}/reference/baseline/"
fi

# Regression weights
if [[ ! -d "${WORKDIR}/reference/weights" ]] || [[ -z "$(ls -A ${WORKDIR}/reference/weights)" ]]; then
    tar -xzf "${LDSCORE_DIR}/1000G_Phase3_weights_hm3_no_MHC.tgz" -C "${WORKDIR}/reference/weights/"
fi

# Create self-referential symlinks required by LDSC's glob.glob() validation
# (must match the prefix string passed to --ref-ld-chr and --w-ld-chr)
[[ ! -L "${WORKDIR}/reference/baseline/baselineLD." ]] && \
    ln -s . "${WORKDIR}/reference/baseline/baselineLD."
[[ ! -L "${WORKDIR}/reference/weights/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC." ]] && \
    ln -s . "${WORKDIR}/reference/weights/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
[[ ! -L "${WORKDIR}/annot/multipathway." ]] && \
    ln -s . "${WORKDIR}/annot/multipathway."
[[ ! -L "${WORKDIR}/annot/othergenes." ]] && \
    ln -s . "${WORKDIR}/annot/othergenes."


# ============================================================
# INLINE: Extract HM3 SNP IDs (needed by LD score jobs)
# ============================================================
echo "[1/4] Extracting HM3 SNP IDs..."
awk 'NR==1{for(i=1;i<=NF;i++) if($i=="SNP") col=i; print "SNP"; next} {print $col}' \
    "${SNPLIST}" > data/hm3_snpids_only.txt


# ============================================================
# INLINE: Gene lists + gene coords
# ============================================================
echo "[2/4] Building gene lists and coords..."
ml R

Rscript - <<EOF
library(data.table)
gmt_file <- "${GMT_FILE}"
lines <- readLines(gmt_file)
genes <- unlist(lapply(lines, function(line) {
  fields <- strsplit(line, "\t")[[1]]
  fields[-(1:2)]
}))
counts <- sort(table(genes))
threshold <- quantile(counts, 0.90)
multipathway_genes <- names(counts[counts >= threshold])
cat(sprintf("Threshold: %d pathways\n", threshold))
cat(sprintf("%d multipathway genes\n", length(multipathway_genes)))
writeLines(multipathway_genes, "data/gene_lists/multipathway_genes.txt")
other_genes <- names(counts[counts < threshold])
writeLines(other_genes, "data/gene_lists/other_genes.txt")
cat(sprintf("%d other genes\n", length(other_genes)))
EOF

Rscript - <<'EOF'
library(data.table)
coords <- fread(
  "data/msigdbgenes.regions",
  col.names = c("ensembl_id", "CHR", "START", "END", "strand", "GENE"),
  colClasses = c("character", "character", "integer", "integer", "character", "character")
)
coords <- coords[CHR %in% c(as.character(1:22))]
fwrite(
  coords[, .(GENE = ensembl_id, CHR, START, END)],
  file = "data/gene_coords/gene_coords.txt",
  sep = "\t"
)
EOF


# ============================================================
# INLINE: Generate SNP annotations [chr 1-22]
# (allele freq step removed — using pre-computed frq from LDSCORE)
# ============================================================
echo "[3/4] Generating SNP annotations (chr 1-22)..."
ml ldsc/1.0.1 bedtools

for chr in {1..22}; do
    make_annot.py \
        --gene-set-file data/gene_lists/multipathway_genes.txt \
        --gene-coord-file data/gene_coords/gene_coords.txt \
        --windowsize 35000 \
        --bimfile ${BIMDIR}/1000G.EUR.QC.${chr}.bim \
        --annot-file annot/multipathway.${chr}.annot.gz

    make_annot.py \
        --gene-set-file data/gene_lists/other_genes.txt \
        --gene-coord-file data/gene_coords/gene_coords.txt \
        --windowsize 35000 \
        --bimfile ${BIMDIR}/1000G.EUR.QC.${chr}.bim \
        --annot-file annot/othergenes.${chr}.annot.gz
done

echo "All inline steps complete. Submitting LD score job..."


# ============================================================
# LSF JOB: LD score computation array [1-22]
# ============================================================
cat > scripts_temp/stage_l2.sh << EOF
#!/bin/bash
#BSUB -J sldsc_l2[1-22]
#BSUB -P ${PROJ}
#BSUB -q ${QUEUE}
#BSUB -n 1
#BSUB -R "rusage[mem=16000]"
#BSUB -W 6:00
#BSUB -o ${WORKDIR}/logs/o.sldsc_l2.%I
#BSUB -e ${WORKDIR}/logs/e.sldsc_l2.%I

cd ${WORKDIR}
ml ldsc/1.0.1
chr=\${LSB_JOBINDEX}

ldsc.py \\
    --l2 \\
    --bfile ${BIMDIR}/1000G.EUR.QC.\${chr} \\
    --ld-wind-cm 1 \\
    --thin-annot \\
    --annot annot/multipathway.\${chr}.annot.gz \\
    --print-snp data/hm3_snpids_only.txt \\
    --out annot/multipathway.\${chr}

ldsc.py \\
    --l2 \\
    --bfile ${BIMDIR}/1000G.EUR.QC.\${chr} \\
    --ld-wind-cm 1 \\
    --thin-annot \\
    --annot annot/othergenes.\${chr}.annot.gz \\
    --print-snp data/hm3_snpids_only.txt \\
    --out annot/othergenes.\${chr}
EOF
JOB_L2=$(submit scripts_temp/stage_l2.sh)

echo ""
echo "============================================================"
echo "Submitted sldsc_l2[1-22]: job ID ${JOB_L2}"
echo "Monitor with: bjobs -w"
echo "Once complete, run run_partitioned_h2.sh for each trait."
echo "============================================================"