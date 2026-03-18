#!/bin/bash
# ============================================================
# Partitioned heritability LDSC — LSF submission script
# Fast steps run inline; only slow steps submitted to LSF.
#
# Dependency graph:
#   [inline] munge, gene_lists, gene_coords, freq[1-22], annot[1-22]
#            → sldsc_l2[1-22] → sldsc_part_h2
# ============================================================

set -euo pipefail

# --- Shared configuration ------------------------------------
PROJ=acc_paul_oreilly
QUEUE=premium

WORKDIR=/sc/arion/projects/psychgen/cotea02_prset/pathway_overlap/exploratory/sldsc
BIMDIR=/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/gsa_mixer/1000G_EUR_Phase3_plink
BASELINE_LDSCORES=/sc/arion/projects/psychgen/cotea02_prset/pathway_overlap/exploratory/sldsc/ldscores/baseline/baselineLD_v2.2.
SNPLIST=/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/results/ldsc/resource/eur_w_ld_chr/w_hm3.snplist
GWAS_FILE=/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/gwas/ibd_gwas.txt
WLDCHR=/sc/arion/projects/psychgen/projects/prs/sample_overlap/analysis/prepare/data/reference/weight-
GMT_FILE=/sc/arion/projects/psychgen/cotea02_prset/pathway_overlap/data/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt

cd "${WORKDIR}"
mkdir -p logs annot data/gene_lists data/gene_coords data/freq results scripts_temp

# Helper: submit a script file and return the numeric job ID
submit() { bsub < "$1" | awk '{print $2}' | tr -d '<>'; }


# ============================================================
# INLINE: Extract HM3 SNP IDs (needed by LD score jobs)
# ============================================================
echo "[1/5] Extracting HM3 SNP IDs..."
awk 'NR==1{for(i=1;i<=NF;i++) if($i=="SNP") col=i; print "SNP"; next} {print $col}' \
    "${SNPLIST}" > data/hm3_snpids_only.txt


# ============================================================
# INLINE: Munge sumstats
# ============================================================
echo "[2/5] Munging sumstats..."
ml ldsc/1.0.1

munge_sumstats.py \
    --sumstats ${GWAS_FILE} \
    --merge-alleles ${SNPLIST} \
    --a1 Allele1 \
    --a2 Allele2 \
    --N-col N \
    --p P.value \
    --snp Rsid \
    --signed-sumstats Effect,0 \
    --chunksize 500000 \
    --out ibd


# ============================================================
# INLINE: Gene lists + gene coords (independent, but fast so sequential)
# ============================================================
echo "[3/5] Building gene lists and coords..."
ml R

Rscript - <<'EOF'
library(data.table)
gmt_file <- "/sc/arion/projects/psychgen/cotea02_prset/pathway_overlap/data/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt"
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
# INLINE: Allele frequency calculation [chr 1-22]
# ============================================================
echo "[4/5] Calculating allele frequencies (chr 1-22)..."
ml plink

for chr in {1..22}; do
    bim=${BIMDIR}/1000G.EUR.QC.${chr}.bim

    awk '$2 in a {print $2; next} !($2 in a) {a[$2]=1; next}' \
        ${bim} > data/freq/dup_snps_chr${chr}.txt

    plink \
        --bfile ${BIMDIR}/1000G.EUR.QC.${chr} \
        --mac 5 \
        --geno 0.01 \
        --exclude data/freq/dup_snps_chr${chr}.txt \
        --make-bed \
        --out data/freq/1000G.EUR.QC.${chr}.filtered \
        --silent

    plink \
        --bfile data/freq/1000G.EUR.QC.${chr}.filtered \
        --freq \
        --out data/freq/1000G.EUR.QC.${chr} \
        --silent
done


# ============================================================
# INLINE: Generate SNP annotations [chr 1-22]
# ============================================================
echo "[5/5] Generating SNP annotations (chr 1-22)..."
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

echo "All inline steps complete. Submitting LSF jobs..."


# ============================================================
# LSF JOB: LD score computation array [1-22]
# No dependency needed — all input files already exist
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
echo "Submitted sldsc_l2[1-22]:  job ID ${JOB_L2}"


# ============================================================
# LSF JOB: Partitioned heritability (single job)
# Depends on: l2[1-22] (all array slots must complete)
# ============================================================
cat > scripts_temp/stage_part_h2.sh << EOF
#!/bin/bash
#BSUB -J sldsc_part_h2
#BSUB -P ${PROJ}
#BSUB -q ${QUEUE}
#BSUB -n 1
#BSUB -R "rusage[mem=16000]"
#BSUB -W 2:00
#BSUB -w "done(${JOB_L2})"
#BSUB -o ${WORKDIR}/logs/o.sldsc_part_h2
#BSUB -e ${WORKDIR}/logs/e.sldsc_part_h2

cd ${WORKDIR}
ml ldsc/1.0.1

ldsc.py \\
    --h2 ibd.sumstats.gz \\
    --ref-ld-chr ${BASELINE_LDSCORES},annot/multipathway.,annot/othergenes. \\
    --frqfile-chr data/freq/1000G.EUR.QC. \\
    --w-ld-chr ${WLDCHR} \\
    --overlap-annot \\
    --print-coefficients \\
    --print-delete-vals \\
    --out results/ibd.multipathway_partitioned_h2
EOF
JOB_H2=$(submit scripts_temp/stage_part_h2.sh)
echo "Submitted sldsc_part_h2:   job ID ${JOB_H2}"


echo ""
echo "============================================================"
echo "Done. Monitor with: bjobs -w"
echo "  sldsc_l2[1-22] (${JOB_L2})  →  sldsc_part_h2 (${JOB_H2})"
echo "============================================================"