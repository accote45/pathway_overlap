# standard GSA-MiXeR pipeline begins here - change only if you know what you're doing

ml singularity
ml python

export MIX_FOLDER="/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/gsa_mixer/"
export MIXER_SIF=${MIX_FOLDER}/"gsa-mixer.sif"
export MIXER_PY="singularity exec --home pwd:/home ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"
export REFERENCE_FOLDER=${MIX_FOLDER}/reference
export BIM_FILE=${MIX_FOLDER}/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim
export SUMSTATS_FOLDER="/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/gwas"
export LOADLIB_FILE=${MIX_FOLDER}/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bin
export ANNOT_FILE=${MIX_FOLDER}/1000G_EUR_Phase3_plink/baseline_v2.2_1000G.EUR.QC.@.annot.gz

# split sumstats by chr
${MIXER_PY} split_sumstats \
    --trait1-file ${SUMSTATS_FOLDER}/${SUMSTATS_FILE}.sumstats.gz
    --out ${SUMSTATS_FOLDER}/${SUMSTATS_FILE}.chr@.sumstats.gz


# compute LD scores to be used for SNP step
${MIXER_PY} ld \
--bfile 1000G.EUR.QC.@ \
--out 1000G.EUR.QC.@.ld \
--r2min 0.05 \
--ldscore-r2min 0.01 \
--ld-window-kb 10000


# generate bin files from LD reference
#.bin file for --load.lib
${MIXER_PY} plsa \
--bim-file 1000G.EUR.QC.@.bim \
--ld-file 1000G.EUR.QC.@.ld \
--use-complete-tag-indices \
--savelib-file 1000G.EUR.QC.@.bin \
--out 1000G.EUR.QC.@


# input gwas header (lower or upper case)
SNP/RSID, CHR, BP/POS, A1/EffectAllele, A2/OtherAllele, N (effect sample size), Z



# split sumstats per chromosome
  ${MIXER_PY} split_sumstats \
      --trait1-file ${SUMSTATS_FOLDER}/Alkaline_phosphatase_gwas.txt \
      --out ${SUMSTATS_FOLDER}/Alkaline_phosphatase.chr@.sumstats.gz


# GSA-MiXeR analysis - baseline
${MIXER_PY} plsa --gsa-base \
        --trait1-file ${SUMSTATS_FOLDER}/Alkaline_phosphatase.chr@.sumstats.gz \
        --out ${OUT_FOLDER}/Alkaline_phosphatase_base \
        --bim-file ${BIM_FILE} --use-complete-tag-indices --loadlib-file ${LOADLIB_FILE} \
        --go-file ${REFERENCE_FOLDER}/gsa-mixer-baseline-annot_10mar2023.csv \
        --annot-file ${ANNOT_FILE} \
        ${EXTRA_FLAGS}

# GSA-MiXeR analysis - enrichment model
${MIXER_PY} plsa --gsa-full \
        --trait1-file ${SUMSTATS_FOLDER}/${SUMSTATS_FILE}.chr@.sumstats.gz \
        --out ${OUT_FOLDER}/${SUMSTATS_FILE}_full \
        --bim-file ${BIM_FILE} --use-complete-tag-indices --loadlib-file ${LOADLIB_FILE} \
        --go-file ${REFERENCE_FOLDER}/gsa-mixer-gene-annot_10mar2023.csv \
        --go-file-test ${REFERENCE_FOLDER}/gsa-mixer-hybridLOO-annot_10mar2023.csv \
        --annot-file ${ANNOT_FILE} \
        --load-params-file ${OUT_FOLDER}/${SUMSTATS_FILE}_base.json \
        ${EXTRA_FLAGS}