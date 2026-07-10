process genotype_filtering{
    // preprocess the 1000G genotypes to ensure there are no duplicated SNPs
    // Also calculate the allele frequency, which is required by LDSC
    publishDir "data/reference/freq", mode: 'copy', pattern: "*frq"
    label 'normal'
    input:
        tuple   val(chr),
                path(bed),
                path(bim),
                path(fam),
                path(pop)
    output:
        tuple   val(chr),
                path("chr${chr}.bed"),
                path("chr${chr}.bim"),
                path("chr${chr}.fam"), emit: geno
        tuple   val(chr),
                path("chr${chr}.frq"), emit: freq
    script:
    name=bed.baseName
    """
    awk '\$2 in a {print \$2; next} !(\$2 in a){a[\$2]=1; next} ' ${bim} > dup.snp
    wc -l dup.snp
    plink \
        --bfile ${name} \
        --mac 5 \
        --geno 0.01 \
        --keep ${pop} \
        --make-bed \
        --out chr${chr} \
        --exclude dup.snp \
        --freq
    """
}

process annotate_baseline_genome{
    // This is to generate the "baseline" LDSC data set that are 
    // recommended by LDSC 
    label 'normal'
    input: 
        tuple   val(chr), 
                path(bed), 
                path(bim), 
                path(fam),
                path("*")
    output:
        tuple   val(chr), 
                path("baseline.snp")
    script:
    prefix=bed.baseName
    """
    beds=""
    files=`ls *bed | grep -v ${bed}`
    function join_by { local IFS="\$1"; shift; echo "\$*"; }
    beds=`join_by , \${files}`
    annotator \
        --target ${prefix} \
        --out baseline \
        --bed \${beds}
    """
}
process trim_annotation{
    
    // Modify the pathway membership SNP file so that we now have roughly 250 gene sets
    // per file, which allow further parallelization for LD score calculation
    label 'more_mem'
    input:
        tuple   val(chr),
                path(annot)
    output:
        tuple   val(chr),
                path("${name}-*snp")
    script:
    name=annot.baseName
    """
    #!/usr/bin/env Rscript
    library(magrittr)
    library(data.table)
    # To speed up LDSC calculation, we limit the number of pathway to 250
    dat <- fread("${annot}")[,-c("Base")]
    setnames(dat, "Background", "Control")
    num <- ncol(dat)-4
    if(num <= 250){
        fwrite(dat, "${name}-0.snp", sep="\\t")
    }else{
        range <- c(seq(4, ncol(dat), 250), ncol(dat)) %>%
            unique
        for(i in 1:(length(range)-1)){
            start <- range[i]+1
            end <- range[i+1]
            dat[,c("CHR","BP","SNP","CM", colnames(dat)[start:end]), with=F] %>%
                fwrite(., paste0("${name}-", i, ".snp"), sep="\\t")
        }
    }
    """
}
process annotate_genome{
    // annotate the genome using the LDSC annotator. This will generate 
    // the required input for LDSC to calculate the set based LD score. 
    label 'normal'
    input: 
        tuple   val(chr), 
                path(bed), 
                path(bim), 
                path(fam), 
                path(gtf), 
                path(gmt), 
                val(type),
                val(wind5),
                val(wind3),
                val(xregion),
                val(out)
    output:
        tuple   val(chr), 
                path("${out}.snp")
    script:
        prefix=bed.baseName
        """
        wind5_com=""
        wind3_com=""
        if [ ! -z "${wind5}" ]; then
            wind5_com="--wind-5 "${wind5}kb
        fi
        if [ ! -z "${wind3}" ]; then
            wind3_com="--wind-3 "${wind3}kb
        fi
        xrange=""
        if [ ! -z "${xregion}" ]; then
            xrange="--x-range ${xregion}"
        fi
        celltype=""
        if [ "${type}" == "specificity" ]; then
            celltype="--celltype"
        fi
        annotator \
            --target ${prefix} \
            --out ${out} \
            --gtf ${gtf} \
            --msigdb ${gmt} \
            \${xrange} \
            \${wind5_com} \
            \${wind3_com}
        """
}


process calculate_ldscore{
    label 'ldsc_calculation'
    input:
        tuple   val(chr), 
                path(annot), 
                path(bed), 
                path(bim), 
                path(fam),  
                path(hm)
    output:
        tuple   path("${chr}.l2.ldscore.gz"), 
                path("${chr}.l2.M_5_50"), 
                path("${chr}.l2.M"), 
                path (annot)
    script:
        base=bed.baseName
        """
        ldsc.py \
            --l2 \
            --bfile ${base} \
            --ld-wind-kb 1000 \
            --annot ${annot} \
            --out ${chr} \
            --print-snp ${hm}
        """
}

process rename_baseline_files{
    // rename the baseline & weight file so that it follows the same format
    // as other files for better downstream integration
    label 'tiny'
    input:
        tuple   path(ldscore), 
                path(m50), 
                path(m), 
                path (annot)
    output:
        tuple   path("${annot_name}.l2.ldscore.gz"), 
                path("${annot_name}.l2.M_5_50"), 
                path("${annot_name}.l2.M"), 
                path (annot)

    script:
    annot_name=annot.baseName
    """
    cp ${ldscore} ${annot_name}.l2.ldscore.gz
    cp ${m50} ${annot_name}.l2.M_5_50
    cp ${m} ${annot_name}.l2.M
    """
}
process unwrap_ldscore{
    // Each LD score file now contained LD scores for around 250 gene sets. 
    // Problem with partition LD score is that it can easily generate singular
    // matrix, leading to failed analyses. To avoid this problem, we want 
    // to run the analysis *one* gene set at a time, therefore minimizing risk
    // of failed anlaysis
    label 'more_mem'
    input:
        path("*")
    output:
        tuple   path ("*l2.ldscore.gz"), 
                path ("*.l2.M*"), 
                path("*annot.gz")
    script:
        """
        #!/usr/bin/env Rscript
        library(magrittr)
        library(data.table)
        ldsc <- list.files(pattern="ldscore.gz")
        m <- list.files(pattern="M\$")
        m50 <- list.files(pattern="M_5_50\$")
        annot <- list.files(pattern="snp\$")

        clean_up <- function(i) {
            gsub("L2\$", "", i)
        }
        gene.sets <- ldsc %>% 
                head(n=1) %>%
                fread(.) %>% 
                .[,!c("CHR", "SNP", "BP")] %>%
                colnames 
        gene.sets.name <- gene.sets %>%
                clean_up
        for(i in ldsc){
            input <- fread(i)
            suffix <- strsplit(i, split="-") %>% unlist %>% tail(n=1)
            for(j in 1:length(gene.sets)){
                input[,c("CHR", "SNP", "BP", gene.sets[j]), with=F] %>%
                    setnames(., gene.sets[j], gene.sets.name[j]) %>%
                    fwrite(., paste(gene.sets.name[j],suffix,sep="-"), sep="\\t")
            }
        }
        print_annot_files <- function(x, suffix, gene.sets.name){
            for(i in x){
                input <- fread(i)
                for(j in gene.sets.name){
                    fwrite(input[,c("CHR", "BP","SNP", "CM", j), with=F], paste(j, suffix, sep="-"), sep="\\t")
                }
            }
        }
        print_m_files <- function(x, gene.sets.name){
            for(i in x){
                input <- fread(i)
                suffix <- strsplit(i, split="-") %>% unlist %>% tail(n=1)
                for(j in 1:length(gene.sets.name)){
                    fwrite(input[,..j], paste(gene.sets.name[j], suffix, sep="-"), col.names=F, row.names=F)
                }
            }
        }
        print_m_files(m, gene.sets.name)
        print_m_files(m50, gene.sets.name)
        suffix <- ldsc %>% 
            head(n=1) %>%
            gsub(".l2.ldscore.gz", ".annot.gz", .)
        print_annot_files(annot, suffix, gene.sets.name)
        """
}

process zip_ldsc{
    // Given the number of files generated for each gene set, it is likely that
    // we will run into limits with pipe limit in downstream analysis. (e.g. cannot 
    // even `ls` within the directory). To avoid this problem, and to reduce storage 
    // space requirement, we will zip all LD score files
    label 'normal'
    publishDir "data/reference/ldsc", mode: 'move'
    input:
        tuple   val(name), 
                path("*")
    output:
        path("${name}.zip")
    script:
    """
    zip ${name}.zip *
    """
}


process munge_sumstats{
    label 'normal'
    input:
        tuple   val(meta),
                path(sumstat)
    output:
        tuple   val(meta),
                path("${name}.sumstats.gz")
    script:
    name = meta instanceof Map ? meta.name : meta
    """ 
    munge_sumstats.py \
        --signed-sumstats BETA,0 \
        --a1 A1 \
        --a2 A2 \
        --snp SNP \
        --p P \
        --sumstats ${sumstat} \
        --out ${name} \
        --chunksize 500000 \
        --N-col N
    """
}

process obtain_baseline_bed{
    // Download the baseline bed file form the LDSC repository
    // TODO: Maybe host this data somewhere within our container to make sure it will always
    // work
    executor 'local'
    label 'tiny'
    output:
        path("*.bed")
    script:
    """
    wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/baselineLD_v2.2_bedfiles.tgz
    tar -xvf baselineLD_v2.2_bedfiles.tgz
    mv bed/*bed .
    """
}

process update_baseline{
    // This process modify the baseline annotation output to better fit requirement of 
    // LDSC. It will also update the weight annotation so that the MHC region is removed
    label 'normal'
    input:
        tuple   val(chr),
                path(annot)
    output:
        tuple   val(chr),
                path("baseline-${chr}.annot"), emit: baseline
        tuple   val(chr),
                path("weight-${chr}.annot"), emit: weight
    script:
    """
    #!/usr/bin/env Rscript
    library(magrittr)
    library(data.table)
    annotate <- fread("${annot}")[, -c("Background")]
    # need to replace the baseline control with a vector of 1
    setnames(annotate, "Base", "baseline")
    fwrite(annotate, "baseline-${chr}.annot", sep="\\t")
    # For weight, we remove MHC as suggested by LDSC
    annotate[, c("CHR", "BP", "SNP", "CM", "baseline")] %>%
        .[(CHR==6 & BP >= 25000000 & BP <= 34000000), baseline := 0] %>%
        fwrite(., "weight-${chr}.annot", sep="\\t")
    """
}



process ldsc_biochemical_analysis{
    // Will generate the lite.results if the run is successful
    label 'normal'
    afterScript 'ls * | grep -v lite.results | xargs rm'
    maxForks '500'
    input:
        tuple   val(pheno), 
                path(sumstat),
                path(biochem), 
                path(baseline), // baseline
                path(weight), // weight
                val(freq),
                path("*")  // freq
    output:
        tuple   val(pheno),
                val("Set"),
                path("${name}.lite.results") optional true
    script:
    name=biochem.baseName
    """
    freq=${freq}
    unzip ${biochem}
    unzip ${baseline}
    unzip ${weight}
    ldsc.py \
        --h2 ${sumstat} \
        --ref-ld-chr baseline-,${name}- \
        --w-ld-chr weight- \
        --frqfile-chr \${freq/.frq} \
        --out ${name} \
        --overlap-annot && \
        sed -e 1b -e '\$!d' ${name}.results |\
        awk 'NR==1{print "Name\\tCoefficient\\tCoefficient_std_error\\tCoefficient_P_value"} NR!=1{gsub("_1\$","",\$1);print \$1"\\t"\$5"\\t"\$6"\\t"\$7}' > ${name}.lite.results \
        || echo "true"
    """
}

process ldsc_specificity_analysis{
    label 'long'
    afterScript 'ls * | grep -v cell_type_results.txt | xargs rm'
    maxForks '100'
    input:
        tuple   val(meta),
                val(sumstat),
                path(specificity),
                path(control),
                path(baseline),
                path(weight)
    output:
        tuple   val(meta),
                val("Specificity"),
                path("${name}.cell_type_results.txt")
    script:
    name=specificity.baseName
    controlName = control.baseName
    """
    unzip ${specificity} 
    unzip ${control}
    unzip ${baseline}
    unzip ${weight}
    echo -e "${name}\\t${name}-,${controlName}-" > ${name}.ldcts
    ldsc.py \
        --h2-cts ${sumstat} \
        --ref-ld-chr baseline- \
        --w-ld-chr weight- \
        --out ${name} \
        --ref-ld-chr-cts ${name}.ldcts
    """
}

process modify_ldsc_result{
    label 'normal'
    input:
        tuple   val(meta),
                val(type),
                path(ldsc)
    output:
        tuple   val(meta),
                val(type),
                path("${name}-${type}-ldsc.csv")
    script:
    name = meta instanceof Map? meta.name: meta
    """
    #!/usr/bin/env Rscript
    library(data.table)
    library(magrittr)
    files <- list.files()
    results <- NULL
    for(i in files){
        results %<>% rbind(., fread(i))
    }
    results %>%
        setnames(., "Name", "Set") %>%
        setnames(., "Coefficient_P_value", "P") %>%
        .[, c("Set", "Coefficient", "P")] %>%
        .[, Software := "LDSC"] %>%
        .[, Type := "${type}"] %>%
        .[, Trait := "${name}"] %>%
        fwrite(., "${name}-${type}-ldsc.csv")
    """
}


process ldsc_group_biochemical_analysis{
    label 'group_ldsc'
    // main reason we need a separated process for this analysis is because when
    // more than one permutation is taken, we will submit 4000 * number of permutation
    // * number of parameter combinations jobs to server, which can easily cause problem.
    // As a result of that, we group multiple pathways in one submission, thus reduce amount 
    // of processes. Otherwise, we would ideally reuse the ldsc_biochemical_analysis processes
    afterScript 'ls * | grep -v lite.results | xargs rm'
    // need to restrict it so that we don't do too much at a same time and hoard all 
    // server resources
    input:
        tuple   val(meta),
                path(sumstat),
                path(baseline), // baseline
                path(weight), // weight
                val(freq),
                path("*"),  // freq
                path(biochem)
    output:
        tuple   val(meta),
                path("*.lite.results")
    script:
    """
    freq=${freq}
    unzip ${baseline}
    unzip ${weight}
    sets=( ${biochem} )
    # Use name of the first entry to avoid name duplication across multiple sets
    name=\${sets[0]/.zip/}
    echo -e "Name\\tCoefficient\\tCoefficient_std_error\\tCoefficient_P_value" > \${name}.lite.results
    for i in \${sets[@]}; 
    do
        unzip \${i}
        {
            ldsc.py \
            --h2 ${sumstat} \
            --ref-ld-chr baseline-,\${i/.zip/}- \
            --w-ld-chr weight- \
            --frqfile-chr \${freq/.frq} \
            --out \${i/.zip/} \
            --overlap-annot && 
            sed -e 1b -e '\$!d' \${i/.zip/}.results | awk 'NR != 1 {gsub("_1\$","",\$1);print \$1"\\t"\$5"\\t"\$6"\\t"\$7}' >> \${name}.lite.results
            } ||
            echo "Error"
        rm \${i/.zip/}-*
    done
    """
}

process clean_freq{
    label 'tiny'
    afterScript "rm ${name}*"
    input:
        tuple   val(chr),
                path(freq),
                path(baseline)
    output:
        tuple   val("chr@-clean.frq"),
                path("chr${chr}-clean.frq")
    script:
    name=baseline.baseName
    """
    unzip ${baseline}
    awk 'NR==FNR{a[\$3]="a"} NR != FNR && (FNR==1 || \$2 in a){print}' ${name}-${chr}.annot ${freq} > chr${chr}-clean.frq
    """
}

process modify_simulated_ldsc_result{
    label 'normal'
    input:
        tuple   val(meta),
                path("*")
    output:
        tuple   val(meta),
                path("*-ldsc.csv")
    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    library(magrittr)
    files <- list.files()
    results <- NULL
    for(i in files){
        results %<>% rbind(., fread(i))
    }
    results %>%
        setnames(., "Name", "Set") %>%
        setnames(., "Coefficient_P_value", "P") %>%
        .[, c("Set", "Coefficient", "P")] %>%
        .[, Software := "LDSC"] %>%
        .[, Perm := "${meta.name}"] %>%
        .[, NumSet := "${meta.numSet}"] %>%
        .[, Heritability := "${meta.herit}"] %>%
        .[, TargetSize := "${meta.targetSize}"] %>%
        .[, BaseSize := "${meta.baseSize}"] %>%
        fwrite(., "${meta.name}-${meta.numSet}-${meta.herit}-${meta.targetSize}-${meta.baseSize}-ldsc.csv")
    """
}
