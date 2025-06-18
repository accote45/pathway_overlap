nextflow.enable.dsl=2

////////////////////////////////////////////////////////////////////
//                  Module inclusion
////////////////////////////////////////////////////////////////////

include {
    prepare_input;
    annotate_genes;
    run_gene_analysis;
    run_real_geneset;
    run_random_sets;
} from './modules/magma/magma.nf'

include { 
    gwas_remove_dup_snps;
    run_random_sets_prset;
} from './modules/prset/prset.nf'

include {
    calc_empirical_pvalues;
    combine_empirical_results;
} from './modules/empirical_pval/empirical_pval.nf'

////////////////////////////////////////////////////////////////////
//                  Setup Channels
////////////////////////////////////////////////////////////////////

import groovy.json.JsonSlurper

// Parse the JSON input into a comprehensive Channel
def jsonSlurper = new JsonSlurper()
def configFile = new File(params.traits_config)
def phenoConfig = jsonSlurper.parseText(configFile.text)

// Create a single comprehensive channel with all trait information
trait_data = Channel.fromList(phenoConfig.collect { content ->
    tuple(
        content.trait,                           // [0] Trait name
        content.gwas_file,                       // [1] GWAS file
        content.rsid_col,                        // [2] RSID column
        content.chr_col,                         // [3] Chromosome column
        content.pos_col,                         // [4] Position column
        content.pval_col,                        // [5] P-value column
        content.n_col,                           // [6] Sample size column
        content.binary_target,                   // [7] Binary trait flag
        content.effect_allele,                   // [8] Effect allele column
        content.other_allele,                    // [9] Other allele column
        content.summary_statistic_name,          // [10] Summary statistic column name
        content.summary_statistic_type           // [11] Summary statistic type (beta/or)
    )
})

////////////////////////////////////////////////////////////////////
//                  Subworkflow: MAGMA analysis
////////////////////////////////////////////////////////////////////

workflow magma {
    take:
    trait_randomization_data  // Now includes randomization method
    
    main:
    // Extract fields needed for MAGMA from trait_randomization_data
    magma_data = trait_randomization_data.map { fullData ->
        def (trait, gwas_file, rsid_col, chr_col, pos_col, pval_col, n_col, binary_target, effect_allele, other_allele, summary_statistic_name, summary_statistic_type, rand_method) = fullData
        tuple(trait, gwas_file, rsid_col, chr_col, pos_col, pval_col, n_col, rand_method)
    }
    
    prepared = prepare_input(magma_data)
    
    // Select gene file based on background setting
    def selected_gene_file = params.gene_files[params.background]

    // Append gene_file path to each tuple
    def snp_loc_with_gene_file = prepared.snp_loc_data.map { it + [selected_gene_file] }

    annotated = annotate_genes(snp_loc_with_gene_file)
    gene_analysis = run_gene_analysis(annotated.gene_annot_data)
    run_real_geneset(gene_analysis.gene_results)
    
    // Create permutations for each trait/randomization combination
    perms_ch = Channel.from(1..params.num_random_sets)
    random_sets_inputs = gene_analysis.gene_results
        .combine(perms_ch)
        .map { trait, gene_result, rand_method, perm -> 
            println "Creating MAGMA job for ${trait} with ${rand_method} permutation ${perm}"
            [trait, gene_result, rand_method, perm] 
        }

    run_random_sets(random_sets_inputs)
    
    // Group completed results by trait and randomization method
    completed_magma = gene_analysis.gene_results
        .map { trait, gene_result, rand_method ->
            [trait, gene_result, rand_method]
        }
    
    emit:
    gene_results = completed_magma
}

////////////////////////////////////////////////////////////////////
//                  Subworkflow: PRSET analysis
////////////////////////////////////////////////////////////////////

workflow prset {
    take:
    trait_randomization_data

    main:
    // Map the comprehensive data to the format needed for PRSet
    prset_data = trait_randomization_data.map { fullData ->
        def (trait, gwas_file, rsid_col, chr_col, pos_col, pval_col, n_col, binary_target, effect_allele, other_allele, summary_statistic_name, summary_statistic_type, rand_method) = fullData
        tuple(trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, summary_statistic_name, summary_statistic_type, rand_method)
    }
    
    // First, remove duplicate SNPs from GWAS files
    deduplicated_gwas = gwas_remove_dup_snps(prset_data)
    
    perms_ch = Channel.from(1..params.num_random_sets)
    
    // Combine prset data with permutation numbers
    prset_random_inputs = deduplicated_gwas
        .combine(perms_ch)
        .map { trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, summary_statistic_name, summary_statistic_type, rand_method, perm ->
            println "Creating PRSET job for ${trait} with ${rand_method} permutation ${perm}"
            [trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, summary_statistic_name, summary_statistic_type, rand_method, perm]
        }
    
    // Run PRSet for random sets
    run_random_sets_prset(prset_random_inputs)
    
    // Create completed channel for empirical calculation
    completed_prset = deduplicated_gwas
        .map { trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, summary_statistic_name, summary_statistic_type, rand_method ->
            [trait, gwas_file, rand_method]
        }
    
    emit:
    summary_results = completed_prset
}

////////////////////////////////////////////////////////////////////
//                  Main Workflow
////////////////////////////////////////////////////////////////////

workflow {
    log.info "Pipeline configuration:"
    log.info "  Run MAGMA: ${params.run_magma}"
    log.info "  Run PRSet: ${params.run_prset}"
    log.info "  Calculate Empirical P-values: ${params.run_empirical}"
    log.info "  Randomization methods: ${params.randomization_methods}"
    
    // Create randomization method channel
    randomization_channel = Channel.fromList(params.randomization_methods)
    
    // Combine trait data with randomization methods
    trait_randomization = trait_data.combine(randomization_channel)
        .map { trait_data_tuple, rand_method ->
            // Append rand_method to the existing tuple
            trait_data_tuple + rand_method
        }
    
    // Store all tool results for empirical p-value calculation
    all_empirical_inputs = Channel.empty()
    
    //////////////////////////////////////////
    // MAGMA WORKFLOW + EMPIRICAL P-VALUES
    //////////////////////////////////////////
    if (params.run_magma) {
        log.info "Running MAGMA analysis for all randomization methods"
        
        // Run MAGMA for each trait/randomization combination
        magma_results = magma(trait_randomization)
        
        // Calculate empirical p-values for MAGMA immediately after completion
        if (params.run_empirical) {
            log.info "Setting up MAGMA empirical p-value calculation"
            
            // Prepare MAGMA results for empirical calculation
            magma_for_empirical = magma_results.gene_results
                .map { trait, gene_result, rand_method ->
                    def real_file = "${params.outdir}/magma_real/${trait}/${trait}_real_set.gsa.out"
                    def random_dir = "${params.outdir}/magma_random/${rand_method}/${params.background}/${trait}"
                    [trait, "magma_${rand_method}", file(real_file), random_dir]
                }
            
            // Calculate empirical p-values for MAGMA
            magma_empirical = calc_empirical_pvalues(magma_for_empirical)
            
            // Add to collection
            all_empirical_inputs = all_empirical_inputs.mix(
                magma_empirical.map { trait, tool, result_file -> result_file }
            )
        }
    }
    
    //////////////////////////////////////////
    // PRSET WORKFLOW + EMPIRICAL P-VALUES  
    //////////////////////////////////////////
    if (params.run_prset) {
        log.info "Running PRSet analysis for all randomization methods"
        
        // Run PRSet for each trait/randomization combination
        prset_results = prset(trait_randomization)
        
        // Calculate empirical p-values for PRSet immediately after completion
        if (params.run_empirical) {
            log.info "Setting up PRSet empirical p-value calculation"
            
            prset_for_empirical = prset_results.summary_results
                .map { trait, summary_result, rand_method ->
                    def real_file = "${params.outdir}/prset_real/${trait}/${trait}_real.summary"
                    def random_dir = "${params.outdir}/prset_random/${rand_method}/${params.background}/${trait}"
                    [trait, "prset_${rand_method}", file(real_file), random_dir]
                }
            
            prset_empirical = calc_empirical_pvalues(prset_for_empirical)
            
            all_empirical_inputs = all_empirical_inputs.mix(
                prset_empirical.map { trait, tool, result_file -> result_file }
            )
        }
    }
    
    //////////////////////////////////////////
    // COMBINE ALL EMPIRICAL RESULTS
    //////////////////////////////////////////
    if (params.run_empirical) {
        log.info "Setting up empirical results combination"
        
        // Collect all empirical inputs and handle empty case
        combined_empirical = combine_empirical_results(
            all_empirical_inputs.collect().ifEmpty([])
        )
    }
}


