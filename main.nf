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
//                  Main Workflow
////////////////////////////////////////////////////////////////////

workflow {
    log.info "Pipeline configuration:"
    log.info "  Run MAGMA: ${params.run_magma}"
    log.info "  Run PRSet: ${params.run_prset}"
    log.info "  Calculate Empirical P-values: ${params.run_empirical}"
    log.info "  Randomization methods: ${params.randomization_methods}"
    
    // Store all tool results for empirical p-value calculation
    all_empirical_inputs = Channel.empty()
    
    //////////////////////////////////////////
    // MAGMA WORKFLOW (OPTIMIZED)
    //////////////////////////////////////////
    if (params.run_magma) {
        log.info "Running MAGMA analysis for all traits"
        
        // Extract fields needed for MAGMA gene analysis (once per trait)
        magma_gene_data = trait_data.map { trait_tuple ->
            def (trait, gwas_file, rsid_col, chr_col, pos_col, pval_col, n_col) = trait_tuple
            tuple(trait, gwas_file, rsid_col, chr_col, pos_col, pval_col, n_col)
        }
        
        // Run gene analysis steps ONCE per trait
        prepared = prepare_input(magma_gene_data)
        
        def selected_gene_file = params.gene_files[params.background]
        def snp_loc_with_gene_file = prepared.snp_loc_data.map { it + [selected_gene_file] }
        
        annotated = annotate_genes(snp_loc_with_gene_file)
        gene_analysis_results = run_gene_analysis(annotated.gene_annot_data)
        
        // Now we have gene results for each trait
        // Combine them with randomization methods for pathway analysis
        gene_results_with_methods = gene_analysis_results.gene_results
            .combine(Channel.fromList(params.randomization_methods))
            .map { trait, gene_result, rand_method -> 
                [trait, gene_result, rand_method]
            }
        
        // Run pathway analysis for real data (once per trait-randomization pair)
        real_geneset_results = run_real_geneset(gene_results_with_methods)
        
        // Run random sets analysis (once per trait-randomization-permutation)
        perms_ch = Channel.from(1..params.num_random_sets)
        random_sets_inputs = gene_results_with_methods
            .combine(perms_ch)
            .map { trait, gene_result, rand_method, perm -> 
                [trait, gene_result, rand_method, perm]
            }
            
        random_sets_results = run_random_sets(random_sets_inputs)
        
        // Group random results by trait and randomization method
        random_results_grouped = random_sets_results
            .groupTuple(by: [0, 2])  // Group by trait and rand_method

        // Calculate empirical p-values - now with explicit dependency on ALL random results
        if (params.run_empirical) {
            log.info "Setting up MAGMA empirical p-value calculation for each randomization method"
            
            // Combine real results with grouped random results
            magma_for_empirical = real_geneset_results
                .join(random_results_grouped)
                .map { tuple -> 
                    // Explicitly extract each element to avoid issues with destructuring
                    def trait = tuple[0]
                    def real_file = tuple[1]
                    def rand_method = tuple[2]
                    def random_files = tuple[3]  // List of random files
                    
                    def random_dir = "${params.outdir}/magma_random/${rand_method}/${params.background}/${trait}"
                    [trait, "magma_${rand_method}", real_file, random_dir]
                }
            
            magma_empirical = calc_empirical_pvalues(magma_for_empirical)
            
            all_empirical_inputs = all_empirical_inputs.mix(
                magma_empirical.map { trait, tool, result_file -> result_file }
            )
        }
    }
    
    //////////////////////////////////////////
    // PRSET WORKFLOW (OPTIMIZED)
    //////////////////////////////////////////
    if (params.run_prset) {
        log.info "Running PRSet analysis for all traits"
        
        // Extract PRSet data without randomization (for deduplication)
        prset_dedup_data = trait_data.map { trait_tuple ->
            def (trait, gwas_file, rsid_col, chr_col, pos_col, pval_col, n_col, 
                 binary_target, effect_allele, other_allele, summary_statistic_name, summary_statistic_type) = trait_tuple
            tuple(trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, 
                 summary_statistic_name, summary_statistic_type)
        }
        
        // Run SNP deduplication ONCE per trait
        deduplicated_gwas = gwas_remove_dup_snps(
            prset_dedup_data.map { it -> it + ["none"] }  // Add placeholder for rand_method
        )
        
        // Now combine with randomization methods
        prset_rand_inputs = deduplicated_gwas
            .map { it -> 
                // Remove the placeholder rand_method  
                def items = it.toList()
                items.removeLast()
                tuple(*items)
            }
            .combine(Channel.fromList(params.randomization_methods))
            .map { trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, 
                  summary_statistic_name, summary_statistic_type, rand_method ->
                tuple(trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, 
                     summary_statistic_name, summary_statistic_type, rand_method)
            }
        
        // Create random set inputs (one per permutation)
        perms_ch = Channel.from(1..params.num_random_sets)
        prset_random_inputs = prset_rand_inputs
            .combine(perms_ch)
            .map { trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, 
                  summary_statistic_name, summary_statistic_type, rand_method, perm ->
                [trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, 
                 summary_statistic_name, summary_statistic_type, rand_method, perm]
            }
        
        // Run PRSet for random sets
        random_prset_results = run_random_sets_prset(prset_random_inputs)
        
        // Group random results by trait and randomization method
        random_prset_grouped = random_prset_results
            .groupTuple(by: [0, 2])  // Group by trait and rand_method

        // Calculate empirical p-values with explicit dependency
        if (params.run_empirical) {
            log.info "Setting up PRSet empirical p-value calculation"
            
            prset_for_empirical = prset_rand_inputs
                .map { tuple -> 
                    def trait = tuple[0]
                    def rand_method = tuple[9]
                    [trait, rand_method]  // Create key for joining
                }
                .join(random_prset_grouped)
                .map { tuple ->
                    def trait = tuple[0]
                    def rand_method = tuple[1]
                    def random_files = tuple[2]  // List of random files
                    
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


