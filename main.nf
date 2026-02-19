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
    run_real_prset;
    run_random_sets_prset;
} from './modules/prset/prset.nf'

include {
    calc_empirical_pvalues;
    calc_empirical_pvalues as calc_empirical_pvalues_gsamixer;
    calc_empirical_pvalues as calc_empirical_pvalues_prset;
} from './modules/empirical_pval/empirical_pval.nf'

include {
    tissue_correlation_analysis
} from './modules/tissuespecificity/tissue_correlation.nf'

include {
    opentargets_stats_correlation;
} from './modules/opentargets/opentargets_stats_correlation.nf'

include {
    malacards_correlation;
} from './modules/malacards/malacards_correlation.nf'

include {
    // GSA-MiXeR
    prepare_gsamixer_sumstats;
    split_gsamixer_sumstats;
    gsamixer_plsa_base;
    gsamixer_plsa_full;
    convert_gmt_for_gsamixer;
} from './modules/gsamixer/gsamixer.nf'

include {
    convert_random_gmt_for_gsamixer;
    gsamixer_plsa_full_random;
} from './modules/gsamixer/gsamixer_random.nf'

include {
    calculate_fpr
} from './modules/fpr/fpr_calculation.nf'

include {
    dorothea_correlation
} from './modules/dorothea/dorothea_correlation.nf'

include {
    generate_birewire_random_gmts;
    generate_keeppathsize_random_gmts;
} from './modules/randomization/generate_random_gmts.nf'

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
    def neff = content.neff_col ?: content.n_col  // Use n_col as fallback if neff_col is missing
    
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
        content.summary_statistic_type,          // [11] Summary statistic type (beta/or)
        content.se_col,                          // [12] Standard error column
        neff                                     // [13] Effective sample size (falls back to n_col)
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
    
    //////////////////////////////////////////
    // RANDOMIZATION SETUP - MUST COMPLETE FIRST
    //////////////////////////////////////////
    def birewire_dir
    def keeppathsize_dir
    
    if (params.generate_random_gmts) {
        log.info "========================================="
        log.info "Randomization Setup"
        log.info "========================================="
        log.info "Generating ${params.num_random_sets} randomized GMT files"
        log.info "This MUST complete before random enrichment analyses"
        log.info "========================================="
        
        // Generate BiRewire random GMTs
        log.info "Starting BiRewire randomization..."
        generate_birewire_random_gmts(
            file(params.geneset_real),
            params.num_random_sets
        )
        birewire_dir = generate_birewire_random_gmts.out.birewire_dir
        
        // Generate KeepPathSize random GMTs
        log.info "Starting KeepPathSize randomization..."
        generate_keeppathsize_random_gmts(
            file(params.geneset_real),
            params.num_random_sets
        )
        keeppathsize_dir = generate_keeppathsize_random_gmts.out.keeppathsize_dir
        
        // Wait for both to complete before proceeding
        Channel.empty()
            .mix(
                generate_birewire_random_gmts.out.gmt_files,
                generate_keeppathsize_random_gmts.out.gmt_files
            )
            .collect()
            .subscribe {
                log.info "========================================="
                log.info "Randomization Complete!"
                log.info "BiRewire files: ${birewire_dir}"
                log.info "KeepPathSize files: ${keeppathsize_dir}"
                log.info "========================================="
            }
    } else {
        log.info "Using existing randomized GMT files"
        birewire_dir = params.gmt_dirs.birewire
        keeppathsize_dir = params.gmt_dirs.keeppathsize
    }

    // Create completion signal that waits for BOTH randomization methods
    gmt_ready_signal = birewire_dir
        .combine(keeppathsize_dir)
        .map { bw_dir, kp_dir -> 
            log.info "All random GMT files ready:"
            log.info "  BiRewire: ${bw_dir}"
            log.info "  KeepPathSize: ${kp_dir}"
            return "ready"
        }
    
    // Initialize channels
    all_empirical_inputs = Channel.empty()
    magma_empirical_results = Channel.empty()
    prset_empirical_results = Channel.empty()
    
    //////////////////////////////////////////
    // MAGMA WORKFLOW - WAIT FOR GMTs
    //////////////////////////////////////////
    if (params.run_magma) {
        log.info "Running MAGMA analysis for all traits"
        
        // Extract fields for gene analysis (independent of GMTs)
        magma_gene_data = trait_data.map { trait_tuple ->
            def (trait, gwas_file, rsid_col, chr_col, pos_col, pval_col, n_col) = trait_tuple
            tuple(trait, gwas_file, rsid_col, chr_col, pos_col, pval_col, n_col)
        }
        
        // Run gene analysis (doesn't need GMTs)
        prepared = prepare_input(magma_gene_data)
        def selected_gene_file = params.gene_files[params.background]
        def snp_loc_with_gene_file = prepared.snp_loc_data.map { it + [selected_gene_file] }
        annotated = annotate_genes(snp_loc_with_gene_file)
        gene_analysis_results = run_gene_analysis(annotated.gene_annot_data)
        
        // Combine gene results with randomization methods
        gene_results_with_methods = gene_analysis_results.gene_results
            .combine(Channel.fromList(params.randomization_methods))
            .map { trait, gene_result, rand_method -> 
                [trait, gene_result, rand_method]
            }
        
        // Run real pathway analysis (uses params.geneset_real, not random GMTs)
        real_geneset_results = run_real_geneset(gene_results_with_methods)
        
        // **KEY FIX**: Wait for GMT generation before starting random sets
        perms_ch = Channel.from(1..params.num_random_sets)
        
        random_sets_inputs = gene_results_with_methods
            .combine(gmt_ready_signal)  // Wait for GMTs to be ready
            .map { trait, gene_result, rand_method, ready_signal -> 
                [trait, gene_result, rand_method]
            }
            .combine(perms_ch)
            .map { trait, gene_result, rand_method, perm -> 
                [trait, gene_result, rand_method, perm]
            }
            
        random_sets_results = run_random_sets(random_sets_inputs)
        
        // Group random results by trait and randomization method
        random_results_grouped = random_sets_results
            .groupTuple(by: [0, 2])

        if (params.run_empirical) {
            magma_for_empirical = real_geneset_results.combine(
                random_results_grouped, 
                by: [0, 2]
            ).map { trait, rand_method, real_result_file, random_files_list ->
                def random_dir = "${params.outdir}/magma_random/${rand_method}/${params.background}/${trait}"
                tuple(trait, "magma_${rand_method}", real_result_file, random_dir)
            }
            
            magma_empirical_results = calc_empirical_pvalues(magma_for_empirical)
            all_empirical_inputs = all_empirical_inputs.mix(magma_empirical_results)
        }
    }
    
    //////////////////////////////////////////
    // PRSet WORKFLOW - WAIT FOR GMTs
    //////////////////////////////////////////
    if (params.run_prset) {
        log.info "Running PRSet analysis for all traits"
        
        prset_dedup_data = trait_data
            .filter { trait_tuple -> 
                def trait = trait_tuple[0].toUpperCase()
                trait != "SCZ" && trait != "IBD" && trait != "AD"
            }
            .map { trait_tuple ->
                def (trait, gwas_file, rsid_col, chr_col, pos_col, pval_col, n_col, 
                     binary_target, effect_allele, other_allele, summary_statistic_name, summary_statistic_type) = trait_tuple
                tuple(trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, 
                     summary_statistic_name, summary_statistic_type, "common")
            }
        
        deduplicated_gwas = gwas_remove_dup_snps(prset_dedup_data)
        
        // Real PRSet (doesn't need random GMTs)
        prset_real_inputs = deduplicated_gwas
            .combine(Channel.fromList(params.randomization_methods))
            .map { trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, 
                  summary_statistic_name, summary_statistic_type, _, rand_method ->
                tuple(trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, 
                      summary_statistic_name, summary_statistic_type, rand_method)
            }
    
        real_prset_results = run_real_prset(prset_real_inputs)
        
        // **KEY FIX**: Wait for GMT generation before PRSet random sets
        prset_rand_inputs = deduplicated_gwas
            .combine(gmt_ready_signal)  // Wait for GMTs
            .map { trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, 
                  summary_statistic_name, summary_statistic_type, _, ready_signal ->
                tuple(trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, 
                      summary_statistic_name, summary_statistic_type)
            }
            .combine(Channel.fromList(params.randomization_methods))
            .map { trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, 
                  summary_statistic_name, summary_statistic_type, rand_method ->
                tuple(trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, 
                      summary_statistic_name, summary_statistic_type, rand_method)
            }
        
        perms_ch = Channel.from(1..params.num_random_sets)
        prset_random_inputs = prset_rand_inputs
            .combine(perms_ch)
            .map { trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, 
                  summary_statistic_name, summary_statistic_type, rand_method, perm ->
                tuple(trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, 
                      summary_statistic_name, summary_statistic_type, rand_method, perm)
            }
        
        random_prset_results = run_random_sets_prset(prset_random_inputs)
        
        // Group random results by trait and randomization method
        random_prset_grouped = random_prset_results
            .map { trait, summary_files, rand_method ->
                tuple(trait, rand_method, summary_files)  // Reorder to put rand_method at index 1
            }
            .groupTuple(by: [0, 1])  // Group by trait AND rand_method

        if (params.run_empirical) {
            log.info "Setting up PRSet empirical p-value calculation"
            
            // Combine real results with grouped random results for empirical p-value calculation
            prset_for_empirical = real_prset_for_combine
                .map { trait, summary_file, rand_method ->
                    [trait, rand_method, summary_file]  // Reorder to match random_prset_grouped structure
                }
                .combine(
                    random_prset_grouped, 
                    by: [0, 1]  // Join by trait (index 0) and rand_method (index 1)
                )
                .map { trait, rand_method, summary_file, random_files ->
                    def random_dir = "${params.outdir}/prset_random/${rand_method}/${params.background}/${trait}"
                    tuple(trait, "prset_${rand_method}", summary_file, random_dir)
                }
            
            // Calculate empirical p-values for PRSet
            prset_empirical_results = calc_empirical_pvalues_prset(prset_for_empirical)
            
            // Add to collection for combined results
            all_empirical_inputs = all_empirical_inputs.mix(prset_empirical_results)
        }
    }
    
    //////////////////////////////////////////
    // GSA-MiXeR WORKFLOW - WAIT FOR GMTs
    //////////////////////////////////////////
    if (params.run_gsamixer && params.run_empirical) {
        log.info "Running GSA-MiXeR for random pathways"
        
        // **KEY FIX**: Wait for GMT generation before converting to GSA-MiXeR format
        random_gmt_files = Channel.fromList(
            params.randomization_methods.collect { rand_method ->
                (1..params.num_random_sets).collect { perm ->
                    def gmt_path = "${params.gmt_dirs[rand_method]}/GeneSet.random${perm}.gmt"
                    def gtf_path = params.gtf_reference
                    tuple(rand_method, perm, file(gmt_path), file(gtf_path))
                }
            }.flatten().collate(4)
        ).combine(gmt_ready_signal)  // Wait for GMTs
         .map { rand_method, perm, gmt_file, gtf_file, ready_signal ->
             tuple(rand_method, perm, gmt_file, gtf_file)
         }
        
        random_gmt_converted = convert_random_gmt_for_gsamixer(random_gmt_files)
        
        // For each trait's base model results, combine with ALL random set files
        gsamixer_trait_random_inputs = gsamixer_base
            .combine(random_gmt_converted)
            .map { trait, base_json, base_weights, base_snps, rand_method, perm, baseline_txt, full_gene_txt, full_gene_set_txt ->
                tuple(trait, 
                      base_json,       
                      base_weights,
                      base_snps,    // Include the weights file
                      rand_method, 
                      perm, 
                      baseline_txt,
                      full_gene_txt, 
                      full_gene_set_txt)
            }
        
        // Run GSA-MiXeR full model for each trait-random set combination
        random_gmt_full_results = gsamixer_plsa_full_random(gsamixer_trait_random_inputs)
        
        // Group random results by trait and randomization method
        random_gmt_full_grouped = random_gmt_full_results
            .map { trait, rand_method, perm, json_file, go_test_enrich ->
                tuple(trait, rand_method, go_test_enrich)
            }
            .groupTuple(by: [0, 1])  // Group by trait and rand_method
        
        // Calculate empirical p-values
        if (params.run_empirical) {
            // Combine real results with grouped random results
            gsamixer_for_empirical = gsamixer_full
                .map { trait, full_json, go_test_enrich -> 
                    tuple(trait, go_test_enrich)
                }
                .combine(
                    random_gmt_full_grouped, 
                    by: 0  // Join by trait
                ).map { trait, go_test_enrich, rand_method, random_jsons ->
                    def random_dir = "${params.outdir}/gsamixer_random/${rand_method}/${trait}"
                    tuple(trait, "gsamixer", go_test_enrich, random_dir)
                }
            
            // Calculate empirical p-values for GSA-MiXeR
            gsamixer_empirical_results = calc_empirical_pvalues_gsamixer(gsamixer_for_empirical)
            
            // Add to collection for combined results
            all_empirical_inputs = all_empirical_inputs.mix(gsamixer_empirical_results)
        }
    }
    
    //////////////////////////////////////////
    // FPR CALCULATION WORKFLOW
    //////////////////////////////////////////
    if (params.run_empirical && params.run_fpr_analysis) {
        log.info "Setting up FPR analysis for random pathway results"
        
        // MAGMA FPR Analysis - wait for random results to complete
        if (params.run_magma) {
            log.info "Calculating FPR for MAGMA random results"
            
            // Use the grouped random results from MAGMA workflow
            magma_fpr_inputs = random_results_grouped
                .map { trait, random_files, rand_method ->  // Fix parameter order
                    def random_dir = "${params.outdir}/magma_random/${rand_method}/${params.background}/${trait}"
                    tuple(trait, "magma", rand_method, random_files, random_dir)
                }
            
            magma_fpr_results = calculate_fpr(magma_fpr_inputs)
        }
        
        // PRSet FPR Analysis - wait for random results to complete  
        if (params.run_prset) {
            log.info "Calculating FPR for PRSet random results"
            
            // Use the grouped random results from PRSet workflow
            prset_fpr_inputs = random_prset_grouped
                .map { trait, rand_method, random_files_list ->  
                    // Flatten the list of file lists and filter for .summary files only
                    def summary_files = random_files_list.flatten().findAll { it.toString().endsWith('.summary') }
                    def random_dir = "${params.outdir}/prset_random/${rand_method}/${params.background}/${trait}"
                    tuple(trait, "prset", rand_method, summary_files, random_dir)
                }
            
            prset_fpr_results = calculate_fpr(prset_fpr_inputs)
        }
        
        // GSA-MiXeR FPR Analysis - wait for random results to complete
        if (params.run_gsamixer) {
            log.info "Calculating FPR for GSA-MiXeR random results"
            
            // Use the grouped random results from GSA-MiXeR workflow
            gsamixer_fpr_inputs = random_gmt_full_grouped
                .map { trait, random_files, rand_method ->  // Fix parameter order
                    def random_dir = "${params.outdir}/gsamixer_random/${rand_method}/${trait}"
                    tuple(trait, "gsamixer", rand_method, random_files, random_dir)
                }
            
            gsamixer_fpr_results = calculate_fpr(gsamixer_fpr_inputs)
        }
    }
    
    //////////////////////////////////////////
    // VALIDATION WORKFLOWS
    //////////////////////////////////////////
    
    if (params.run_empirical) {
        log.info "Setting up validation workflows"
        
        // Helper function to group empirical results by trait and base tool
        def group_by_trait_tool = { ch ->
            ch.map { trait, tool, emp_file ->
                // Extract randomization method and base tool
                def base_tool = tool.replaceAll(/_.*$/, '')  // Remove everything after first underscore
                def rand_method = tool.replaceAll(/^.*_/, '') // Get everything after last underscore
                [trait, base_tool, rand_method, emp_file]
            }
            .groupTuple(by: [0, 1]) // Group by trait and base_tool
            .map { trait, base_tool, rand_methods, result_files ->
                // Find indices for each randomization method
                def birewire_idx = rand_methods.findIndexOf { it == 'birewire' }
                def keeppathsize_idx = rand_methods.findIndexOf { it == 'keeppathsize' }
                
                // Only proceed if both randomization methods exist
                if (birewire_idx != -1 && keeppathsize_idx != -1) {
                    [trait, base_tool, result_files[birewire_idx], result_files[keeppathsize_idx]]
                } else {
                    log.warn "Missing randomization method for ${trait} with ${base_tool}"
                    return null
                }
            }
            .filter { it != null }
        }
        
        // Convert malacards_traits string to list for filtering
        def malacards_trait_list = params.malacards_traits.tokenize(',').collect { it.trim().toLowerCase() }
        
        //////////////////////////////////////////
        // OPENTARGETS CORRELATION
        //////////////////////////////////////////
        if (params.run_ot_correlation) {
            log.info "Running OpenTargets correlation analysis"
            
            // MAGMA
            if (params.run_magma) {
                log.info "OpenTargets correlation for MAGMA"
                
                magma_for_opentargets = group_by_trait_tool(magma_empirical_results)
                    .filter { trait, tool, birewire, keeppathsize -> 
                        params.opentargets_supported_traits.contains(trait)
                    }
                
                magma_opentargets_correlation = opentargets_stats_correlation(magma_for_opentargets)
            }
            
            // PRSet
            if (params.run_prset) {
                log.info "OpenTargets correlation for PRSet"
                
                prset_for_opentargets = group_by_trait_tool(prset_empirical_results)
                    .filter { trait, tool, birewire, keeppathsize -> 
                        params.opentargets_supported_traits.contains(trait)
                    }
                
                prset_opentargets_correlation = opentargets_stats_correlation(prset_for_opentargets)
            }
        }
        
        //////////////////////////////////////////
        // TISSUE CORRELATION
        //////////////////////////////////////////
        if (params.run_tissue_correlation) {
            log.info "Running Tissue Specificity correlation analysis"
            
            // MAGMA
            if (params.run_magma) {
                log.info "Tissue correlation for MAGMA"
                
                magma_for_tissue_corr = group_by_trait_tool(magma_empirical_results)
                tissue_correlation_analysis(magma_for_tissue_corr)
            }
            
            // PRSet
            if (params.run_prset) {
                log.info "Tissue correlation for PRSet"
                
                prset_for_tissue_corr = group_by_trait_tool(prset_empirical_results)
                tissue_correlation_analysis(prset_for_tissue_corr)
            }
        }
        
        //////////////////////////////////////////
        // MALACARDS CORRELATION
        //////////////////////////////////////////
        if (params.run_malacards_correlation) {
            log.info "Running MalaCards correlation analysis"
            
            // MAGMA
            if (params.run_magma) {
                log.info "MalaCards correlation for MAGMA"
                
                magma_for_malacards_corr = group_by_trait_tool(magma_empirical_results)
                    .filter { trait, tool, birewire, keeppathsize -> 
                        malacards_trait_list.contains(trait.toLowerCase())
                    }
                
                malacards_correlation(magma_for_malacards_corr)
            }
            
            // PRSet
            if (params.run_prset) {
                log.info "MalaCards correlation for PRSet"
                
                prset_for_malacards_corr = group_by_trait_tool(prset_empirical_results)
                    .filter { trait, tool, birewire, keeppathsize -> 
                        malacards_trait_list.contains(trait.toLowerCase())
                    }
                
                malacards_correlation(prset_for_malacards_corr)
            }
        }
        
        //////////////////////////////////////////
        // DOROTHEA CORRELATION
        //////////////////////////////////////////
        if (params.run_dorothea_correlation) {
            log.info "Running DoRothEA correlation analysis"
            
            // MAGMA
            if (params.run_magma) {
                log.info "DoRothEA correlation for MAGMA"
                
                magma_for_dorothea_corr = group_by_trait_tool(magma_empirical_results)
                dorothea_correlation(magma_for_dorothea_corr)
            }
            
            // PRSet
            if (params.run_prset) {
                log.info "DoRothEA correlation for PRSet"
                
                prset_for_dorothea_corr = group_by_trait_tool(prset_empirical_results)
                dorothea_correlation(prset_for_dorothea_corr)
            }
        }
    }
}