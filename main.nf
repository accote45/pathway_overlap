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
    delta_rank_ot_correlation
} from './modules/opentargets/delta_rank_ot.nf'

include {
    delta_rank_malacards_correlation
} from './modules/malacards/delta_rank_malacards.nf'

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
    log.info "  Randomization methods: ${params.randomization_methods}"
    log.info "  Run Tissue Specificity Analysis: ${params.run_tissue_correlation}"
    
    // Initialize empty channels
    all_empirical_inputs = Channel.empty()
    magma_empirical_results = Channel.empty()
    prset_empirical_results = Channel.empty()
    
    //////////////////////////////////////////
    // MAGMA WORKFLOW
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
        
        // Combine gene results with randomization methods for pathway analysis
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

        // Calculate empirical p-values - with explicit dependency on ALL random results
        if (params.run_empirical) {
            log.info "Setting up MAGMA empirical p-value calculation for each randomization method"
            
            // Combine real results with grouped random results for empirical p-value calculation
            magma_for_empirical = real_geneset_results.combine(
                random_results_grouped, 
                by: [0, 2]  // Join by trait and rand_method
            ).map { trait, result_file, rand_method, random_files ->
                def random_dir = "${params.outdir}/magma_random/${rand_method}/${params.background}/${trait}"
                tuple(trait, "magma", result_file, random_dir)
            }
            
            // Calculate empirical p-values for MAGMA
            magma_empirical_results = calc_empirical_pvalues(magma_for_empirical)
            
            // Add to collection for combined results
            all_empirical_inputs = all_empirical_inputs.mix(magma_empirical_results)
            
            // Define channel for trait/method combinations from empirical results
            magma_by_trait_method = magma_empirical_results.map { trait, tool, emp_file ->
                tuple(trait, tool)
            }
        }
    }
    
    //////////////////////////////////////////
    // PRSET WORKFLOW
    //////////////////////////////////////////
    if (params.run_prset) {
        log.info "Running PRSet analysis for all traits"
        
        // Extract PRSet data without randomization (for deduplication)
        prset_dedup_data = trait_data
            .filter { trait_tuple ->
                // Filter out SCZ from PRSet analysis
                return trait_tuple[0].toUpperCase() != "SCZ" 
            }
            .map { trait_tuple ->
                def (trait, gwas_file, rsid_col, chr_col, pos_col, pval_col, n_col, 
                     binary_target, effect_allele, other_allele, summary_statistic_name, summary_statistic_type) = trait_tuple
                tuple(trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, 
                     summary_statistic_name, summary_statistic_type, "common")  // Add a placeholder for rand_method
            }
        
        // Run SNP deduplication ONCE per trait
        deduplicated_gwas = gwas_remove_dup_snps(prset_dedup_data)
        
        // Run real PRSet analysis for each randomization method
        prset_real_inputs = deduplicated_gwas
            .combine(Channel.fromList(params.randomization_methods))
            .map { trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, 
                  summary_statistic_name, summary_statistic_type, _, rand_method ->
                tuple(trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, 
                      summary_statistic_name, summary_statistic_type, rand_method)
            }
    
        // Run PRSet for real pathway sets
        real_prset_results = run_real_prset(prset_real_inputs)
        
        // Now combine with randomization methods for the random analysis
        prset_rand_inputs = deduplicated_gwas
            .combine(Channel.fromList(params.randomization_methods))
            .map { trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, 
                  summary_statistic_name, summary_statistic_type, _, rand_method ->
                tuple(trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, 
                      summary_statistic_name, summary_statistic_type, rand_method)
            }
        
        // Create random set inputs (one per permutation)
        perms_ch = Channel.from(1..params.num_random_sets)
        prset_random_inputs = prset_rand_inputs
            .combine(perms_ch)
            .map { trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, 
                  summary_statistic_name, summary_statistic_type, rand_method, perm ->
                tuple(trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, 
                      summary_statistic_name, summary_statistic_type, rand_method, perm)
            }
        
        // Run PRSet for random sets
        random_prset_results = run_random_sets_prset(prset_random_inputs)
        
        // Group random results by trait and randomization method
        random_prset_grouped = random_prset_results
            .groupTuple(by: [0, 2])  // Group by trait and rand_method

        // Calculate empirical p-values with explicit dependency
        if (params.run_empirical) {
            log.info "Setting up PRSet empirical p-value calculation"
            
            // Combine real results with grouped random results for empirical p-value calculation
            prset_for_empirical = real_prset_results.combine(
                random_prset_grouped,
                by: [0, 2]  // Join by trait and rand_method
            ).map { trait, result_files, rand_method, random_files ->
                def random_dir = "${params.outdir}/prset_random/${rand_method}/${params.background}/${trait}"
                tuple(trait, "prset", result_files, random_dir)
            }
            
            // Calculate empirical p-values for PRSet
            prset_empirical_results = calc_empirical_pvalues_prset(prset_for_empirical)
            
            // Add to collection for combined results
            all_empirical_inputs = all_empirical_inputs.mix(prset_empirical_results)
        }
    }
    
    //////////////////////////////////////////
    // OPENTARGETS COMPARISON WORKFLOW
    //////////////////////////////////////////
    if (params.run_empirical && params.run_opentargets) {
        log.info "Setting up OpenTargets comparison with multiple ranking methods"
        
        // Run OpenTargets comparison for MAGMA (if enabled)
        if (params.run_magma) {
            // Filter traits supported by OpenTargets
            magma_for_opentargets = magma_by_trait_method
                .groupTuple(by: [0, 1])
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
                .filter { trait, tool, birewire, keeppathsize -> 
                    params.opentargets_supported_traits.contains(trait)
                }
            
            // Step 3: Run correlation statistics analysis (new)
            if (params.run_ot_correlation) {
                magma_opentargets_correlation = opentargets_stats_correlation(magma_for_opentargets)
            }
        }
        
        // Run OpenTargets comparison for PRSet (if enabled)
        if (params.run_prset) {
            // Filter traits supported by OpenTargets  
            prset_for_opentargets = prset_by_trait_method
                .groupTuple(by: [0, 1])
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
                .filter { trait, tool, birewire, keeppathsize -> 
                    params.opentargets_supported_traits.contains(trait)
                }
            
            // Step 3: Run correlation statistics analysis (new)
            if (params.run_ot_correlation) {
                prset_opentargets_correlation = opentargets_stats_correlation(prset_for_opentargets)
            }
        }
    }
    
    //////////////////////////////////////////
    // TISSUE CORRELATION WORKFLOW - Independent of Tissue Specificity
    //////////////////////////////////////////
    if (params.run_empirical && params.run_tissue_correlation) {
        log.info "Setting up Tissue Correlation analysis (independent of Tissue Specificity)"
        
        // Run tissue correlation analysis for MAGMA
        if (params.run_magma) {
            log.info "Setting up Tissue Correlation analysis for MAGMA"
            
            // Create a channel for tissue correlation, identical setup to tissue specificity
            magma_for_tissue_corr = magma_by_trait_method
                .groupTuple(by: [0, 1])
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
            
            // Run tissue correlation analysis for all MAGMA traits
            tissue_correlation_analysis(magma_for_tissue_corr)
        }
        
        // Run tissue correlation analysis for PRSet
        if (params.run_prset) {
            log.info "Setting up Tissue Correlation analysis for PRSet"
            
            // Create a channel for tissue correlation, identical setup to tissue specificity
            prset_for_tissue_corr = prset_by_trait_method
                .groupTuple(by: [0, 1])
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
            
            // Run tissue correlation analysis for PRSet
            tissue_correlation_analysis(prset_for_tissue_corr)
        }
    }
    
    //////////////////////////////////////////
    // MALACARDS CORRELATION WORKFLOW
    //////////////////////////////////////////
    if (params.run_empirical && params.run_malacards_correlation) {
        log.info "Setting up MalaCards correlation analysis"

        // MAGMA
        if (params.run_magma) {
            log.info "MalaCards correlation for MAGMA"

            magma_for_malacards_corr = magma_by_trait_method
                .groupTuple(by: [0, 1]) // [trait, base_tool]
                .map { trait, base_tool, rand_methods, result_files ->
                    def birewire_idx = rand_methods.findIndexOf { it == 'birewire' }
                    def keeppathsize_idx = rand_methods.findIndexOf { it == 'keeppathsize' }
                    if (birewire_idx != -1 && keeppathsize_idx != -1) {
                        [trait, base_tool, result_files[birewire_idx], result_files[keeppathsize_idx]]
                    } else {
                        log.warn "Missing randomization method for ${trait} with ${base_tool}"
                        return null
                    }
                }
                .filter { it != null }

            malacards_correlation(magma_for_malacards_corr)
        }

        // PRSet
        if (params.run_prset) {
            log.info "MalaCards correlation for PRSet"

            prset_for_malacards_corr = prset_by_trait_method
                .groupTuple(by: [0, 1]) // [trait, base_tool]
                .map { trait, base_tool, rand_methods, result_files ->
                    def birewire_idx = rand_methods.findIndexOf { it == 'birewire' }
                    def keeppathsize_idx = rand_methods.findIndexOf { it == 'keeppathsize' }
                    if (birewire_idx != -1 && keeppathsize_idx != -1) {
                        [trait, base_tool, result_files[birewire_idx], result_files[keeppathsize_idx]]
                    } else {
                        log.warn "Missing randomization method for ${trait} with ${base_tool}"
                        return null
                    }
                }
                .filter { it != null }

            malacards_correlation(prset_for_malacards_corr)
        }
    }
    //////////////////////////////////////////
    // DELTA-RANK CORRELATION WORKFLOW
    //////////////////////////////////////////
    if (params.run_empirical && (params.run_delta_rank_ot || params.run_delta_rank_malacards)) {
        log.info "Setting up delta-rank correlation analyses"

        // Helper channel (BireWire-only result per trait/tool)
        def birewire_only = { ch ->
            ch.groupTuple(by: [0, 1])
              .map { trait, base_tool, rand_methods, result_files ->
                  def idx = rand_methods.findIndexOf { it == 'birewire' }
                  if (idx != -1) {
                      [trait, base_tool, result_files[idx]]
                  } else {
                      log.warn "Missing BireWire results for ${trait} with ${base_tool}"
                      return null
                  }
              }
              .filter { it != null }
        }

        // MAGMA
        if (params.run_magma) {
            def magma_bw = birewire_only(magma_by_trait_method)

            // Restrict OT delta-rank to whitelist only
            def magma_bw_ot = magma_bw.filter { trait, tool_base, birewire_file ->
                params.opentargets_supported_traits.contains(trait)
            }

            if (params.run_delta_rank_ot) {
                log.info "Delta-rank OT correlation for MAGMA (whitelist only)"
                delta_rank_ot_correlation(magma_bw_ot)
            }
            if (params.run_delta_rank_malacards) {
                log.info "Delta-rank MalaCards correlation for MAGMA"
                delta_rank_malacards_correlation(magma_bw)  // unfiltered
            }
        }

        // PRSet
        if (params.run_prset) {
            def prset_bw = birewire_only(prset_by_trait_method)

            // Restrict OT delta-rank to whitelist only
            def prset_bw_ot = prset_bw.filter { trait, tool_base, birewire_file ->
                params.opentargets_supported_traits.contains(trait)
            }

            if (params.run_delta_rank_ot) {
                log.info "Delta-rank OT correlation for PRSet (whitelist only)"
                delta_rank_ot_correlation(prset_bw_ot)
            }
            if (params.run_delta_rank_malacards) {
                log.info "Delta-rank MalaCards correlation for PRSet"
                delta_rank_malacards_correlation(prset_bw)  // unfiltered
            }
        }
    }
    
    //////////////////////////////////////////
    // GSA-MiXeR WORKFLOW
    //////////////////////////////////////////
    if (params.run_gsamixer) {
      log.info "Running GSA-MiXeR for all traits"

      // Build required input tuple (must include neff_col because module expects it)
      gsamixer_inputs = trait_data.map { trait, gwas_file, rsid_col, chr_col, pos_col, pval_col,
                                         n_col, binary_target, effect_allele, other_allele,
                                         summary_statistic_name, summary_statistic_type, se_col, neff_col ->
        tuple(
          trait,
          file(gwas_file),
          rsid_col,
          chr_col,
          pos_col,
          pval_col,
          n_col,
          binary_target,
          effect_allele,
          other_allele,
          summary_statistic_name,
          summary_statistic_type,
          se_col,
          neff_col
        )
      }

      // One-time reference generation
      ch_refs_result = convert_gmt_for_gsamixer(
        file(params.geneset_real),
        file(params.gtf_reference)
      )

      // Create individual channels from the outputs using a more explicit approach
      ch_baseline = Channel.fromPath("${params.outdir}/gsamixer_reference/baseline.txt")
      ch_full_gene = Channel.fromPath("${params.outdir}/gsamixer_reference/full_gene.txt") 
      ch_full_gene_set = Channel.fromPath("${params.outdir}/gsamixer_reference/full_gene_set.txt")

      // Prepare per-trait sumstats and split by chromosome (your existing logic)
      gsamixer_prepared = prepare_gsamixer_sumstats(gsamixer_inputs)
      gsamixer_split    = split_gsamixer_sumstats(gsamixer_prepared)

      // Wire baseline into every base job
      ch_base_in = gsamixer_split
        .combine(ch_baseline)                                    // adds the singleton baseline file to each tuple
        .map { trait, chrom_sumstats, baseline -> tuple(trait, chrom_sumstats, baseline) }

      gsamixer_base = gsamixer_plsa_base(ch_base_in)

      // Wire full_gene + full_gene_set into every full job
      ch_full_in = gsamixer_base
        .combine(ch_full_gene)                               // adds full_gene.txt
        .combine(ch_full_gene_set)                           // adds full_gene_set.txt
        .map { trait, base_json, base_log, base_weights, base_snps, full_gene, full_gene_set ->
            tuple(trait, base_json, base_log, base_weights, full_gene, full_gene_set, base_snps)
        }

      gsamixer_full = gsamixer_plsa_full(ch_full_in)
    }
    
    if (params.run_gsamixer && params.run_empirical) {
        log.info "Running GSA-MiXeR for random pathways - generating random gene sets once for all traits"
        
        // Create a channel of random GMT files - trait-independent
        random_gmt_files = Channel.fromList(
            params.randomization_methods.collect { rand_method ->
                (1..params.num_random_sets).collect { perm ->
                    def gmt_path = "${params.gmt_dirs[rand_method]}/GeneSet.random${perm}.gmt"
                    def gtf_path = params.gtf_reference
                    tuple(rand_method, perm, file(gmt_path), file(gtf_path))
                }
            }.flatten().collate(4)
        )
        
        // Convert random GMTs to GSA-MiXeR format - ONCE, independent of traits
        random_gmt_converted = convert_random_gmt_for_gsamixer(random_gmt_files)
        
        // For each trait's base model results, combine with ALL random set files
        gsamixer_trait_random_inputs = gsamixer_base
          .combine(random_gmt_converted)
          .map { trait, base_json, base_log, base_weights, base_snps, rand_method, perm, full_gene_txt, full_gene_set_txt ->
              // Make sure all elements are properly ordered in the tuple
              tuple(trait, 
                    file("${params.outdir}/gsamixer/${trait}/${trait}.chr*.sumstats.gz"), 
                    rand_method, 
                    perm, 
                    full_gene_txt, 
                    full_gene_set_txt,
                    base_json, 
                    base_log, 
                    base_weights,
                    base_snps)
          }
        
        // Run GSA-MiXeR full model for each trait-random set combination
        random_gmt_full_results = gsamixer_plsa_full_random(gsamixer_trait_random_inputs)
        
        // Group random results by trait and randomization method
        random_gmt_full_grouped = random_gmt_full_results
            .map { trait, rand_method, perm, json_file, log_file ->
                tuple(trait, rand_method, json_file)
            }
            .groupTuple(by: [0, 1])  // Group by trait and rand_method
        
        // Calculate empirical p-values
        if (params.run_empirical) {
            // Combine real results with grouped random results
            gsamixer_for_empirical = gsamixer_full
                .map { trait, full_json, full_log -> 
                    tuple(trait, full_json)
                }
                .combine(
                    random_gmt_full_grouped, 
                    by: 0  // Join by trait
                ).map { trait, real_json, rand_method, random_jsons ->
                    def random_dir = "${params.outdir}/gsamixer_random/${rand_method}/${trait}"
                    tuple(trait, "gsamixer", real_json, random_dir)
                }
            
            // Calculate empirical p-values for GSA-MiXeR
            gsamixer_empirical_results = calc_empirical_pvalues_gsamixer(gsamixer_for_empirical)
            
            // Add to collection for combined results
            all_empirical_inputs = all_empirical_inputs.mix(gsamixer_empirical_results)
        }
    }
}
