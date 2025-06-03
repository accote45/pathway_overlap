//////nextflow.nf

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

////////////////////////////////////////////////////////////////////
//                  Setup Channels
////////////////////////////////////////////////////////////////////

// ---------------------------------------------------
// ----  GWAS  ------ //
// ---------------------------------------------------

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
        content.binary_target,                    // [7] Binary trait flag
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
    trait_data

    main:
    // Extract only the fields needed for MAGMA
    magma_data = trait_data.map { fullData ->
        def (trait, gwas_file, rsid_col, chr_col, pos_col, pval_col, n_col, binary_target, effect_allele, other_allele, summary_statistic_name, summary_statistic_type) = fullData
        tuple(trait, gwas_file, rsid_col, chr_col, pos_col, pval_col, n_col)
    }
    
    prepared = prepare_input(magma_data)
    // Select gene file based on background setting
    def selected_gene_file = params.gene_files[params.background]

    // Append gene_file path to each tuple in snp_loc_data
    def snp_loc_with_gene_file = prepared.snp_loc_data.map { it + [selected_gene_file] }

    annotated = annotate_genes(snp_loc_with_gene_file)
    gene_analysis = run_gene_analysis(annotated.gene_annot_data)
    run_real_geneset(gene_analysis.gene_results)
    
    perms_ch = Channel.from(1..params.num_random_sets)
    random_sets_inputs = gene_analysis.gene_results
        .combine(perms_ch)
        .map { trait, gene_result, perm -> 
            println "Creating MAGMA job for ${trait} with permutation ${perm}"
            [trait, gene_result, perm] 
        }

    run_random_sets(random_sets_inputs)
    
    emit:
    gene_results = gene_analysis.gene_results
}

////////////////////////////////////////////////////////////////////
//                  Subworkflow: PRSET analysis
////////////////////////////////////////////////////////////////////

workflow prset {
    take:
    trait_data

    main:
    // Map the comprehensive data to the format needed for PRSet
    prset_data = trait_data.map { fullData ->
        def (trait, gwas_file, rsid_col, chr_col, pos_col, pval_col, n_col, binary_target, effect_allele, other_allele, summary_statistic_name, summary_statistic_type) = fullData
        tuple(trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, summary_statistic_name, summary_statistic_type)
    }
    
    // First, remove duplicate SNPs from GWAS files
    deduplicated_gwas = gwas_remove_dup_snps(prset_data)
    
    perms_ch = Channel.from(1..params.num_random_sets)
    
    // Combine prset data with permutation numbers
    prset_random_inputs = deduplicated_gwas
        .combine(perms_ch)
        .map { trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, summary_statistic_name, summary_statistic_type, perm ->
            println "Creating PRSET job for ${trait} with permutation ${perm}"
            [trait, gwas_file, binary_target, effect_allele, other_allele, rsid_col, pval_col, summary_statistic_name, summary_statistic_type, perm]
        }
    
    // Run PRSet for random sets
    run_random_sets_prset(prset_random_inputs)
}

////////////////////////////////////////////////////////////////////
//                  Main workflow
////////////////////////////////////////////////////////////////////

workflow {
    log.info "Pipeline configuration:"
    log.info "  Run MAGMA: ${params.run_magma}"
    log.info "  Run PRSet: ${params.run_prset}"
    
    if (params.run_magma) {
        log.info "Running MAGMA analysis"
        magma_results = magma(trait_data)

        if (params.run_prset) {
            // Force collection of all MAGMA results before PRSet
            all_magma_done = magma_results.gene_results.collect()
        }

    } else {
        log.info "Skipping MAGMA analysis"
    }

    if (params.run_prset) {
        log.info "Running PRSet analysis"
        prset(trait_data)
    } else {
        log.info "Skipping PRSet analysis"
    }
        // If neither is enabled, show warning
    if (!params.run_magma && !params.run_prset) {
        log.warn "Neither MAGMA nor PRSet is enabled. Set --run_magma true or --run_prset true"
    }
}








    prset(trait_data)
}


