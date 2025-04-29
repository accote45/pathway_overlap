//////nextflow.nf

nextflow.enable.dsl=2

////////////////////////////////////////////////////////////////////
//                  Module inclusion
////////////////////////////////////////////////////////////////////

include {
    prepare_input
    annotate_genes
    run_gene_analysis
    run_real_geneset
    run_random_sets
} from './modules/magma/magma.nf'

////////////////////////////////////////////////////////////////////
//                  Setup Channels
////////////////////////////////////////////////////////////////////

// ---------------------------------------------------
// ----  GWAS  ------ //
// ---------------------------------------------------

import groovy.json.JsonSlurper

// Parse the JSON input into a Channel
def jsonSlurper = new JsonSlurper()
def configFile = new File(params.traits_config)
def phenoConfig = jsonSlurper.parseText(configFile.text)

sumstat = Channel.fromList(phenoConfig.collect { content ->
    tuple(
        content.trait,
        content.gwas_file,
        content.rsid_col,
        content.chr_col,
        content.pos_col,
        content.pval_col,
        content.n_col
    )
})

//Channel
//    .from(1..params.num_random_sets)
//    .set { rand_set_ids }

////////////////////////////////////////////////////////////////////
//                  Subworkflow: MAGMA analysis
////////////////////////////////////////////////////////////////////

workflow magma {
    take:
    sumstat

    main:
    prepared = prepare_input(sumstat)
    // Select gene file based on background setting
    def selected_gene_file = params.gene_files[params.background]

    // Append gene_file path to each tuple in snp_loc_data
    def snp_loc_with_gene_file = prepared.snp_loc_data.map { it + [selected_gene_file] }

    annotated = annotate_genes(snp_loc_with_gene_file)
    gene_analysis = run_gene_analysis(annotated.gene_annot_data)
    run_real_geneset(gene_analysis.gene_results)
    gene_analysis.gene_results

//    run_random_sets(gene_analysis.gene_results)

perms_ch = Channel.from(1..params.num_random_sets)
random_sets_inputs = gene_analysis.gene_results
    .combine(perms_ch)
    .map { trait, gene_result, perm -> 
        println "Creating job for ${trait} with permutation ${perm}"
        [trait, gene_result, perm] 
    }

run_random_sets(random_sets_inputs)
}

////////////////////////////////////////////////////////////////////
//                  Main workflow
////////////////////////////////////////////////////////////////////

workflow {
    magma(sumstat)
}
