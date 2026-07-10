#!/usr/bin/env python3
"""
PascalX Pathway Enrichment Wrapper
Performs pathway enrichment using pre-computed gene scores.
Output: CSV file with columns: Pathway_Name, Gene_List, Values, Pvalue
"""

import sys
import csv
import numpy as np
import os
import gc
sys.path.insert(0, '/opt/PascalX/python/')
import PascalX
from PascalX import genescorer
from PascalX import pathway

def main():
    if len(sys.argv) != 8:
        print("Usage: run_pascalx_pathways.py <trait> <gene_scores_file> <gwas_file> <gmt_file> <genome_annot> <ref_panel> <output_suffix>")
        sys.exit(1)

    trait = sys.argv[1]
    gene_scores_file = sys.argv[2]
    gwas_file = sys.argv[3]        # REQUIRED for meta-gene re-scoring (see load_GWAS below)
    gmt_file = sys.argv[4]
    genome_file = sys.argv[5]
    ref_panel = sys.argv[6]        # e.g. "/pascalx_ref/EUR.1KG.GRCh37"
    output_suffix = sys.argv[7]
    
    try:
        print(f"Starting pathway enrichment for trait: {trait}")
        print(f"Output suffix: {output_suffix}")
        print(f"GMT file: {gmt_file}")
        print(f"Genome annotation: {genome_file}")
        
        # Initialize gene scorer (needed for pathway analysis)
        Scorer = genescorer.chi2sum()
        print("Gene scorer initialized")

        # Load reference panel (REQUIRED: pathway scoring re-scores fused/meta-genes on the fly)
        Scorer.load_refpanel(ref_panel, parallel=1)
        print(f"Reference panel loaded: {ref_panel}")

        # Load GWAS summary stats (REQUIRED). Pathway scoring fuses overlapping genes
        # into meta-genes and re-scores them from SNP-level z-scores; load_scores only
        # restores single-gene scores, so without the GWAS every meta-gene fails with
        # "0 genes scored / can not be scored (check annotation)". Columns match the
        # prepare_pascalx_gwas output and the gene-scoring step.
        if not os.path.exists(gwas_file):
            raise FileNotFoundError(f"GWAS file not found: {gwas_file}")
        Scorer.load_GWAS(gwas_file, rscol=0, a1col=1, a2col=2, pcol=3, bcol=4, header=True)
        print(f"GWAS loaded: {gwas_file}")

        # Load genome annotation (required before loading scores)
        if not os.path.exists(genome_file):
            raise FileNotFoundError(f"Genome annotation not found: {genome_file}")
        Scorer.load_genome(genome_file, ccol=1, cid=0, csymb=0, cstx=2, cetx=3, cs=4, header=False)
        print("Genome annotation loaded")
        
        # Load pre-computed gene scores
        if not os.path.exists(gene_scores_file):
            raise FileNotFoundError(f"Gene scores file not found: {gene_scores_file}")
        Scorer.load_scores(gene_scores_file)
        print(f"Gene scores loaded from: {gene_scores_file}")
        
        # Initialize pathway scorer
        Pscorer = pathway.chi2rank(Scorer)
        print("Pathway scorer initialized")
        
        # Load gene set (GMT file)
        if not os.path.exists(gmt_file):
            raise FileNotFoundError(f"GMT file not found: {gmt_file}")
        # GMT format: name, description, gene1, gene2, ...
        # ncol=0 (pathway name in first column), fcol=2 (genes start from 3rd column, 0-indexed)
        M = Pscorer.load_modules(gmt_file, ncol=0, fcol=2)
        n_pathways = len(M)
        print(f"Loaded {n_pathways} pathways from GMT file")
        
        # Force garbage collection before scoring
        gc.collect()
        
        # Score pathways
        print("Scoring pathways (this may take several minutes)...")
        RESULT = Pscorer.score(M)
        print("Pathway scoring completed")
        
        # Save results to CSV
        output_file = f"{trait}_{output_suffix}.csv"
        with open(output_file, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['Pathway_Name', 'Gene_List', 'Values', 'Pvalue'])
            
            pathway_count = 0
            for item in RESULT[0]:
                if item:
                    pathway_name = item[0]
                    gene_list = ', '.join(item[1])
                    values = ', '.join(['nan' if np.isnan(val) else str(val) for val in item[2]])
                    pvalue = item[3]
                    writer.writerow([pathway_name, gene_list, values, pvalue])
                    pathway_count += 1
        
        print(f'Results saved to: {output_file}')
        print(f'Successfully processed {pathway_count} pathways')
        print("SUCCESS: Pathway enrichment complete")
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
