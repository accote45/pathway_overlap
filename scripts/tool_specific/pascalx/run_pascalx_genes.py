#!/usr/bin/env python3
"""
PascalX Gene Scoring Wrapper
Scores all genes using GWAS summary statistics and reference panel.
Output: gene_scores.txt file that can be reused for real/random pathway enrichment.
"""

import sys
import os
sys.path.insert(0, '/opt/PascalX/python/')
import PascalX
from PascalX import genescorer

def main():
    if len(sys.argv) != 5:
        print("Usage: run_pascalx_genes.py <trait> <gwas_file> <ref_panel> <genome_annot>")
        sys.exit(1)
    
    trait = sys.argv[1]
    gwas_file = sys.argv[2]
    ref_panel = sys.argv[3]  # Now expects mounted path like "/data/pascalx_reference/EUR.1KG.GRCh37"
    genome_file = sys.argv[4]  # Now expects "/data/pascalx_reference/msigdbgenes.regions"
    
    # Verify files exist at mounted paths
    print(f"Checking reference files at container paths:")
    print(f"  Reference panel prefix: {ref_panel}")
    print(f"  Genome annotation: {genome_file}")
    
    required_extensions = ['.bed', '.bim', '.fam']
    missing_files = []
    
    for ext in required_extensions:
        full_path = f"{ref_panel}{ext}"
        if not os.path.exists(full_path):
            missing_files.append(full_path)
        else:
            print(f"  Found: {full_path}")
    
    if missing_files:
        raise FileNotFoundError(
            f"Reference panel files not found at mounted paths. Missing: {', '.join(missing_files)}\n"
            f"Expected prefix: {ref_panel}\n"
            f"Container working directory: {os.getcwd()}\n"
            f"Available files: {os.listdir(os.path.dirname(ref_panel)) if os.path.exists(os.path.dirname(ref_panel)) else 'directory not found'}"
        )
    
    print(f"Reference panel validated: {ref_panel}")
    
    try:
        print(f"Starting gene scoring for trait: {trait}")
        
        # Initialize gene scorer with window size and MAF filter
        Scorer = genescorer.chi2sum(window=35000, MAF=0.01)
        print("Gene scorer initialized (window=35kb, MAF>=0.01)")
        
        # Load reference panel with parallelization
        Scorer.load_refpanel(ref_panel, parallel=15)
        print(f"Reference panel loaded: {ref_panel}")
        
        # Load GWAS file (already at correct mounted path from Nextflow)
        if not os.path.exists(gwas_file):
            raise FileNotFoundError(f"GWAS file not found: {gwas_file}")
        Scorer.load_GWAS(gwas_file, rscol=0, a1col=1, a2col=2, pcol=3, bcol=4, header=True)
        print(f"GWAS loaded: {gwas_file}")
        
        # Load genome annotation
        if not os.path.exists(genome_file):
            raise FileNotFoundError(f"Genome annotation not found: {genome_file}")
        Scorer.load_genome(genome_file, ccol=1, cid=0, csymb=0, cstx=2, cetx=3, cs=4, header=False)
        print("Genome annotation loaded")
        
        # Score all genes
        print("Scoring genes (this may take several minutes)...")
        R = Scorer.score_all(parallel=15)
        n_genes = len(R) if R else 0
        print(f"Gene scoring completed: {n_genes} genes scored")
        
        # Save gene scores for reuse
        output_file = 'gene_scores.txt'
        Scorer.save_scores(output_file)
        print(f"Gene scores saved to: {output_file}")
        print("SUCCESS: Gene scoring complete")
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
