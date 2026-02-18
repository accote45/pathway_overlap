#!/usr/bin/env python3

import argparse
import sys
import os
import pandas as pd

# Add the current directory to path so we can import pascalx.py
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

def main():
    parser = argparse.ArgumentParser(description='Run PascalX pathway enrichment with standardized GWAS format')
    parser.add_argument('--gwas-file', required=True, help='Path to standardized GWAS file')
    parser.add_argument('--trait', required=True, help='Trait name')
    parser.add_argument('--geneset-file', required=True, help='Path to gene set file (GMT)')
    parser.add_argument('--output', required=True, help='Output file path')
    
    args = parser.parse_args()
    
    # Since we're using standardized format, column indices are fixed:
    # RSID=0, CHR=1, POS=2, PVAL=3, N=4, A1=5, A2=6, BETA=7, SE=8, NEFF=9
    fixed_columns = {
        'rsid_col_idx': 0,
        'chr_col_idx': 1, 
        'pos_col_idx': 2,
        'pval_col_idx': 3,
        'n_col_idx': 4,
        'a1_col_idx': 5,
        'a2_col_idx': 6,
        'beta_col_idx': 7,
        'se_col_idx': 8,
        'neff_col_idx': 9
    }
    
    print(f"Processing trait: {args.trait}")
    print(f"GWAS file: {args.gwas_file}")
    print(f"Gene set file: {args.geneset_file}")
    print(f"Using standardized column indices: {fixed_columns}")
    
    try:
        # Import and call your existing PascalX functions
        # Note: You'll need to adapt this to match your actual pascalx.py interface
        from pascalx import run_pathway_enrichment  # Adjust this import as needed
        
        result = run_pathway_enrichment(
            gwas_file=args.gwas_file,
            geneset_file=args.geneset_file,
            trait=args.trait,
            **fixed_columns
        )
        
        # Write results to output file
        if isinstance(result, pd.DataFrame):
            result.to_csv(args.output, sep='\t', index=False)
        else:
            # Handle other result formats as needed
            with open(args.output, 'w') as f:
                f.write(str(result))
                
        print(f"Results written to: {args.output}")
        
    except Exception as e:
        print(f"Error running PascalX: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()