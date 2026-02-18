process standardize_gwas_format {
    executor 'local'
    tag "${trait}_standardize"
    
    publishDir "${params.outdir}/standardized_gwas", mode: 'copy', overwrite: true
    
    input:
    tuple val(trait),
          path(gwas_file),
          val(rsid_col),
          val(chr_col),
          val(pos_col),
          val(pval_col),
          val(n_col),
          val(binary_target),
          val(effect_allele),
          val(other_allele),
          val(summary_statistic_name),
          val(summary_statistic_type),
          val(se_col),
          val(neff)
    
    output:
    tuple val(trait),
          path("${trait}_standardized.txt")
    
    script:
    """
    awk -v trait="${trait}" \\
        -v rsid="${rsid_col}" \\
        -v chr="${chr_col}" \\
        -v pos="${pos_col}" \\
        -v pval="${pval_col}" \\
        -v n="${n_col}" \\
        -v a1="${effect_allele}" \\
        -v a2="${other_allele}" \\
        -v beta="${summary_statistic_name}" \\
        -v se="${se_col}" \\
        -v neff="${neff}" \\
        'BEGIN {OFS="\\t"} 
         NR==1 {
            # Find column indices
            for (i=1; i<=NF; i++) {
                if (\$i == rsid) rsid_idx = i
                if (\$i == chr) chr_idx = i  
                if (\$i == pos) pos_idx = i
                if (\$i == pval) pval_idx = i
                if (\$i == n) n_idx = i
                if (\$i == a1) a1_idx = i
                if (\$i == a2) a2_idx = i
                if (\$i == beta) beta_idx = i
                if (\$i == se) se_idx = i
                if (\$i == neff) neff_idx = i
            }
            # Print standardized header
            print "RSID", "CHR", "POS", "PVAL", "N", "A1", "A2", "BETA", "SE", "NEFF"
         } 
         NR>1 {
            # Extract values using found indices, with fallbacks
            rsid_val = (rsid_idx ? \$(rsid_idx) : ".")
            chr_val = (chr_idx ? \$(chr_idx) : ".")
            pos_val = (pos_idx ? \$(pos_idx) : ".")
            pval_val = (pval_idx ? \$(pval_idx) : ".")
            n_val = (n_idx ? \$(n_idx) : ".")
            a1_val = (a1_idx ? \$(a1_idx) : ".")
            a2_val = (a2_idx ? \$(a2_idx) : ".")
            beta_val = (beta_idx ? \$(beta_idx) : ".")
            se_val = (se_idx ? \$(se_idx) : ".")
            neff_val = (neff_idx ? \$(neff_idx) : n_val)  # Use N as fallback for NEFF
            
            print rsid_val, chr_val, pos_val, pval_val, n_val, a1_val, a2_val, beta_val, se_val, neff_val
         }' ${gwas_file} > ${trait}_standardized.txt
    """
}