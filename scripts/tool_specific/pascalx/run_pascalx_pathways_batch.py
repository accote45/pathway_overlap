#!/usr/bin/env python3
"""
PascalX Pathway Enrichment Wrapper -- BATCHED (random sets only)

Loads the reference panel, genome annotation, and pre-computed gene scores ONCE,
then scores a contiguous range of randomized GMT files, writing one CSV per
permutation with EXACTLY the same filename the per-permutation script produced:

    {trait}_random{perm}.{rand_method}.csv

This is behaviour-preserving for everything downstream: only the number of
process launches (and redundant reference-panel loads) changes. The reference
panel is identical across all permutations of a trait, so loading it once per
batch instead of once per permutation is the whole point of this script.
"""

import sys
import csv
import os
import gc
import numpy as np
sys.path.insert(0, '/opt/PascalX/python/')
from PascalX import genescorer
from PascalX import pathway


def write_result(result, output_file):
    """Write one PascalX pathway RESULT to CSV (same format as the single script)."""
    pathway_count = 0
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Pathway_Name', 'Gene_List', 'Values', 'Pvalue'])
        for item in result[0]:
            if item:
                pathway_name = item[0]
                gene_list = ', '.join(item[1])
                values = ', '.join(['nan' if np.isnan(val) else str(val) for val in item[2]])
                pvalue = item[3]
                writer.writerow([pathway_name, gene_list, values, pvalue])
                pathway_count += 1
    return pathway_count


def main():
    if len(sys.argv) != 9:
        print("Usage: run_pascalx_pathways_batch.py <trait> <gene_scores_file> "
              "<gwas_file> <genome_annot> <ref_panel> <rand_method> <start_perm> <end_perm>")
        sys.exit(1)

    trait = sys.argv[1]
    gene_scores_file = sys.argv[2]
    gwas_file = sys.argv[3]        # REQUIRED for meta-gene re-scoring (see load_GWAS below)
    genome_file = sys.argv[4]
    ref_panel = sys.argv[5]
    rand_method = sys.argv[6]
    start_perm = int(sys.argv[7])
    end_perm = int(sys.argv[8])

    # Same container-internal layout the per-permutation process used:
    #   /randomized_gene_sets/random_<rand_method>/GeneSet.random<perm>.gmt
    gmt_dir = f"/randomized_gene_sets/random_{rand_method}"

    try:
        print(f"Batch pathway enrichment for trait: {trait}")
        print(f"  Randomization method: {rand_method}")
        print(f"  Permutation range:    {start_perm}..{end_perm}")
        print(f"  GMT directory:        {gmt_dir}")

        # ---- Heavy setup: performed ONCE for the whole batch ----
        Scorer = genescorer.chi2sum()
        print("Gene scorer initialized")

        # Reference panel is required: pathway scoring re-scores fused/meta-genes
        # on the fly. It is identical for every permutation, so load it just once.
        Scorer.load_refpanel(ref_panel, parallel=1)
        print(f"Reference panel loaded: {ref_panel}")

        # GWAS summary stats are REQUIRED for that fused/meta-gene re-scoring: it needs
        # SNP-level z-scores, which load_scores does NOT restore. Identical across every
        # permutation, so load once. Without this every meta-gene fails with
        # "0 genes scored / can not be scored (check annotation)".
        if not os.path.exists(gwas_file):
            raise FileNotFoundError(f"GWAS file not found: {gwas_file}")
        Scorer.load_GWAS(gwas_file, rscol=0, a1col=1, a2col=2, pcol=3, bcol=4, header=True)
        print(f"GWAS loaded: {gwas_file}")

        if not os.path.exists(genome_file):
            raise FileNotFoundError(f"Genome annotation not found: {genome_file}")
        Scorer.load_genome(genome_file, ccol=1, cid=0, csymb=0, cstx=2, cetx=3, cs=4, header=False)
        print("Genome annotation loaded")

        if not os.path.exists(gene_scores_file):
            raise FileNotFoundError(f"Gene scores file not found: {gene_scores_file}")
        Scorer.load_scores(gene_scores_file)
        print(f"Gene scores loaded from: {gene_scores_file}")

        Pscorer = pathway.chi2rank(Scorer)
        print("Pathway scorer initialized")

        # ---- Loop over permutations, reusing the loaded reference/scores ----
        n_done = 0
        for perm in range(start_perm, end_perm + 1):
            gmt_file = os.path.join(gmt_dir, f"GeneSet.random{perm}.gmt")
            if not os.path.exists(gmt_file):
                raise FileNotFoundError(f"GMT file not found: {gmt_file}")

            M = Pscorer.load_modules(gmt_file, ncol=0, fcol=2)
            print(f"[perm {perm}] Loaded {len(M)} pathways from {gmt_file}")

            gc.collect()
            RESULT = Pscorer.score(M, parallel=int(os.environ.get("PASCALX_PARALLEL", "1")))

            output_file = f"{trait}_random{perm}.{rand_method}.csv"
            n = write_result(RESULT, output_file)
            print(f"[perm {perm}] Wrote {n} pathways -> {output_file}")

            del M, RESULT
            gc.collect()
            n_done += 1

        print(f"SUCCESS: Batch complete ({n_done} permutations scored)")

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
