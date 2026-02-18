



###### run pathway enrichment analysis

import sys
import os
import csv
import numpy as np

rand = sys.argv[2]
path = sys.argv[1]

sys.path.append('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/software/PascalX/python/')
sys.path.append('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/software/PascalX/core')

# change working direcotry
os.chdir('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/')

# import genescorer class, initialize sum of chisq based gene scorer
from PascalX import genescorer
Scorer = genescorer.chi2sum(window=35000,MAF=0.01)
	# default window=50KB, varcutoff=0.99

# set reference panel
#Scorer.load_refpanel('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/results/genegene_correlation_test/1KG_filtered/Genome_phase3_ref_panel_EUR',keepfile=None)
Scorer.load_refpanel('./software/PascalX/misc/temp/bmi.EUR.1KG.GRCh37',parallel=10)

# ensembl
from PascalX.genome import genome

#import gene annotation file
Scorer.load_genome('./data/msigdb.regions',ccol=1,cid=5,csymb=0,cstx=2,cetx=3,cs=4,header=True)

# load GWAS
Scorer.load_GWAS('./data/gwas/SNP_gwas_mc_merge_nogc.tbl.uniq.gz',rscol=0,a1col=1,a2col=2,pcol=6,bcol=4,header=True)
#a1col=1,a2col=2
#Scorer.plot_genesnps('AOAH',mark_window=True,show_correlation=True);
	#to visualize GWAS SNPs gene-wise


# score all genes in annotation
#R = Scorer.score_all()
# score genes on chr21-22
#R = Scorer.score_chr(chrs=[21,22])

# to save gene symbol pvals in a tab separated text file
#Scorer.save_scores(f'./results/pathwaydb_enrichment/pascal_real/gene_scores_{ref_file_number}.txt')



###### pathway analysis

# import saved gene scores
Scorer.load_scores(f'./results/pathwaydb_enrichment/pascal_real/bmi/gene_scores_msigdb.txt')

from PascalX import pathway
Pscorer = pathway.chi2rank(Scorer) # rank based scoring (recommended)
#Pscorer = pathway.chi2perm(Scorer) # mote carlo based scoring

# load tab separated file of pathways
M = Pscorer.load_modules(f'/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/data/randomized_gene_sets/random_{path}/GeneSet.random{rand}.gmt',ncol=0,fcol=2)
###ncol= is the column with the name of the module and fcol= the first column with a gene symbol. It is assumed that other member genes follow in subsequent columns.

# scoring
RESULT = Pscorer.score(M)

# Define the file name for saving the results
output_file = f'./results/pathwaydb_enrichment_msigdbgenes/pascal_{path}/bmi/randset_{rand}.csv'

# Open the output file in write mode
with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    # Write the header row
	writer.writerow(['Pathway Name','Gene List','Values','Pvalue'])
    # Iterate over the results
	for item in RESULT[0]:
        if item:  # Make sure the item is not empty
            pathway_name = item[0]  # Pathway name
            gene_list = ', '.join(item[1])  # Join the gene list into a single string
            values = ', '.join(['nan' if np.isnan(val) else str(val) for val in item[2]])  # Convert array to string, handling NaN
            score = item[3]  # Score
            # Write each row
            writer.writerow([pathway_name, gene_list, values, score])

print(f'Results saved to {output_file}')






















