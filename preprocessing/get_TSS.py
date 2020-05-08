"""Subsets *patched_contigs_TSS.bed file (created by process.reference.files.sh)
by the genes of interest (lincRNA and protein coding)
"""

import os
import operator

dir = os.environ["RAREVARDIR"]
dir = dir + '/reference/'
output = dir + 'v8.genes.TSS.bed'

genes_file = dir + 'v8.genes.lincRNA.protein.txt'
tss_file = dir + 'gencode.v26.genes.v8.patched_contigs_TSS.bed'

with open (genes_file) as f:
    genes = [line.rstrip('\n') for line in f]

with open (tss_file) as f:
    tss_all = [line.rstrip('\n').split('\t') for line in f]

# find matching lines
genes_set = set(genes)
tss_all_genes = [line[3] for line in tss_all]
idx = [i for i, item in enumerate(tss_all_genes) if item in genes_set]

# get TSS for target genes and write to file
tss = [tss_all[x] for x in idx]

with open(output, 'w') as f:
    for line in tss:
        f.write('\t'.join(line))
        f.write('\n')
