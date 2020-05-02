import os

dir = os.environ["RAREVARDIR"]
dir = dir + '/reference/'

genes_file = dir + 'v8.genes.lincRNA.protein.txt'
tss_file = dir + 'gencode.v26.genes.v8.patched_contigs_TSS.bed'

with open (genes_file) as f:
    genes = [line.rstrip('\n') for line in f]

with open (tss_file) as f:
    tss_all = [line.rstrip('\n').split('\t') for line in f]

# find matching lines
idx = []

for i in range(len(genes)):
    gene = genes[i]
    for j in range(len(tss_all)):
        tss_gene = tss_all[j][3]
        if gene == tss_gene:
            idx.append(j)

tss = [tss_all[x] for x in idx]

with open(dir + 'v8.tss.txt', 'w') as f:
    for line in tss:
        for val in line:
            f.write(val + '\t')

        f.write('\n')
