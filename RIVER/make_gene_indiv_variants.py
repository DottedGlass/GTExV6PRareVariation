"""Makes file with following columns
Col1: Gene
Col2: Indiv
Col3: Comma separated list of variants

Col1 is restricted to genes that have at least one multi-tissue outlier
    individual
Col2 is all individuals that have at least one rare variant within 10kb TSS
    of the gene in that row
Col3 has each variant of the form $chrom:$position:$major_allele:$variant_allele
"""

import os
import argparse
import pandas as pd
import allel

print('hard coded files')
dir = os.environ["RAREVARDIR"] + '/' # upper-level directory
outliers_file = dir + 'data/v8/outliers_medz_picked.txt'
vcf_file = dir + 'download/gtex8/GTEx_AFA_10kb_TSS.vcf.gz'
TSS_file = dir + 'reference/v8.genes.TSS_minus10k.bed'
rare_var_file = dir + 'reference/GTEx_AFA_10kb_TSS_AF_rare.frq'

# if __name__ == '__main__':
#
#     # parse arguments from command line
#     parser = argparse.ArgumentParser(description='Make file with outlier gene and individual pairs with list of rare variants')
#     parser.add_argument('outliers', type=str, help='Input dir for multi-tissue outliers file')
#     parser.add_argument('vcf', type=str, help='Input dir for vcf file')
#     parser.add_argument('TSS', type=str, help='Input dir for BED file with 10kb TSS of genes')
#     parser.add_argument('rare_var', type=str, help='Input dir for allele frequency file with only rare variants')
#     parser.add_argument('outfile', type=str, help='output file')
#     args = parser.parse_args()
#
#     dir = os.environ["RAREVARDIR"] + '/' # upper-level directory
#
#     outliers_file = dir + args.outliers     # data/v8/outliers_medz_picked.txt
#     vcf_file = dir + args.vcf               # download/gtex8/GTEx_AFA_10kb_TSS.vcf
#     TSS_file = dir + args.TSS               # reference/v8.genes.TSS_minus10k.bed
#     rare_var_file = dir + args.rare_var     # reference/GTEx_AFA_10kb_TSS_AF_rare.frq
#     outfile = dir + args.outfile

# read in files
outliers = pd.read_csv(outliers_file, sep='\t')
outliers.columns = map(str.lower, outliers.columns)

vcf = allel.read_vcf(vcf_file, numbers={'ALT': 1}) # assumes biallelic

TSS = pd.read_csv(TSS_file, sep='\t', names=['chrom', 'start', 'stop', 'gene'])

rare_var = pd.read_csv(rare_var_file, sep='\t', skiprows=1, names=['chrom', 'pos', 'n_alleles', 'n_chr', 'ref', 'alt'])
rare_var_chrom_pos = rare_var['chrom'] + ':' + rare_var['pos'].astype(str)
rare_var.insert(rare_var.shape[1],'chrom:pos', rare_var_chrom_pos)
rare_var.drop_duplicates('chrom:pos',inplace=True) # remove duplicates

# filter TSS for genes with outliers
TSS = TSS[TSS['gene'].isin(outliers['gene'])]

# filter rare_var for variants that are within a 10kb TSS window of outlier gene
# note that TSS has 0-based indexing whereas rare_var has 1-based indexing
gene_dict = dict()      # dictionary with key=index, value=gene
rare_var['pos0'] = rare_var['pos'] - 1   #get 0-based indexing of rare_var
chr_list = list(TSS.drop_duplicates('chrom')['chrom'])
for chr in chr_list[:2]:
    TSS_chr = TSS[TSS['chrom'] == chr]
    rare_var_chr = rare_var[rare_var['chrom'] == chr]

    for _, row in TSS_chr.iterrows():
        # find rows in rare_var that are within a 10kb TSS of the outlier gene
        start = row['start']; stop = row['stop'] + 1
        idx_bool = rare_var_chr['pos0'].between(start,stop)
        rare_var_chr_idx = list(rare_var_chr[idx_bool].index)

        # keep track of the gene associated with these variants
        for idx in rare_var_chr_idx:
            if idx not in gene_dict.keys():
                gene_dict[idx] = row['gene']

rare_var_idx = list(gene_dict.keys())
rare_var_idx.sort()
rare_var = rare_var.iloc[rare_var_idx,:]

# add column for the gene associated with the position in rare_var
genes = [gene_dict[x] for x in rare_var_idx]
rare_var.insert(0, 'gene', genes)

# build dataframe from vcf file with columns
# chrom pos ref alt chrom:pos
vcf_df = pd.DataFrame({'chrom':vcf['variants/CHROM'], 'pos':vcf['variants/POS'], 'ref':vcf['variants/REF'], 'alt':vcf['variants/ALT']})
vcf_df_chrom_pos = vcf_df['chrom'] + ':' + vcf_df['pos'].astype(str)
vcf_df.insert(vcf_df.shape[1], 'chrom:pos', vcf_df_chrom_pos)
dup_index = vcf_df.drop_duplicates('chrom:pos').index # remove duplicates
vcf_df = vcf_df.iloc[dup_index,:]

# get genotype data from vcf file
genotype = vcf['calldata/GT']


# for each rare variant, find individuals that have the rare variant
col1 = []
col2 = []
col3 = []

# find matching chrom and pos in rare_var and vcf_df
rare_var = rare_var[rare_var['chrom:pos'].isin(vcf_df['chrom:pos'])]
vcf_df = vcf_df[vcf_df['chrom:pos'].isin(rare_var['chrom:pos'])]

print(rare_var)
print(vcf_df)

# $chrom:$position:$major_allele:$variant_allele
