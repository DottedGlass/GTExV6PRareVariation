"""Filter for rare variants based on allele frequency <0.01 in both GTEx and 1KG
"""
import os
import argparse
import numpy as np
import pandas as pd

def read_frq(frq_file):
    """Read in frequency file.
    Returns
    - dataframe of original frequency file
    - dictionary with key=position, value=row index in dataframe
    """
    # read frequency file as dataframe
    df = pd.read_csv(frq_file, sep='\t', skiprows=1, header=None)

    # create list of positions
    pos = df[0] + ":" + df[1].astype(str)
    pos = list(pos.values)

    # create dictionary
    pos_idx = list(range(len(pos)))
    pos_dict = dict(zip(pos,pos_idx))

    return df, pos_dict

def filter_rare(gtex_df, gtex_pos_dict, oneKG_df, oneKG_dict, MAF_thres = 0.01):
    """Filter for rare variants based on allele frequncy < 0.01 in both GTEx and 1KG
    Returns gtex_rare_df
    """

    filter_idx = []
    # loop through each variant in GTEx
    for i in range(len(gtex_df)):
        if i % 100000 == 99999:
            print("Checked " + str(i+1) + " positions")

        # get allele frequency as dictionary (removes NaN entries)
        allele_list = list(gtex_df.iloc[i,:][4:]) # list causes NaN to be cast as float nan
        frq = {x.split(':')[0]:float(x.split(':')[1]) for x in allele_list if isinstance(x,str)}

        # find rare alleles (nonzero and less than MAF_thres)
        rare_alleles = [x for x in frq.keys() if frq[x] > 0 and frq[x] < MAF_thres]

        # continue only if there are rare alleles in GTEx
        if len(rare_alleles) > 0:

            # check if GTEx variant position exists in 1KG
            gtex_pos = gtex_df.iloc[i,0] + ":" + gtex_df.iloc[i,1].astype(str)
            if gtex_pos in oneKG_dict:

                # get row idx for oneKG_df
                j = oneKG_dict[gtex_pos]

                # get allele freqency from 1KG
                allele_list_KG = list(oneKG_df.iloc[j,:][4:]) # list causes NaN to be cast as float nan
                frq_KG = {x.split(':')[0]:float(x.split(':')[1]) for x in allele_list_KG if isinstance(x,str)}

                # check if alleles in GTEx are present in 1KG
                if all(x in frq_KG for x in rare_alleles):
                    # check if allele in GTEx has >= MAF_thres in 1KG
                    # this step filters out a variant if all the GTEx rare alleles are
                    # not rare in 1KG
                    common_alleles = [x for x in rare_alleles if frq_KG[x] >= MAF_thres]

                    if len(rare_alleles) > len(common_alleles):
                        # shared position, allele in GTEx is rare (0 < AF < MAF_thres), and
                        # allele in 1KG is rare (AF < MAF_thres)
                        filter_idx.append(i)
                else:
                    # shared position, allele in GTEx but not in 1KG
                    filter_idx.append(i)

            else:
                # position in GTEx but not in 1KG
                filter_idx.append(i)

    # filter for rare variants in GTEx
    gtex_rare_df = gtex_df.iloc[filter_idx,:]
    print("Keeping " + str(len(filter_idx) + " rare variants"))
    return(gtex_rare_df)

def save_filtered_frq(gtex_file,gtex_rare_df,outfile):
    """Saves gtex_rare_df as tab delimited file with same headers as gtex_file
    """

    # get header from original gtex_file
    with open(gtex_file,'r') as f:
        # header = f.readline().rstrip().split('\t')
        header = f.readline()

    # write filted gtex file with header only
    with open(outfile,'w') as f:
        f.write(header)

    # write the data to filtered gtex file
    # gtex_rare_df.to_csv(outfile, sep='\t', index=False, header=header)
    gtex_rare_df.to_csv(outfile, sep='\t', index=False, header=None, mode='a')


if __name__ == '__main__':

    # parse arguments from command line
    parser = argparse.ArgumentParser(description='Filter for rare variants based on allele frequency <0.01 in both GTEx and 1KG')
    parser.add_argument('GTEx', type=str, help='Input dir for GTEx frequency file')
    parser.add_argument('oneKG', type=str, help='Input dir for 1KG frequency file')
    parser.add_argument('outfile', type=str, help='output file for filtered list of rare variants')
    args = parser.parse_args()

    gtex_file = args.GTEx
    oneKG_file = args.oneKG
    outfile = args.outfile

    # parameters
    MAF_thres = 0.01 #threshold to be considered rare variant

    # read in gtex and 1KG frequency files
    gtex_df, gtex_pos_dict = read_frq(gtex_file)
    oneKG_df, oneKG_dict = read_frq(oneKG_file)

    # filter for rare variants
    gtex_rare_df = filter_rare(gtex_df, gtex_pos_dict, oneKG_df, oneKG_dict, MAF_thres=MAF_thres)

    # save filtered file
    save_filtered_frq(gtex_file,gtex_rare_df,outfile)
