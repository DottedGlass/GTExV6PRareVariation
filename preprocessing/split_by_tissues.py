#!/usr/bin/env python

import os
import sys
import fileinput
import argparse
import gzip
import numpy as np
import scipy as sp

############## FUNCTIONS

def sample_dict_maker(SAMPLE_fh):
	'''
	Reads in sample annotation file and makes dict. Keys are tissue IDs, and values are sample IDs.
	'''
	SAMPLE = open(SAMPLE_fh)
	SAMPLE.readline()
	sample_dict = {}
	missing_data  = []
	for line in SAMPLE:
		split_line = line.rstrip().split('\t')
		if len(split_line) == 2:
			sample, tissue = split_line
			if tissue != '':
				sample_dict[sample] = tissue
		else:
			missing_data.append(line)
	SAMPLE.close()
	if missing_data != []:
		print("Missing data on these line(s)")
		print(missing_data)
	return sample_dict

def split_by_tissue(sample_dict, GTEX_fh, OUT_DIR, END):
	'''
	Read in GTEx RPKM file and split by tissues.
	'''
	print("Splitting by tissue...")
	GTEX = open(GTEX_fh)
	for line in GTEX:
		line = line.rstrip().split('\t')
		if 'Name' in line:
			header = line
			break

	tissue_dict = {}
	for sample in header[2:]:
		if sample in sample_dict:
			tissue = sample_dict[sample]
			if tissue not in tissue_dict:
				tissue_dict[tissue] = []
			tissue_dict[tissue].append(header.index(sample))

	header = np.array(header)
	file_dict = {}
	for tissue in tissue_dict:
		file_dict[tissue] = open(OUT_DIR + '/' + tissue + END, 'w')
		out_names = header[tissue_dict[tissue]]
		for i in range(len(out_names)):
			out_names[i] = '-'.join(out_names[i].split('-')[0:2])
		file_dict[tissue].write('\t'.join(['Gene'] + out_names.tolist()) + '\n')

	for line in GTEX:
		line = np.array(line.rstrip().split())
		for tissue in tissue_dict:
			file_dict[tissue].write('\t'.join([line[0]] + line[tissue_dict[tissue]].tolist()) + '\n')

	for tissue in file_dict:
		file_dict[tissue].close()
	
	print("Done")



############## MAIN
USAGE = """
	Splits combined RPKM file from GTEx into expression matrices by tissue.
    """

parser = argparse.ArgumentParser(description = USAGE)
parser.add_argument('--GTEX', dest = 'GTEX', required = True, help = 'GTEx combined RPKM file')
parser.add_argument('--OUT', dest = 'OUT', required = True, help = 'directory to hold output RPKM matrices for each tissue')
parser.add_argument('--SAMPLE', dest = 'SAMPLE', required = True, help = 'sample annotation file')
parser.add_argument('--END', dest = 'END', required = True, help = 'file ending')

args = parser.parse_args()

GTEX_fh = args.GTEX
SAMPLE_fh = args.SAMPLE
OUT_DIR = args.OUT
END = args.END

## Make sample dict and get list of tissues
sample_dict = sample_dict_maker(SAMPLE_fh)

## Split GTEx combined file by tissue
split_by_tissue(sample_dict, GTEX_fh, OUT_DIR, END)
