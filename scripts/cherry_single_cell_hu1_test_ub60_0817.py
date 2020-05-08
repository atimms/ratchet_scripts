#!/usr/bin/env python
import numpy as np
import sys
import os

##parameters
delim = '\t'

##working dir
working_dir = '/data/atimms/cherry_single_cell_0717/hu1_test'
os.chdir(working_dir)

def transpose_file(infile, outfile):
	with open(infile, 'r') as file, open(outfile, 'w') as final:
		#make list for each line
		a = [x.split() for x in file]
		#print a
		#turn into array
		b = np.array(a)
		#print b
		#transpose array
		c = b.T
		#convert array back to list of list
		d = c.tolist()
		#print d
		#convert each list back to 
		for e in d:
			line = delim.join(e) + '\n'
			final.write(line)

def get_cell_ids(infile, outfile):
	with open(infile, 'r') as in_file, open(outfile, 'w') as final:
		for line in in_file:
			line = line.split(delim)
			cell_id = line[0]
			final.write(delim.join([cell_id, cell_id + '\n']))

def get_gene_ids(infile, outfile):
	with open(infile, 'r') as in_file, open(outfile, 'w') as final:
		lc = 0
		for line in in_file:
			lc += 1
			if lc == 1:
				final.write('gene_short_name' + '\n')
			else:
				line = line.split(delim)
				gene_id = line[0]
				final.write(delim.join([gene_id, gene_id + '\n']))

##just looking at hu1_1 data
# '''
##make three input files
project_name = 'Hu1_retina_1'
original_file = 'Hu1_retina_1.counts.tsv'
raw_martix = project_name + '.raw_martix.txt'
sample_sheet = project_name + '.sample_sheet.txt'
gene_annotation = project_name + '.gene_ann.txt'
##expression matrix
transpose_file(original_file, raw_martix)
##sample cd
get_cell_ids(original_file, sample_sheet)
##gene annotation
get_gene_ids(raw_martix, gene_annotation)
# '''


