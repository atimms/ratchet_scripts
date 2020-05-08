#!/usr/bin/env python
import numpy as np
import sys
import os
import math

##parameters
delim = '\t'

##working dir
working_dir = '/data/atimms/cherry_single_cell_0717/mouse_human_retina_seurat_0719'
os.chdir(working_dir)


def get_shared_genes(file1, file2):
	gene_dict = {}
	##read and get genenames in dict
	for infile in [file1, file2]:
		with open(infile, 'r') as in_fh:
			lc = 0
			for line in in_fh:
				lc += 1
				if lc >1:
					line = line.rstrip().split(delim)
					gene = line[0]
					if gene in gene_dict:
						gene_dict[gene] += 1
					else:
						gene_dict[gene] = 1
	##read and make new files for genes in both files
	for infile in [file1, file2]:
		outfile = infile.rsplit('.', 1)[0] + '.shared_genes.txt'
		with open(infile, 'r') as in_fh, open(outfile, 'w') as out_fh:
			lc = 0
			for line in in_fh:
				lc += 1
				if lc == 1:
					out_fh.write(line)
				else:
					line = line.split(delim)
					gene = line[0]
					if gene_dict[gene] == 2:
						out_fh.write(delim.join(line))


mouse_matrix = 'GSE63472_P14Retina_merged_digital_expression.txt'
human_matrix = 'cherry_hu148_0817.raw_martix.txt'

get_shared_genes(mouse_matrix, human_matrix)

