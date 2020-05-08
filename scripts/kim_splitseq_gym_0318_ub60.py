#!/usr/bin/env python
import numpy as np
import sys

##parameters
delim = ','


def transpose_file(infile, outfile):
	with open(infile, 'r') as file, open(outfile, 'w') as final:
		#make list for each line
		a = [x.split(',') for x in file]
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


original_matrix = 'digital_gene_expression.csv'
transposed_matrix = 'splitseq_0318.tp_matrix.csv'

transpose_file(original_matrix, transposed_matrix)