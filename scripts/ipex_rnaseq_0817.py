#!/usr/bin/env python
import subprocess
import os
from operator import itemgetter

##parameters
delim = '\t'

##working dir
working_dir = '/data/atimms/ipex_rnaseq_0817'
os.chdir(working_dir)

##file names
all_counts = 'ipex_gene_edit_0716.htseq_counts.txt'
all_metadata = 'ipex_gene_edit_0716.metadata.modified_names.txt'
counts_suffix = '.htseq_counts.txt'
metadata_suffix = '.metadata.txt'
file_prefix = 'ipex_rnaseq_0817'



def get_subset_of_data(name, group_req):
	samples_to_keep = []
	metadata_out_file = file_prefix + '.' + name + metadata_suffix
	counts_out_file = file_prefix + '.' + name + counts_suffix
	##metadata
	with open(all_metadata, "r") as min_fh, open(metadata_out_file, "w") as mout_fh:
		lc = 0
		for line in min_fh:
			lc += 1
			if lc == 1:
				mout_fh.write(line)
			else:
				line = line.rstrip().split(delim)
				sample = line[0]
				group = line[1]
				if group in group_req:
					samples_to_keep.append(sample)
					mout_fh.write(delim.join(line) + '\n')
	# print samples_to_keep
	##count data
	with open(all_counts, "r") as cin_fh, open(counts_out_file, "w") as cout_fh:
		lc = 0
		for line in cin_fh:
			lc += 1
			line = line.rstrip().split(delim)
			if lc == 1:
				all_samples = line
				indices_to_keep = [i for i, x in enumerate(line) if x in samples_to_keep]
				##add genename
				indices_to_keep = [0] + indices_to_keep
				print indices_to_keep
				header = itemgetter(*indices_to_keep)(line)
				cout_fh.write(delim.join(header) + '\n')
			else:
				line_out = itemgetter(*indices_to_keep)(line)
				cout_fh.write(delim.join(line_out) + '\n')





##run methods
##experiment1
exp_name = 'exp1'
names_req = ['Activated_Teffector', 'Activated_Treg', 'Edited_Treg']
get_subset_of_data(exp_name, names_req)




