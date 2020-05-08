#!/usr/bin/env python
import numpy as np
import sys
import os
import math
import glob

##parameters
delim = '\t'

##working dir
working_dir = '/data/atimms/cherry_single_cell_0717/seurat_de_1217'
os.chdir(working_dir)


def combine_scs_in_sample_sheet(in_file, out_file, cluster_info_dict):
	with open(in_file, 'r') as in_fh, open(out_file, 'w') as out_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc == 1:
				out_fh.write(line)
			else:
				line = line.rstrip().split(delim)
				sample = line[0]
				cluster = line[3]
				if cluster in cluster_info_dict:
					out_fh.write(delim.join(line[:3] + [cluster_info_dict[cluster]]) + '\n')
					print cluster, cluster_info_dict[cluster]
				else:
					out_fh.write(delim.join(line) + '\n')


def make_dict_from_expression_data(infile):
	exp_dict = {}
	with open(infile, 'r') as in_file:
		lc = 0
		for line in in_file:
			lc += 1
			line = line.rstrip().split(',')
			if lc == 1:
				cluster_indeces = [i.replace('"', '') for i in line[1:]]
				print cluster_indeces
			else:
				gene = line[0].replace('"', '')
				exp_values = [i.replace('"', '') for i in line[1:]]
				exp_values = [float(i) for i in exp_values]
				if gene in exp_dict:
					print gene, 'seen multiple times'
				else:
					exp_dict[gene] = exp_values
					print gene, exp_values, len(exp_values)
	return cluster_indeces, exp_dict

def combine_de_results_exp_values_two_clusters(infiles, sclus_index, exp_dict, file_prefix):
	for infile in infiles:
		clusters = infile.split('.')[1].split('_')
		outfile = file_prefix + '.' + '_'.join(clusters) + '.s_cluster.de.xls'
		cluster_index = [sclus_index.index(i) for i in clusters]
		print clusters, cluster_index
		with open(infile, 'r') as in_fh, open(outfile, 'w') as out_fh:
			lc = 0
			for line in in_fh:
				lc += 1
				line = line.rstrip().split(',')
				line = [i.replace('"', '') for i in line]
				if lc == 1:
					header = line + clusters + ['log2_fc']
					out_fh.write(delim.join(header) + '\n')
				else:
					gene = line[0]
					line_out = line
					exp_values = []
					for cluster in clusters:
						cluster_index = sclus_index.index(cluster)
						exp_value = exp_dict[gene][cluster_index]
						# print gene, cluster, cluster_index, exp_value
						exp_values.append(exp_value)
					if exp_values[0] == 0.0 or exp_values[1] == 0.0:
						fold_change = 'na'
					else:
						fold_change = exp_values[1]/ exp_values[0]
						# print fold_change
						fold_change = math.log(fold_change, 2)
					exp_values.append(fold_change)
					exp_values = [str(i) for i in exp_values]
					out_fh.write(delim.join(line_out + exp_values) + '\n')


def get_genes_increased_in_a_cluster(clusters):
	for cluster in clusters:
		infiles = glob.glob('retina_seurat_1217.' + cluster + '_*xls')
		outfile = 'retina_seurat_1217.' + cluster + '.sig_against_all_others.txt'
		# print cluster, len(infiles)

		sig_dict = {}
		for infile in infiles:
			other_cluster = infile.split('.')[1].split('_')[1]
			with open(infile, 'r') as in_fh:
				lc = 0
				for line in in_fh:
					lc += 1
					if lc >1:
						line = line.rstrip().split(delim)
						if line[4] != 'NA':
							qvalue = float(line[4])
							log2_fc = line[8]
							if log2_fc == 'na' or float(log2_fc) <0:							
								gene = line[0]
								info = [other_cluster] + line[6:]
								if qvalue <= 0.05:
									if gene in sig_dict:
										sig_dict[gene][0] += 1
										sig_dict[gene][1].extend(info)
									else:
										sig_dict[gene] = [1, info]
		with open(outfile, 'w') as out_fh:
			for g in sig_dict:
				if sig_dict[g][0] >=8:
					print cluster, g
					out_fh.write(delim.join([cluster, g] + sig_dict[g][1] + ['\n']))






##run methods

##parameters 
project_name = 'retina_seurat_1217'
cluster_dict = {'sc0':'sc034', 'sc3':'sc034', 'sc4':'sc034', 'sc2':'sc2911', 'sc9':'sc2911', 'sc11':'sc2911'}
original_sample_sheet = 'retina_seurat_0817.sample_sheet.txt'
combined_sample_sheet = 'retina_seurat_1217.sample_sheet.txt'
all_clusters = ['sc034', 'sc1', 'sc10', 'sc12', 'sc13', 'sc2911', 'sc5', 'sc6', 'sc7', 'sc8']

##before r: change cluster names i.e. combine 0,3,4 and 2, 9, 11
# combine_scs_in_sample_sheet(original_sample_sheet, combined_sample_sheet, cluster_dict)


##after r:
##combine de results with expression matrix
gene_expression_by_cluster = 'retina_seurat_1217.mean_expression_per_cluster.csv'
##make dict from expression values
# s_cluster_index, gebc_dict = make_dict_from_expression_data(gene_expression_by_cluster)
##1vs1 combination
# expression_files = glob.glob('retina_seurat_1217*de.csv' )
##add expression values to csv files
# combine_de_results_exp_values_two_clusters(expression_files, s_cluster_index, gebc_dict, project_name)
##then compare results: so find clusters which differ from all other clusters
get_genes_increased_in_a_cluster(all_clusters)

