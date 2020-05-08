#!/usr/bin/env python
import numpy as np
import sys
import os
import math

##parameters
delim = '\t'

##working dir
working_dir = '/data/atimms/cherry_single_cell_0717/seurat_de'
os.chdir(working_dir)

##methods
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

def make_sample_dict(infile):
	s_dict = {}
	with open(infile, 'r') as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc >1:
				line = line.rstrip().split(',')
				sample_id = line[0].strip('"')
				donor = sample_id.split('_')[-4]
				seurat_cluster = line[1].strip('"')
				# print sample_id, donor, seurat_cluster
				if sample_id in s_dict:
					print sample_id, 'is seen twice !!'
				else:
					s_dict[sample_id] = [donor, seurat_cluster]
	return s_dict

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

def write_counts_sample_files(sample_info_dict, count_files, final_count_file, final_sample_sheet):
	new_sample_dict = {}
	sample_count = 0
	for count_file in count_files:
		with open(count_file, 'r') as c_file:
			lc = 0
			for line in c_file:
				lc += 1
				if lc == 1:
					header = line
				else:
					line = line.split(delim)
					sample = count_file.split('.')[0] + '_' + line[0]
					# print sample
					if sample in sample_info_dict:
						sample_count += 1
						print sample, sample_count
						new_sample_dict[sample] = [sample_info_dict[sample], line[1:]]
	with open(final_count_file, 'w') as fc_file, open(final_sample_sheet, 'w') as fss_file:
		fc_file.write(header)
		fss_file.write(delim.join(['sample', 'sample', 'donor', 's_cluster']) + '\n')
		for s in new_sample_dict:
			fc_file.write(delim.join([s] + new_sample_dict[s][1]))
			fss_file.write(delim.join([s, s] + new_sample_dict[s][0]) + '\n')

def select_clusters(clusters_req, sample_sheet, temp_matrix, prefix):
	for cluster_req in clusters_req:
		final_sample_sheet = prefix + '.'+ '_'.join(cluster_req) + '.sample_sheet.txt'
		temp_matrix_file = prefix + '.'+ '_'.join(cluster_req) + '.matrix_temp.txt'
		final_matrix_file = prefix + '.'+ '_'.join(cluster_req) + '.matrix.txt'

		samples_to_keep = []
		##new sample sheet
		with open(sample_sheet, 'r') as samp_fh, open(final_sample_sheet, 'w') as fsamp_fh:
			lc = 0
			for line in samp_fh:
				lc += 1
				if lc == 1:
					fsamp_fh.write(line)
				else:
					line = line.rstrip().split(delim)
					sample = line[0]
					cluster = line[3][2:]
					if cluster in cluster_req:
						samples_to_keep.append(sample)
						fsamp_fh.write(delim.join(line) + '\n')
						print cluster, sample
		print len(samples_to_keep)
		##new matrix file
		with open(temp_matrix, 'r') as mat_fh, open(temp_matrix_file, 'w') as tmat_fh:
			lc = 0
			for line in mat_fh:
				lc += 1
				if lc == 1:
					tmat_fh.write(line)
				else:
					line = line.rstrip().split(delim)
					sample = line[0]
					if sample in samples_to_keep:
						tmat_fh.write(delim.join(line) + '\n')
		##then transpoe
		transpose_file(temp_matrix_file, final_matrix_file)

def select_clusters_vs_rest(clusters_req, sample_sheet, prefix):
	for cluster_req in clusters_req:
		final_sample_sheet = prefix + '.'+ '_'.join(cluster_req) + '.vs_rest.sample_sheet.txt'
		##new sample sheet
		with open(sample_sheet, 'r') as samp_fh, open(final_sample_sheet, 'w') as fsamp_fh:
			lc = 0
			for line in samp_fh:
				lc += 1
				line = line.rstrip().split(delim)
				if lc == 1:
					fsamp_fh.write(delim.join(line + ['in_cluster']) + '\n')
				else:
					cluster = line[3][2:]
					if cluster in cluster_req:
						fsamp_fh.write(delim.join(line + ['yes']) + '\n')
					else:
						fsamp_fh.write(delim.join(line + ['no']) + '\n')


def select_multiple_clusters(clusters_req, sample_sheet, temp_matrix, prefix):
	for cluster_req in clusters_req:
		final_sample_sheet = prefix + '.'+ '_'.join(cluster_req[0]) + '_vs_' + '_'.join(cluster_req[1]) + '.sample_sheet.txt'
		temp_matrix_file = prefix + '.'+ '_'.join(cluster_req[0]) + '_vs_' + '_'.join(cluster_req[1]) + '.matrix_temp.txt'
		final_matrix_file = prefix + '.'+ '_'.join(cluster_req[0]) + '_vs_' + '_'.join(cluster_req[1]) + '.matrix.txt'
		samples_to_keep = []
		##new sample sheet
		with open(sample_sheet, 'r') as samp_fh, open(final_sample_sheet, 'w') as fsamp_fh:
			lc = 0
			for line in samp_fh:
				lc += 1
				line = line.rstrip().split(delim)
				if lc == 1:
					fsamp_fh.write(delim.join(line + ['in_cluster']) + '\n')
				else:
					sample = line[0]
					cluster = line[3][2:]
					if cluster in cluster_req[0]:
						fsamp_fh.write(delim.join(line + ['yes']) + '\n')
						samples_to_keep.append(sample)
					if cluster in cluster_req[1]:
						fsamp_fh.write(delim.join(line + ['no']) + '\n')
						samples_to_keep.append(sample)
						print cluster, sample
		print len(samples_to_keep)
		##new matrix file
		with open(temp_matrix, 'r') as mat_fh, open(temp_matrix_file, 'w') as tmat_fh:
			lc = 0
			for line in mat_fh:
				lc += 1
				if lc == 1:
					tmat_fh.write(line)
				else:
					line = line.rstrip().split(delim)
					sample = line[0]
					if sample in samples_to_keep:
						tmat_fh.write(delim.join(line) + '\n')
		##then transpoe
		transpose_file(temp_matrix_file, final_matrix_file)

def make_dict_from_expression_data(infile):
	exp_dict = {}
	with open(infile, 'r') as in_file:
		lc = 0
		for line in in_file:
			lc += 1
			line = line.rstrip().split(',')
			if lc == 2:
				cluster_indeces = [i.replace('"', '') for i in line[1:]]
				print cluster_indeces
			elif lc >= 2:
				gene = line[0].replace('"', '')
				exp_values = [i.replace('"', '') for i in line[1:]]
				exp_values = [float(i) for i in exp_values]
				if gene in exp_dict:
					print gene, 'seen multiple times'
				else:
					exp_dict[gene] = exp_values
					# print gene, exp_values, len(exp_values)
	return cluster_indeces, exp_dict

def combine_de_results_exp_values_two_clusters(tests, sclus_index, exp_dict, file_prefix):
	for test in tests:
		infile = file_prefix + '.' + '_'.join(test) + '.s_cluster.de_with_hu.csv'
		outfile = file_prefix + '.' + '_'.join(test) + '.s_cluster.de_with_hu.xls'
		clusters = ['sc' + i for i in test]
		# cluster_index = [sclus_index.index(i) for i in clusters]
		# print clusters, cluster_index
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
						print gene, cluster, cluster_index, exp_value
						exp_values.append(exp_value)
					if exp_values[0] == 0.0 or exp_values[1] == 0.0:
						fold_change = 'na'
					else:
						fold_change = exp_values[1]/ exp_values[0]
						print fold_change
						fold_change = math.log(fold_change, 2)
					exp_values.append(fold_change)
					exp_values = [str(i) for i in exp_values]
					out_fh.write(delim.join(line_out + exp_values) + '\n')

def combine_de_results_exp_values_vs_rest(tests, sclus_index, exp_dict, file_prefix):
	for test in tests:
		infile = file_prefix + '.' + '_'.join(test) + '.vs_rest.s_cluster.de_with_hu.csv'
		outfile = file_prefix + '.' + '_'.join(test) + '.vs_rest.s_cluster.de_with_hu.xls'
		clusters = ['sc' + i for i in test]
		all_cluster_index = [sclus_index.index(i) for i in clusters]
		print clusters, all_cluster_index
		with open(infile, 'r') as in_fh, open(outfile, 'w') as out_fh:
			lc = 0
			for line in in_fh:
				lc += 1
				line = line.rstrip().split(',')
				line = [i.replace('"', '') for i in line]
				if lc == 1:
					header = line + ['_'.join(clusters), 'the_rest'] + ['log2_fc']
					out_fh.write(delim.join(header) + '\n')
				else:
					gene = line[0]
					line_out = line
					exp_values, expr_values_the_rest = [], []
					i_count = 0
					for i in exp_dict[gene]:
						if i_count in all_cluster_index:
							exp_values.append(i)
						else:
							expr_values_the_rest.append(i)
						i_count += 1
					# print exp_values, len(exp_values)
					# print expr_values_the_rest, len(expr_values_the_rest)
					cluster_mean = sum(exp_values) / len(exp_values)
					the_rest_mean = sum(expr_values_the_rest) / len(expr_values_the_rest)
					if cluster_mean == 0.0 or the_rest_mean == 0.0:
						fold_change = 'na'
					else:
						fold_change = the_rest_mean/ cluster_mean
						# print fold_change
						fold_change = math.log(fold_change, 2)
					values_out = [cluster_mean, the_rest_mean, fold_change]
					values_out = [str(i) for i in values_out]
					out_fh.write(delim.join(line_out + values_out) + '\n')

def combine_de_results_exp_values_vs_multiple_groups(tests, sclus_index, exp_dict, file_prefix):
	for test in tests:
		infile = file_prefix + '.' + '_'.join(test[0]) + '_vs_' + '_'.join(test[1]) + '.s_cluster.de_with_hu.csv'
		outfile = file_prefix + '.' + '_'.join(test[0]) + '_vs_' + '_'.join(test[1]) + '.s_cluster.de_with_hu.xls'
		clusters_1 = ['sc' + i for i in test[0]]
		cluster_index_1 = [sclus_index.index(i) for i in clusters_1]
		clusters_2 = ['sc' + i for i in test[1]]
		cluster_index_2 = [sclus_index.index(i) for i in clusters_2]
		print clusters_1, cluster_index_1, clusters_2, cluster_index_2
		with open(infile, 'r') as in_fh, open(outfile, 'w') as out_fh:
			lc = 0
			for line in in_fh:
				lc += 1
				line = line.rstrip().split(',')
				line = [i.replace('"', '') for i in line]
				if lc == 1:
					header = line + ['_'.join(clusters_1), '_'.join(clusters_2)] + ['log2_fc']
					out_fh.write(delim.join(header) + '\n')
				else:
					gene = line[0]
					line_out = line
					exp_values_1, exp_values_2 = [], []
					i_count = 0
					for i in exp_dict[gene]:
						if i_count in cluster_index_1:
							exp_values_1.append(i)
						if i_count in cluster_index_2:
							exp_values_2.append(i)
						i_count += 1
					# print exp_values, len(exp_values)
					# print expr_values_the_rest, len(expr_values_the_rest)
					cluster_mean_1 = sum(exp_values_1) / len(exp_values_1)
					cluster_mean_2 = sum(exp_values_2) / len(exp_values_2)
					if cluster_mean_1 == 0.0 or cluster_mean_2 == 0.0:
						fold_change = 'na'
					else:
						fold_change = cluster_mean_2/ cluster_mean_1
						# print fold_change
						fold_change = math.log(fold_change, 2)
					values_out = [cluster_mean_1, cluster_mean_2, fold_change]
					values_out = [str(i) for i in values_out]
					out_fh.write(delim.join(line_out + values_out) + '\n')




##run methods

##look at cells from seurat data
##make three input files
project_name = 'retina_seurat_0817'
temp_matrix = project_name + '.temp_martix.txt'
raw_martix = project_name + '.raw_martix.txt'
sample_sheet = project_name + '.sample_sheet.txt'
gene_annotation = project_name + '.gene_ann.txt'
seurat_clusters = 'seuratClusters.csv'
count_files = ['Hu1_retina_1.counts.tsv', 'Hu1_retina_2.counts.tsv', 'Hu4_retina_1.counts.tsv', 
		'Hu4_retina_2.counts.tsv', 'Hu4_retina_3.counts.tsv', 'Hu4_retina_4.counts.tsv', 
		'Hu8_retina_1.counts.tsv', 'Hu8_retina_2.counts.tsv', 'TJC_Hu1_ret_1.counts.tsv', 
		'TJC_Hu1_ret_2.counts.tsv']

##parameters for choosing clusters
##select clusters one vs one ie 0 vs 1
# clusters_req = [['0','6'], ['0','1'], ['0','3'], ['0','4']]
clusters_req = [['0','4'], ['3','4'], ['4','5']]
##select clusters vs the rest i.e. 0+3 vs. everything else
# clusters_vs_rest_req = [['0', '3'], ['1'], ['2', '9', '11'], 
# 		['4'], ['5'], ['6'], ['7'], ['8'], ['10'], ['12'], ['13']]
clusters_vs_rest_req = [['0', '3', '4']]
##multiples vs multiples
multiple_clusters = [[['0'], ['1', '2', '5', '6', '7', '8', '9', '10', '11', '12', '13']]]


'''
##prep files for monocle analysis

##make dict from seurat cluster file
# sample_dict = make_sample_dict(seurat_clusters)
##add gene counts to sample dict and write our counts and sample sheet
# write_counts_sample_files(sample_dict, count_files, temp_matrix, sample_sheet)
##transpose count file
# transpose_file(temp_matrix, raw_martix)
##get gene ids and make gene annotation file from transposed counts
# get_gene_ids(raw_martix, gene_annotation)
##test on one
##select clusters one vs one ie 0 vs 1
# select_clusters(clusters_req, sample_sheet, temp_matrix, project_name)
##select clusters vs the rest
# select_clusters_vs_rest(clusters_vs_rest_req, sample_sheet, project_name)
##two groups i.e. 0 and 1 vs 5, 6 and 7
select_multiple_clusters(multiple_clusters, sample_sheet, temp_matrix, project_name)
'''

# '''
##post analysis

##combine de results with expression matrix
gene_expression_by_cluster = 'retina_seurat_0817.mean_expression_per_cluster.csv'
##make dict from expression values
s_cluster_index, gebc_dict = make_dict_from_expression_data(gene_expression_by_cluster)
##1vs1
# combine_de_results_exp_values_two_clusters(clusters_req, s_cluster_index, gebc_dict, project_name)
##vs the rest
# combine_de_results_exp_values_vs_rest(clusters_vs_rest_req, s_cluster_index, gebc_dict, project_name)
##two groups i.e. 0 and 1 vs 5, 6 and 7
combine_de_results_exp_values_vs_multiple_groups(multiple_clusters, s_cluster_index, gebc_dict, project_name)
# '''

















