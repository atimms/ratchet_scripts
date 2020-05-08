#!/usr/bin/env python
import numpy as np
import sys
import os

##parameters
delim = '\t'

##working dir
working_dir = '/data/atimms/cherry_single_cell_0717/hu148'
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

def write_counts_sample_files(count_files, final_count_file, final_sample_sheet):
	sample_count = 0
	fc = 0
	with open(final_count_file, 'w') as fc_file, open(final_sample_sheet, 'w') as fss_file:
		for count_file in count_files:
			fc +=1
			with open(count_file, 'r') as c_file:
				lc = 0
				for line in c_file:
					lc += 1
					if lc == 1:
						header = line
						if fc == 1:
							fc_file.write(header)
							fss_file.write(delim.join(['sample', 'sample', 'donor']) + '\n')
					else:
						line = line.split(delim)
						sample = count_file.split('.')[0] + '_' + line[0]
						if 'Hu1' in count_file:
							donor = 'Hu1'
						elif 'Hu4' in count_file:
							donor = 'Hu4'
						elif 'Hu8' in count_file:
							donor = 'Hu8'
						else:
							print 'donor name not apparent in file:', count_file
						# print sample
						sample_count += 1
						print sample, sample_count
						fc_file.write(delim.join([sample] + line[1:]))
						fss_file.write(delim.join([sample, sample, donor]) + '\n')

def add_seurat_clusters_info(infile, seurat_data, outfile):
	##make dict with seurat data
	seurat_dict = {}
	with open(seurat_data, 'r') as seu_fh:
		lc = 0
		for line in seu_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(',')
				cell = line[0]
				seu_cluster = line[1]
				if cell not in seurat_dict:
					seurat_dict[cell] = seu_cluster
				else:
					print 'cell %s already seen'%cell
	# for c in seurat_dict:
	# 	print c, seurat_dict[c]
	with open(infile, 'r') as in_fh, open(outfile, 'w') as out_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.rstrip().split(',')
			line = [i.replace('"', '') for i in line]
			if lc == 1:
				out_fh.write(delim.join(line + ['seurat_cluster']) + '\n')
			else:
				cell = line[0]
				print cell
				if cell in seurat_dict:
					s_cluster = seurat_dict[cell]
				else:
					s_cluster = 'na'
				out_fh.write(delim.join(line + [s_cluster]) + '\n')

def compare_cell_charateristics(infile, outfile_prefix):
	cluster_list, s_cluster_list, cell_type_list = [],[],[]
	seurat_cell_type_dict, cluster_cell_type_dict, cluster_scluster_dict = {}, {}, {}
	##get lists information for 3 data file
	with open(infile, 'r') as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc >1:
				line = line.rstrip().split(delim)
				# print line
				cluster = line[6]
				scluster = line[11]
				cell_type = line[5]
				##make list of all types
				if cluster not in cluster_list:
					cluster_list.append(cluster)
				if scluster not in s_cluster_list:
					s_cluster_list.append(scluster)
				if cell_type not in cell_type_list:
					cell_type_list.append(cell_type)
				##populate comparison dicts
				if scluster in seurat_cell_type_dict:
					seurat_cell_type_dict[scluster].append(cell_type)
				else:
					seurat_cell_type_dict[scluster] = [cell_type]
				if cluster in cluster_cell_type_dict:
					cluster_cell_type_dict[cluster].append(cell_type)
				else:
					cluster_cell_type_dict[cluster] = [cell_type]
				if cluster in cluster_scluster_dict:
					cluster_scluster_dict[cluster].append(scluster)
				else:
					cluster_scluster_dict[cluster] = [scluster]

	print cluster_list, s_cluster_list, cell_type_list
	# '''
	##compare seurat clusters with cell types
	outfile = outfile_prefix + '.seurat_cluster_vs_cell_type.xls'
	with open(outfile, 'w') as out_fh:
		sc_count = 0
		for sc in seurat_cell_type_dict:
			sc_count += 1
			##make dict to get counts
			cell_type_count_dict = dict.fromkeys(cell_type_list, 0)
			print sc
			##get counts
			cell_types = seurat_cell_type_dict[sc]
			for cell_type in cell_types:
				cell_type_count_dict[cell_type] += 1

			##check and then print
			for i in cell_type_count_dict:
				print i, cell_type_count_dict[i]
			if sc_count ==1:
				out_fh.write(delim.join(['seurat_cluster'] + cell_type_count_dict.keys()) + '\n')
			counts = [str(i) for i in cell_type_count_dict.values()]
			out_fh.write(delim.join([sc] + counts) + '\n')
	##compare monocle clusters with cell types
	outfile = outfile_prefix + '.monocle_cluster_vs_cell_type.xls'
	with open(outfile, 'w') as out_fh:
		mc_count = 0
		for mc in cluster_cell_type_dict:
			mc_count += 1
			##make dict to get counts
			cell_type_count_dict = dict.fromkeys(cell_type_list, 0)
			print mc
			##get counts
			cell_types = cluster_cell_type_dict[mc]
			for cell_type in cell_types:
				cell_type_count_dict[cell_type] += 1

			##check and then print
			for i in cell_type_count_dict:
				print i, cell_type_count_dict[i]
			if mc_count ==1:
				out_fh.write(delim.join(['monocle_cluster'] + cell_type_count_dict.keys()) + '\n')
			counts = [str(i) for i in cell_type_count_dict.values()]
			out_fh.write(delim.join([mc] + counts) + '\n')
	# '''
	##compare monocle clusters with seurat clusters
	outfile = outfile_prefix + '.monocle_cluster_vs_seurat_clusters.xls'
	with open(outfile, 'w') as out_fh:
		mc2_count = 0
		for mc2 in cluster_scluster_dict:
			mc2_count += 1
			##make dict to get counts
			scluster_count_dict = dict.fromkeys(s_cluster_list, 0)
			print mc2
			##get counts
			s_clusters = cluster_scluster_dict[mc2]
			for s_cluster in s_clusters:
				scluster_count_dict[s_cluster] += 1

			##check and then print
			for i in scluster_count_dict:
				print i, scluster_count_dict[i]
			if mc2_count ==1:
				out_fh.write(delim.join(['monocle_cluster'] + scluster_count_dict.keys()) + '\n')
			counts = [str(i) for i in scluster_count_dict.values()]
			out_fh.write(delim.join([mc2] + counts) + '\n')
##run methods

##look at cells from seurat data
##make three input files
project_name = 'cherry_hu148_0817'
temp_matrix = project_name + '.temp_martix.txt'
raw_martix = project_name + '.raw_martix.txt'
sample_sheet = project_name + '.sample_sheet.txt'
gene_annotation = project_name + '.gene_ann.txt'
count_files = ['Hu1_retina_1.counts.tsv', 'Hu1_retina_2.counts.tsv', 'Hu4_retina_1.counts.tsv', 
		'Hu4_retina_2.counts.tsv', 'Hu4_retina_3.counts.tsv', 'Hu4_retina_4.counts.tsv', 
		'Hu8_retina_1.counts.tsv', 'Hu8_retina_2.counts.tsv', 'TJC_Hu1_ret_1.counts.tsv', 
		'TJC_Hu1_ret_2.counts.tsv']

##add gene counts to sample dict and write our counts and sample sheet
# write_counts_sample_files(count_files, temp_matrix, sample_sheet)
##transpose count file
# transpose_file(temp_matrix, raw_martix)
##get gene ids and make gene annotation file from transposed counts
# get_gene_ids(raw_martix, gene_annotation)


##run monocle in r and get cell information, so compare the cluser with makers and seurat clusters
test_name = 'cell_info1'
seurat_clusters = 'seuratClusters.csv'
cell_type_info = 'cherry_hu148_0817.' + test_name + '.csv'
cell_type_with_seurat_info = 'cherry_hu148_0817.' + test_name + '_with_seurat.xls'
comparision_prefix = 'cherry_hu148_0817.' + test_name
##add the seurat cluster info
add_seurat_clusters_info(cell_type_info, seurat_clusters, cell_type_with_seurat_info)
compare_cell_charateristics(cell_type_with_seurat_info, comparision_prefix)
##test2
test_name = 'cell_info2'
seurat_clusters = 'seuratClusters.csv'
cell_type_info = 'cherry_hu148_0817.' + test_name + '.csv'
cell_type_with_seurat_info = 'cherry_hu148_0817.' + test_name + '_with_seurat.xls'
comparision_prefix = 'cherry_hu148_0817.' + test_name
##add the seurat cluster info
add_seurat_clusters_info(cell_type_info, seurat_clusters, cell_type_with_seurat_info)
compare_cell_charateristics(cell_type_with_seurat_info, comparision_prefix)






