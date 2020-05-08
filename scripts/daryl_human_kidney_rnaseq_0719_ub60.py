#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import shutil


##parameters
delim = '\t'
threads = '16'

##setup working directory where results will be
working_dir = '/data/atimms/daryl_human_kidney_rnaseq_0719'
os.chdir(working_dir)


##methods
def compare_cluster_genes_to_de(c_files, d_files, out_file, p_req):
	print(p_req)
	de_dict = {}
	##make dict with which gene are differentially expressed
	for d_file in d_files:
		de_gene_count = 0
		with open(d_file, "r") as d_ph:
			lc = 0
			de_test = d_file.split('.')[1]
			for line in d_ph:
				lc += 1
				if lc >1:
					line = line.rstrip().split(',')
					gene = line[0].strip('"')
					adj_p = line[6]
					if adj_p == 'NA':
						# print(adj_p)
						adj_p = 1
					else:
						adj_p = float(adj_p)
					# print(adj_p)
					if adj_p <= p_req:
						de_gene_count += 1
						# print(line, gene, adj_p)
						if gene in de_dict:
							de_dict[gene].append(de_test)
						else:
							de_dict[gene] = [de_test]
		print(de_test, de_gene_count)
	print(len(de_dict))
	##get all DE info from de_dict for genes in a cluster
	cluster_dict = {}
	for c_file in c_files:
		with open(c_file, "r") as c_ph:
			lc = 0
			c_name = c_file.split('.')[0]
			for line in c_ph:
				lc += 1
				if lc >1:
					genename = line.rstrip()
					if genename in de_dict:
						# print(genename, de_dict[genename])
						de_results = de_dict[genename]
						if c_name in cluster_dict:
							cluster_dict[c_name].extend(de_results)
						else:
							cluster_dict[c_name] = de_results
	##get results
	with open(out_file, "w") as out_ph:
		d_names = [d.split('.')[1] for d in d_files]
		# print(d_names)
		header = ['cluster'] + d_names
		out_ph.write(delim.join(header) + '\n')
		for c in cluster_dict:
			# print(c, cluster_dict[c])
			c_de_list = cluster_dict[c]
			results = []		
			for d_name in d_names:
				de_count = c_de_list.count(d_name)
				# print(c, d_name, de_count)
				results.append(de_count)
			line_out = [c] + [str(r) for r in results]
			out_ph.write(delim.join(line_out) + '\n')






##run methods
cluster_genelists = ['day2_1.txt', 'day2_2.txt', 'day2_3.txt', 'day2_4.txt', 'day2_5.txt', 'day2_6.txt', 'day5_1.txt', 
	'day5_2.txt', 'day5_3.txt', 'day5_4.txt', 'day5_5.txt', 'day5_6.txt']
de_files = ['GSE126805_0719.3months_vs_1year.csv', 'GSE126805_0719.post_vs_1year.csv', 'GSE126805_0719.post_vs_3months.csv', 
	'GSE126805_0719.pre_vs_1year.csv', 'GSE126805_0719.pre_vs_3months.csv', 'GSE126805_0719.pre_vs_post.csv']
results_file_05 = 'GSE126805_0719.clusters_de.p05.xls'
results_file_01 = 'GSE126805_0719.clusters_de.p01.xls'
results_file_005 = 'GSE126805_0719.clusters_de.p005.xls'
# de_files = ['GSE126805_0719.3months_vs_1year.csv']

compare_cluster_genes_to_de(cluster_genelists, de_files, results_file_05, 0.05)
compare_cluster_genes_to_de(cluster_genelists, de_files, results_file_01, 0.01)
compare_cluster_genes_to_de(cluster_genelists, de_files, results_file_005, 0.005)
