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
working_dir = '/data/atimms/vishal_rnaseq_0618'
working_dir = '/data/atimms/lan_pig_cpb_rnaseq_0719'
os.chdir(working_dir)

##ensemble db files, generated from table browser on 
##one ens id to anouther
ens_id_file = 'pig_ens_id.txt'
ens_id_gene_file = 'pig_ens_to_gene.txt'

def add_gene_id_to_de_results(infile, outfile):
	##make 2 dicts to query
	ens_id_dict, ens_gene_dict = {}, {}
	with open(ens_id_file, "r") as ens_ph:
		lc = 0
		for line in ens_ph:
			lc += 1
			if lc >1:
				line = line.split(delim)
				id1 = line[12].split('.')[0]
				id2 = line[1].split('.')[0]
				if id1 in ens_id_dict:
					ens_id_dict[id1].append(id2)
				else:
					ens_id_dict[id1] = [id2]
	with open(ens_id_gene_file, "r") as eg_ph:
		lc = 0
		for line in eg_ph:
			lc += 1
			if lc >1:
				line = line.rstrip().split(delim)
				ens_id = line[0].split('.')[0]
				gene_id = line[1]
				if ens_id in ens_gene_dict:
					ens_gene_dict[ens_id].append(gene_id)
				else:
					ens_gene_dict[ens_id] = [gene_id]
	##query file 
	with open(infile, "r") as in_ph, open(outfile, "w") as out_ph:
		lc = 0
		for line in in_ph:
			lc += 1
			line = line.split(',')
			if lc ==1:
				header = [line[0],'gene name'] + line[1:]
				out_ph.write(','.join(header))
			else:
				de_ens_id = line[0].replace('"', '')
				
				if de_ens_id[:2] != 'NA':
					# print(de_ens_id)
					if de_ens_id in ens_id_dict:
						new_de_ens_ids = ens_id_dict[de_ens_id]
						gene_ids = []
						for new_de_ens_id in new_de_ens_ids:
							if new_de_ens_id in ens_gene_dict:
								gene_id = ens_gene_dict[new_de_ens_id]
								gene_ids.extend(gene_id)
					else:
						gene_ids = ['']
					gene_ids = list(set(gene_ids))
					gene_ids = '_'.join(gene_ids)
					# print(de_ens_id, gene_ids)
					line_out = [line[0], gene_ids] + line[1:]
					out_ph.write(','.join(line_out))		

##run methods
de_files = ['visal_0618_hippocampus.CPB_DHCA_vsCPB_DHCAG_CSF.csv', 'visal_0618_hippocampus.sham_vs_CPB_DHCA.csv', 
		'visal_0618_hippocampus.sham_vsCPB_DHCAG_CSF.csv', 'visal_0618_striatum.CPB_DHCA_vs_CPB_DHCA_G_CSF.csv', 
		'visal_0618_striatum_reduced.CPB_DHCA_vs_CPB_DHCA_G_CSF.csv', 'visal_0618_striatum_reduced.sham_vs_CPB_DHCA.csv', 
		'visal_0618_striatum_reduced.sham_vs_CPB_DHCA_G_CSF.csv', 'visal_0618_striatum.sham_vs_CPB_DHCA.csv', 
		'visal_0618_striatum.sham_vs_CPB_DHCA_G_CSF.csv']
# de_files = ['visal_0618_hippocampus.CPB_DHCA_vsCPB_DHCAG_CSF.csv']
de_files = ['lan_pig_cpb_0719.Hypoxia_vs_Control.csv']


for de_file in de_files:
	new_de_file = de_file.rsplit('.', 1)[0] + '.gene_id_added.csv'
	add_gene_id_to_de_results(de_file, new_de_file)







