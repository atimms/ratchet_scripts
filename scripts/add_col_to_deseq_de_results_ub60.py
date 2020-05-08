#!/usr/bin/env python
import sys
import subprocess
import os


##parameters
delim = ','

def add_if_in_genelist_deseq_results(infile, outfile, genes):
	genelist = []
	##get list of genes
	with open(genes, "U") as gene_fh:
		for line in gene_fh:
			line = line.rstrip()
			genelist.append(line)
	print(genelist)
	##update file
	with open(infile, "r") as in_fh, open(outfile, "w") as out_fh:
		line_count = 0
		for line in in_fh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				out_fh.write(delim.join(line + ['genelist', '\n']))
			else:
				gene = line[0].strip('"')
				# print line[0], gene
				if gene in genelist:
					out_fh.write(delim.join(line + ['in genelist', '\n']))
					print('gene found:', gene)
				else:
					out_fh.write(delim.join(line + ['not in genelist', '\n']))



# add_if_in_genelist_deseq_results('kim_rnaseq_1016.de.egl_vs_pk.csv', 'kim_rnaseq_1016.de.egl_vs_pk.asd176.csv', 'ASD_GeneList176.txt')
# add_if_in_genelist_deseq_results('kim_rnaseq_1016.de.egl_vs_pk.csv', 'kim_rnaseq_1016.de.egl_vs_pk.chd.csv', 'chd_genelist_0517.txt')
# add_if_in_genelist_deseq_results('kim_rnaseq_1016.de.egl_vs_pk.csv', 'kim_rnaseq_1016.de.egl_vs_pk.epilesy.csv', 'epilesy_genelist_0517.txt')
# add_if_in_genelist_deseq_results('kim_rnaseq_0517.pcl_egl.de.csv', 'kim_rnaseq_0517.pcl_egl.de.asd162.csv', 'ASD_162gene_list.txt')
# add_if_in_genelist_deseq_results('kim_rnaseq_0717.pcl_egl.de.csv', 'kim_rnaseq_0717.pcl_egl.de.asd161.csv', 'ASD_161gene_list.txt')
# add_if_in_genelist_deseq_results('kim_rnaseq_0717.pcl_egl.de.csv', 'kim_rnaseq_0717.pcl_egl.de.scz156.csv', 'SCZ_156genelist.txt')

##first set of gene list
# gene_files = ['Arvantis_2d_8d.txt', 'Arvantis_2d_SO.txt', 'Arvantis_8d_2d.txt', 'Arvantis_8d_SO.txt', 'Arvantis_SO_2d.txt', 'Arvantis_SO_8d.txt']
##second set
# gene_files = ['Grgic_d0_d2.txt', 'Grgic_d0_d5.txt', 'Grgic_d2_d0.txt', 'Grgic_d5_d0.txt']
# de_files = ['acomy_rnaseq_0618_mouse.sham_vs_day2.csv', 'acomy_rnaseq_0618_mouse.sham_vs_day5.csv', 
# 	'acomy_rnaseq_0618_acomy.sham_vs_day2.genename.csv', 'acomy_rnaseq_0618_acomy.sham_vs_day5.genename.csv']
# for genelist in gene_files:
# 	for de_file in de_files:
# 		new_file = de_file.split('.')[0] + '.' + de_file.split('.')[1] + '.' + genelist.split('.')[0] + '.volcano.csv'
# 		add_if_in_genelist_deseq_results(de_file, new_file, genelist)

##for laura 0119
# gene_files = ['human_known_genes.txt', 'GSE1987.txt', 'GSE33532.txt', 'GSE44077.txt','GSE1987_down.txt', 'GSE33532_down.txt', 'GSE44077_down.txt']
# de_files = ['laura_0119.results.tumor_vs_normal.csv']
# for genelist in gene_files:
# 	for de_file in de_files:
# 		new_file = de_file.split('.')[0] + '.' + de_file.split('.')[1] + '.' + genelist.split('.')[0] + '.volcano.csv'
# 		add_if_in_genelist_deseq_results(de_file, new_file, genelist)

##tim 1019


add_if_in_genelist_deseq_results('cherry_rnaseq_1019.ret_rpe.retina_vs_rpe.csv', 'cherry_rnaseq_1019.ret_rpe.retina_vs_rpe.volcano.csv', 'genelist.txt')