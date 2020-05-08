#!/usr/bin/env python
import sys
import subprocess
import os

'''
module load local_python/3.6.5
'''
##parameters
delim = ','
# delim = '\t'
working_dir = '/data/atimms/milena_rnaseq_0818'
os.chdir(working_dir)

##methods
def split_dataframe_into_groups_get_average(in_file, outfile):
	mean_dict = {}
	with open(in_file, 'r') as in_fh:
		lc = 0
		sample_group_list = ['gene']
		sample_groups = []
		for line in in_fh:
			lc += 1
			line = line.rstrip().split(',')

			if lc == 1:
				samples = line[1:]
				# print(samples)
				for s in samples:
					# print(s)
					s = s.strip('"')
					if s.startswith('Y'):
						s = s.split('_')[0] + '_' + s.split('_')[-2]
					else:
						s = s.split('_')[0] + '_' + s.split('_')[2] + '_' + s.split('_')[-2]
					# print(s)
					sample_group_list.append(s)
					if s not in sample_groups:
						sample_groups.append(s)
				# print(sample_group_list, sample_groups)
			else:
				gene = line[0]
				# print(sample_groups)
				gene_means = []
				for sg in sample_groups:
					indices = [i for i, x in enumerate(sample_group_list) if x == sg]
					# print(line)
					values = [ float(line[i]) for i in indices]
					value_mean = sum(values) / len(values)
					# print(sg,indices, values, value_mean)
					gene_means.append(str(value_mean))
				if gene not in mean_dict:
					mean_dict[gene] = gene_means
				else:
					print(gene, ' gene seen twice')

	with open(outfile, 'w') as out_fh:
		out_fh.write(','.join(['gene'] + sample_groups) + '\n')
		for g in mean_dict:
			out_fh.write(','.join([g] + mean_dict[g]) + '\n')

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
				out_fh.write(delim.join(line + ['label', '\n']))
			else:
				gene = line[0].strip('"')
				log_fc = abs(float(line[2]))
				padj = line[6]
				if padj == 'NA':
					padj = 1.0
				else:
					padj = float(padj)
				print(line, gene, log_fc, padj)
				if gene in genelist:
					out_fh.write(delim.join(line + ['in gene set', '\n']))
					print('gene found:', gene)
				else:
					if log_fc >= 1 and padj <= 0.001:
						out_fh.write(delim.join(line + ['significant', '\n']))
						print(line, 'sig..........')
					else:
						out_fh.write(delim.join(line + ['not significant', '\n']))

def make_ind_genelists(infile):
	out_files = []
	with open(infile, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split('\t')
			outfile = line[0] + '.txt'
			genes = line[1:]
			with open(outfile, "w") as out_fh:
				out_files.append(outfile)
				for g in genes:
					out_fh.write(g + '\n')

	return(out_files)



##get data frames from rnaseq analysis
vst_counts = 'milena_rnaseq_0818_all_samples.deseq.vst_counts.csv'
vst_counts_by_group = 'milena_rnaseq_0818_all_samples.deseq.vst_counts_by_group.csv'

##average different groups for making heatmaps
# split_dataframe_into_groups_get_average(vst_counts, vst_counts_by_group)

##make file to graph volcano plots
##params
all_genelists = 'all_genelists.txt'
de_files = ['milena_rnaseq_0818_day22.LIG4_Mutant_vs_WT.csv', 'milena_rnaseq_0818_day4.LIG4_Mutant_vs_WT.csv', 
		'milena_rnaseq_0818_day6.LIG4_Mutant_vs_WT.csv']

##make individual genelist
genelist_files = make_ind_genelists(all_genelists)
print(genelist_files)


# de_files = ['milena_rnaseq_0818_day22.LIG4_Mutant_vs_WT.csv']
# genelist_files = ['GO_Cell_proliferation_in_forebrain.txt']

for genelist in genelist_files:
	for de_file in de_files:
		new_file = de_file.split('.')[0] + '.' + de_file.split('.')[1] + '.' + genelist.split('.')[0] + '.volcano.csv'
		add_if_in_genelist_deseq_results(de_file, new_file, genelist)





