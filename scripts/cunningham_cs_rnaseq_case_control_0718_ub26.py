#!/usr/bin/python
import filtering_annotated
import os
import subprocess
import glob
import shutil
from scipy import stats

##set input variables and parameters
delim = '\t'
##working directory
working_dir = '/data/atimms/cs_rnaseq_case_ctl_0718'
os.chdir(working_dir)



# rvis_lite_exac_out_unfiltered_suffix = '.exac.unfiltered.rvis_lite.xls'
# rvis_lite_exac_out_suffix = '.exac.rvis_lite.xls'
# exac_unfiltered_prefix = 'exac.all'
# exac_passed_prefix = 'exac.passed'
# comb_case_prefix = 'cs_rnaseq_1116.c1_and_2.cases'
# comb_ctls_prefix = 'cs_rnaseq_1116.c1_and_2.ctls'
# refgene_genename_col = 16
# count_col = 88 #start from 0
# esp6500_info_col = 89 #start from 1
# cases_passed_annotated = 'cs_rnaseq_1116.c1_and_2.cases.annotated.xls'
# refgene = '/Users/atimms/Desktop/ngs/references/annovar/hg19/hg19_refGene.txt'


##methods

def split_var_file_by_ind(infile):
	outfile = infile.rsplit('.', 1)[0] + '.split.xls'
	with open(infile, "r") as infh, open(outfile, "w") as outfh:
		line_count = 0
		for line in infh:
			line_count += 1
			line = line.rstrip().split(delim)
			if line_count == 1:
				header = line[:88] + ['sample', '\n']
				outfh.write(delim.join(header))
			else:
				samples = line[89].split(', ')
				for sample in samples:
					line_out = line[:88] + [sample, '\n']
					outfh.write(delim.join(line_out))

def make_gene_dict_from_var_file(var_file):
	gene_dict = {}
	with open(var_file, "r") as infh:
		line_count = 0
		for line in infh:
			line_count += 1
			if line_count > 1:
				line = line.rstrip().split(delim)
				gene = line[15]
				sample = line[88]
				if gene in gene_dict:
					if sample not in gene_dict[gene]:
						gene_dict[gene].append(sample)
				else:
					gene_dict[gene] = [sample]
	# for g in gene_dict:
	# 	print g, gene_dict[g]
	return gene_dict

def make_gene_pair_dict_from_var_dict(gene_dict):
	pair_dict = {}
	for gene1 in gene_dict:
		g1_samples = gene_dict[gene1]
		for g1_sample in g1_samples:
			for gene2 in gene_dict:
				if gene1 != gene2:
					g2_samples = gene_dict[gene2]
					if g1_sample in g2_samples:
						# print gene1, gene2, g1_sample
						gene_pair = gene1 + '_' + gene2
						if gene_pair in pair_dict:
							pair_dict[gene_pair].append(g1_sample)
						else:
							pair_dict[gene_pair] = [g1_sample]
	# for g in pair_dict:
	# 	print g, pair_dict[g]
	return pair_dict

def combine_case_ctl_pair_dict(case_dict, ctl_dict):
	combined_dict = {}
	for gp in case_dict:
		if gp in ctl_dict:
			combined_dict[gp] = [case_dict[gp], ctl_dict[gp]]
		else:
			combined_dict[gp] = [case_dict[gp], []]
	# for g in combined_dict:
	# 	print g, combined_dict[g]
	return combined_dict


def get_gene_pairs_cases_controls(case_file, control_file, out_file, case_total, ctl_total):
	##make dict with all samples for each gene 
	case_gene_dict = make_gene_dict_from_var_file(case_file)
	ctl_gene_dict = make_gene_dict_from_var_file(control_file)
	##get gene pairs 
	case_pair_dict = make_gene_pair_dict_from_var_dict(case_gene_dict)
	ctl_pair_dict = make_gene_pair_dict_from_var_dict(ctl_gene_dict)
	##combine case and ctl dict
	combined_pair_dict = combine_case_ctl_pair_dict(case_pair_dict, ctl_pair_dict)
	##get results, do fischers exact test and write
	with open(out_file, "w") as outfh:
		header = ['gene 1', 'gene 2', 'case count', 'control count', 'p value', 'cases', 'controls', '\n']
		outfh.write(delim.join(header))
		for gene_pair in combined_pair_dict:
			genes = gene_pair.split('_')
			case_count = len(combined_pair_dict[gene_pair][0])
			##only if there are 3 or more samples with the same gene 
			if case_count >2:
				ctl_count = len(combined_pair_dict[gene_pair][1])
				cases = '_'.join(combined_pair_dict[gene_pair][0])
				ctls = '_'.join(combined_pair_dict[gene_pair][1])
				case_counts = [case_count, case_total - case_count]
				ctl_counts = [ctl_count, ctl_total - ctl_count]
				oddsratio, pvalue = stats.fisher_exact([case_counts,ctl_counts])
				# print genes, case_counts, ctl_counts, pvalue
				final_line = genes + [str(case_count), str(ctl_count), str(pvalue), cases, ctls, '\n']
				outfh.write(delim.join(final_line))






##run methods
##parameters
number_of_cases = 391
number_of_ctls = 85
frequencies = [0.01, 0.02, 0.04, 0.001]
comb_case_prefix = 'cs_rnaseq_1116.c1_and_2.cases'
comb_ctls_prefix = 'cs_rnaseq_1116.c1_and_2.ctls'
gp_results_prefix = 'cs_rnaseq_1116.c1_and_2.gene_pairs'
file_prefixes = [comb_case_prefix + '.', comb_ctls_prefix + '.']
analysis_groups = ['.damaging_all', '.damaging_any', '.protein_changing']
gene_groups = ['.TWIST1_IGF1', '.RUNX2_IGF1', '.TWIST1_RUNX2']
##split varaint file so individual line per individual/variant
for freq in frequencies:
	for analysis_group in analysis_groups:
		for gene_group in gene_groups:
			for file_prefix in file_prefixes:
				infile = file_prefix + str(freq) + analysis_group + gene_group + '.xls'
				split_var_file_by_ind(infile)
			case_infile = comb_case_prefix + '.' + str(freq) + analysis_group + gene_group + '.split.xls'
			control_infile = comb_ctls_prefix + '.' + str(freq) + analysis_group + gene_group + '.split.xls'
			gene_pair_outfile = gp_results_prefix + '.' + str(freq) + analysis_group + gene_group + '.xls'
			get_gene_pairs_cases_controls(case_infile, control_infile, gene_pair_outfile, number_of_cases, number_of_ctls)
##get gene pairs per analysis
# for freq in frequencies:
# 	for file_suffix in ['.damaging_all.split.xls', '.damaging_any.split.xls', '.protein_changing.split.xls']:
# 		gene_pair_file = gp_results_prefix + '.' + str(freq) + '.' + file_suffix.split('.')[1] + '.xls'
# 		get_gene_pairs_cases_controls(comb_case_prefix + '.' + str(freq) + file_suffix, comb_ctls_prefix + '.' + str(freq) + file_suffix, gene_pair_file, number_of_cases, number_of_ctls)


# get_gene_pairs_cases_controls('cs_rnaseq_1116.c1_and_2.cases.0.001.damaging_all.split.xls', 'cs_rnaseq_1116.c1_and_2.ctls.0.001.damaging_all.split.xls', 'temp_results.xls', number_of_cases, number_of_ctls)

