#!/usr/bin/env python
import subprocess
import glob
import shutil
import os
import math


delim = '\t'


def rank_deseq_results_for_gsea(working_dir, deseq_de_results, ranked_de_results):
	os.chdir(working_dir)
	with open(deseq_de_results, "r") as inph, open(ranked_de_results, "w") as outph:
		outph.write(delim.join(['gene', 'signed_log10_pvalue', '\n']))
		line_count = 0
		for line in inph:
			line_count += 1
			line = line.strip().split(',')
			if line_count > 1:
				gene = line[0].replace('"', '')
				log2fc = line[2]
				pvalue = line[5]
				if pvalue != 'NA':
					if pvalue == '0':
						pvalue = '1e-300'
					# print gene, log2fc, pvalue
					log10_p = abs(math.log10(float(pvalue)))
					if float(log2fc) >= 0:
						final = log10_p/1
					else:
						final = log10_p/-1
					print(gene, log2fc, pvalue, log10_p, final)
					outph.write(delim.join([gene, str(final), '\n']))

def rank_deseq_results_for_gsea_from_microarray(working_dir, deseq_de_results, ranked_de_results):
	os.chdir(working_dir)
	gene_list = []
	with open(deseq_de_results, "r") as inph, open(ranked_de_results, "w") as outph:
		outph.write(delim.join(['gene', 'signed_log10_pvalue', '\n']))
		line_count = 0
		for line in inph:
			line_count += 1
			line = line.strip().split(delim)
			if line_count > 1:
				gene = line[0].replace('"', '')
				##some lines with no genename
				if len(line) == 5:
					if gene not in gene_list:
						# print line
						gene_list.append(gene)
						log2fc = line[1]
						pvalue = line[4]
						if pvalue != 'NA':
							if pvalue == '0':
								pvalue = '1e-300'
							# print gene, log2fc, pvalue
							log10_p = abs(math.log10(float(pvalue)))
							if float(log2fc) >= 0:
								final = log10_p/1
							else:
								final = log10_p/-1
							print(gene, log2fc, pvalue, log10_p, final)
							outph.write(delim.join([gene, str(final), '\n']))

def rank_deseq_results_for_gsea_make_genes_upper(working_dir, deseq_de_results, ranked_de_results):
	os.chdir(working_dir)
	with open(deseq_de_results, "r") as inph, open(ranked_de_results, "w") as outph:
		outph.write(delim.join(['gene', 'signed_log10_pvalue', '\n']))
		line_count = 0
		for line in inph:
			line_count += 1
			line = line.strip().split(',')
			if line_count > 1:
				gene = line[0].replace('"', '').upper()
				log2fc = line[2]
				pvalue = line[5]
				if pvalue != 'NA':
					if pvalue == '0':
						pvalue = '1e-300'
					# print gene, log2fc, pvalue
					log10_p = abs(math.log10(float(pvalue)))
					if float(log2fc) >= 0:
						final = log10_p/1
					else:
						final = log10_p/-1
					print(gene, log2fc, pvalue, log10_p, final)
					outph.write(delim.join([gene, str(final), '\n']))


def rank_deseq_results_for_gsea_single_cell(working_dir, deseq_de_results, ranked_de_results):
	os.chdir(working_dir)
	with open(deseq_de_results, "r") as inph, open(ranked_de_results, "w") as outph:
		outph.write(delim.join(['gene', 'signed_log10_pvalue', '\n']))
		line_count = 0
		for line in inph:
			line_count += 1
			line = line.strip().split(',')
			if line_count > 1:
				# print line
				gene = line[7].replace('"', '').upper()
				wt = float(line[11])
				ko =  float(line[12])
				pvalue = line[3]
				if pvalue != 'NA':
					if pvalue == '0':
						pvalue = '1e-300'
					# print gene, log2fc, pvalue
					log10_p = abs(math.log10(float(pvalue)))
					if wt > ko:
						final = log10_p/1
					else:
						final = log10_p/-1
					print(gene, wt, ko, pvalue, log10_p, final)
					outph.write(delim.join([gene, str(final), '\n']))


wd = '/data/atimms/acomy_rnaseq_0618'
names = ['acomy_rnaseq_0618_acomy_reduced.sham_vs_day2.genename', 'acomy_rnaseq_0618_acomy_reduced.sham_vs_day5.genename', 
		'acomy_rnaseq_0618_acomy.sham_vs_day2.genename', 'acomy_rnaseq_0618_acomy.sham_vs_day5.genename', 
		'acomy_rnaseq_0618_acomy_treatment.sham_vs_treatment.genename', 'acomy_rnaseq_0618_mouse.sham_vs_day2.acomy_genes',
		'acomy_rnaseq_0618_mouse.sham_vs_day5.acomy_genes', 'acomy_rnaseq_0618_mouse_treatment.sham_vs_treatment.acomy_genes',
		'acomy_rnaseq_0618_mouse.sham_vs_day2', 'acomy_rnaseq_0618_mouse.sham_vs_day5', 'acomy_rnaseq_0618_mouse_treatment.sham_vs_treatment']

# names = ['acomy_rnaseq_0618_acomy_treatment.sham_vs_treatment.genename_unique', 'acomy_rnaseq_0618_acomy.sham_vs_day2.genename_unique', 
# 		'acomy_rnaseq_0618_acomy.sham_vs_day5.genename_unique']

for name in names:
	deseq_de_results = name + '.csv'
	ranked_de_results = name + '.rnk'
	ranked_de_results_upper = name + '.upper.rnk'
	rank_deseq_results_for_gsea(wd, deseq_de_results, ranked_de_results)
	rank_deseq_results_for_gsea_make_genes_upper(wd, deseq_de_results, ranked_de_results_upper)
