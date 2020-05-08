#!/usr/bin/env python
import subprocess
import os
import glob
import math

##general parameters
delim = '\t'
human_mouse_file = '/data/atimms/references/human_mouse_hcop_fifteen_column_0517.txt'
gsea = '/home/atimms/programs/GSEA/gsea-3.0.jar'
kegg_gene_sets = 'gseaftp.broadinstitute.org://pub/gsea/gene_sets/c2.cp.kegg.v5.2.symbols.gmt'
biocarta_gene_sets = 'gseaftp.broadinstitute.org://pub/gsea/gene_sets/c2.cp.biocarta.v5.2.symbols.gmt'
go_gene_sets = 'gseaftp.broadinstitute.org://pub/gsea/gene_sets/c5.all.v5.2.symbols.gmt'
##methods
def rank_deseq_results_for_gsea_from_geo_microarray(deseq_de_results, ranked_de_results):
	gene_list = []
	with open(deseq_de_results, "r") as inph, open(ranked_de_results, "w") as outph:
		outph.write(delim.join(['gene', 'signed_log10_pvalue', '\n']))
		line_count = 0
		for line in inph:
			line_count += 1
			line = line.strip().split(delim)
			if line_count > 1:
				gene = line[6].translate(None, '"').split("///")[0]
				# print line
				if gene in gene_list:
					print 'gene %s seen multiple times'%gene
				else:
					if gene != '':
						gene_list.append(gene)
						log2fc = line[5].translate(None, '"')
						pvalue = line[1].translate(None, '"')
						if pvalue != 'NA':
							if pvalue == '0':
								pvalue = '1e-300'
							# print gene, log2fc, pvalue
							log10_p = abs(math.log10(float(pvalue)))
							if float(log2fc) >= 0:
								final = log10_p/1
							else:
								final = log10_p/-1
							# print gene, log2fc, pvalue, log10_p, final
							outph.write(delim.join([gene, str(final), '\n']))


def rank_deseq_results_for_gsea(deseq_de_results, ranked_de_results):
	with open(deseq_de_results, "r") as inph, open(ranked_de_results, "w") as outph:
		outph.write(delim.join(['gene', 'signed_log10_pvalue', '\n']))
		line_count = 0
		for line in inph:
			line_count += 1
			line = line.strip().split(',')
			if line_count > 1:
				gene = line[0].translate(None, '"')
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
					print gene, log2fc, pvalue, log10_p, final
					outph.write(delim.join([gene, str(final), '\n']))

def convert_mouse_to_human_ranked_file(infile, outfile, hcop_file):
	##make dict from mouse and human gene names
	gename_dict = {}
	with open(hcop_file, 'r') as hcop_fh:
		line_count = 0
		for line in hcop_fh:
			line_count += 1
			line = line.rstrip().split(delim)
			if line_count >1:
				human_symbol = line[4]
				mouse_symbol = line[11]
				# print human_symbol, mouse_symbol
				if mouse_symbol in gename_dict:
					gename_dict[mouse_symbol].append(human_symbol)
				else:
					gename_dict[mouse_symbol] = [human_symbol]
	##check dict
	# for g in gename_dict:
	# 	print g, gename_dict[g]
	##update file
	# '''
	with open(infile, 'r') as in_fh, open(outfile, 'w') as out_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			if line_count == 1:
				out_fh.write(line)
			else:
				line = line.rstrip().split(delim)
				m_gene = line[0]
				if m_gene in gename_dict:
					h_gene = gename_dict[m_gene][0]
					signed_value = line[1]
					# print h_gene, m_gene
					out_fh.write(delim.join([h_gene,signed_value]) + '\n')
				else:
					print 'gene %s not in file'%m_gene
	# '''

def run_gsea(ranked_file, gene_sets, test_id, sets_to_plot, geneset_min_max):
	##java -Xmx512m xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c2.cp.biocarta.v5.2.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk /Users/atimms/gsea_home/ranked_files/greb1l_rnaseq_0817.wt_ko.de.cap.rnk -scoring_scheme weighted -rpt_label test -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /Users/atimms/gsea_home/output/sep21 -gui false
	##java -cp ~/programs/gsea-3.0.jar -Xmx512m xtools.gsea.GseaPreranked
	run_gsea = subprocess.Popen(['java', '-cp', gsea, '-Xmx5000m', 'xtools.gsea.GseaPreranked', '-gmx', gene_sets,
			'-norm', 'meandiv', '-nperm', '1000', '-rnk', ranked_file, '-scoring_scheme', 'weighted', '-rpt_label', test_id,
			'-make_sets', 'true', '-plot_top_x', sets_to_plot, '-rnd_seed', 'timestamp' '-set_max', geneset_min_max[1],
			'-set_min', geneset_min_max[0], '-zip_report', 'false', '-out', working_dir, '-gui', 'false'])
	run_gsea.wait()























##run methods

##greb1l rnaseq 1017

##day4
'''
##parmeters
name = 'greb1l_rnaseq_d4_1017'
deseq_de_results = name + '.wt_ko.de.csv'
ranked_de_results = name + '.rnk'
ranked_de_results_human = name + '.human.rnk'
working_dir = '/data/atimms/dave_greb1l_rnaseq_1017'
##run stuff
os.chdir(working_dir)
rank_deseq_results_for_gsea(deseq_de_results, ranked_de_results)
convert_mouse_to_human_ranked_file(ranked_de_results, ranked_de_results_human, human_mouse_file)
run_gsea(ranked_de_results_human, biocarta_gene_sets, name + '.biocarta', '200', ['15','500'])
run_gsea(ranked_de_results_human, kegg_gene_sets, name + '.kegg', '200', ['15','500'])
run_gsea(ranked_de_results_human, go_gene_sets, name + '.gene_ontolology', '200', ['15','500'])
'''
##day0
##parmeters
name = 'greb1l_rnaseq_d0_1017'
deseq_de_results = name + '.wt_ko.de.csv'
ranked_de_results = name + '.rnk'
ranked_de_results_human = name + '.human.rnk'
working_dir = '/data/atimms/dave_greb1l_rnaseq_1017'
##run stuff
os.chdir(working_dir)
# rank_deseq_results_for_gsea(deseq_de_results, ranked_de_results)
# convert_mouse_to_human_ranked_file(ranked_de_results, ranked_de_results_human, human_mouse_file)
# run_gsea(ranked_de_results_human, biocarta_gene_sets, name + '.biocarta', '200', ['15','500'])
# run_gsea(ranked_de_results_human, kegg_gene_sets, name + '.kegg', '200', ['15','500'])
run_gsea(ranked_de_results_human, go_gene_sets, name + '.gene_ontolology', '200', ['15','500'])

