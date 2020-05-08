#!/usr/bin/env python
import sys
import subprocess
import os


##parameters
delim = '\t'

def make_file_for_graphing(metadata_file, deseq_file, outfile, genelist):
	##make dict from normalized counts file
	count_dict = {}
	for gene in genelist:
		with open(deseq_file, "r") as norm_fh:
			line_count = 0
			for line in norm_fh:
				line = line.rstrip().split(',')
				line_count += 1
				if line_count == 1:
					samples = line[1:]
					real_samples = [s.strip('"') for s in samples]
					# print samples, real_samples
				else:
					nc_gene = line[0].strip('"')
					if nc_gene == gene:
						print line
						ncs = line[1:]
						index_no = 0
						for nc in ncs:
							sample = real_samples[index_no]
							print sample, nc, index_no
							index_no += 1
							if sample in count_dict:
								count_dict[sample].append(nc)
							else:
								count_dict[sample] = [nc]
							
	# print count_dict.keys()
	# print len(count_dict.keys())
	##update file
	with open(metadata_file, "r") as in_fh, open(outfile, "w") as out_fh:
		line_count = 0
		for line in in_fh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				out_fh.write(delim.join([line[0], line[2]] + genelist + ['\n']))
			else:
				sample_name = line[0]
				tissue = line[2]
				normcounts = count_dict[sample_name]
				out_fh.write(delim.join([sample_name, tissue] + normcounts + ['\n']))
				# gene = line[0].strip('"')
				# # print line[0], gene
				# if gene in genelist:
				# 	out_fh.write(delim.join(line + ['in genelist', '\n']))
				# 	print 'gene found:', gene
				# else:
				# 	out_fh.write(delim.join(line + ['not in genelist', '\n']))

def make_file_for_graphing_facetted(metadata_file, deseq_file, outfile, genelist):
	##update file
	with open(metadata_file, "r") as in_fh, open(outfile, "w") as out_fh:
		line_count = 0
		for line in in_fh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				out_fh.write(delim.join([line[0], line[2], 'gene', 'count','\n']))
			else:
				sample_name = line[0]
				tissue = line[2]

				for gene in genelist:
					count_dict = {}
					with open(deseq_file, "r") as norm_fh:
						line_count = 0
						for line in norm_fh:
							line = line.rstrip().split(',')
							line_count += 1
							if line_count == 1:
								samples = line[1:]
								real_samples = [s.strip('"') for s in samples]
								# print samples, real_samples
							else:
								nc_gene = line[0].strip('"')
								if nc_gene == gene:
									print line
									ncs = line[1:]
									index_no = 0
									for nc in ncs:
										sample = real_samples[index_no]
										print sample, nc, index_no
										index_no += 1
										if sample in count_dict:
											count_dict[sample].append(nc)
										else:
											count_dict[sample] = [nc]

					normcounts = count_dict[sample_name]
					out_fh.write(delim.join([sample_name, tissue, gene] + normcounts + ['\n']))
				# gene = line[0].strip('"')
				# # print line[0], gene
				# if gene in genelist:
				# 	out_fh.write(delim.join(line + ['in genelist', '\n']))
				# 	print 'gene found:', gene
				# else:
				# 	out_fh.write(delim.join(line + ['not in genelist', '\n']))
def make_file_for_linegraph_facetted(metadata_file, deseq_file, outfile, genelist):
	##update file
	with open(metadata_file, "r") as in_fh, open(outfile, "w") as out_fh:
		line_count = 0
		for line in in_fh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				out_fh.write(delim.join([line[0], line[2], line[4], 'gene', 'count','\n']))
			else:
				sample_name = line[0]
				tissue = line[2]
				age = line[4]
				for gene in genelist:
					count_dict = {}
					with open(deseq_file, "r") as norm_fh:
						line_count = 0
						for line in norm_fh:
							line = line.rstrip().split(',')
							line_count += 1
							if line_count == 1:
								samples = line[1:]
								real_samples = [s.strip('"') for s in samples]
								# print samples, real_samples
							else:
								nc_gene = line[0].strip('"')
								if nc_gene == gene:
									# print line
									ncs = line[1:]
									index_no = 0
									for nc in ncs:
										sample = real_samples[index_no]
										# print sample, nc, index_no
										index_no += 1
										if sample in count_dict:
											count_dict[sample].append(nc)
										else:
											count_dict[sample] = [nc]

					normcounts = count_dict[sample_name]
					out_fh.write(delim.join([sample_name, tissue, age, gene] + normcounts + ['\n']))
				# gene = line[0].strip('"')
				# # print line[0], gene
				# if gene in genelist:
				# 	out_fh.write(delim.join(line + ['in genelist', '\n']))
				# 	print 'gene found:', gene
				# else:
				# 	out_fh.write(delim.join(line + ['not in genelist', '\n']))


##kim's rnaseq_1016 analzed 0417
'''
genes = ['ATOH1', 'PAX6', 'SKOR2', 'CALB1']
md_file = 'kim_rnaseq_1016.star_fc.metadata2.txt'
norm_counts = 'kim_rnaseq_1016.norm_counts.csv'
file_for_graphing = 'kim_rnaseq_1016.4_genes_0517.r_boxplot.txt'
file_for_graphing_facet = 'kim_rnaseq_1016.4_genes_0517.r_boxplot_facet.txt'
# make_file_for_graphing(md_file, norm_counts, file_for_graphing, genes)
make_file_for_graphing_facetted(md_file, norm_counts, file_for_graphing_facet, genes)
'''

##kim's rnaseq_1016 analzed 0717
genes = ['ATOH1', 'PAX6', 'SKOR2', 'CALB1']
md_file = 'kim_rnaseq_1016.star_fc.3_tissues.metadata_0717.txt'
norm_counts = 'kim_rnaseq_0717.norm_counts.4_tissues.csv'
r_log_counts = 'kim_rnaseq_0717.rlog_counts.4_tissues.csv'
file_for_graphing = 'kim_rnaseq_0717.4_genes.r_boxplot.txt'
file_for_graphing_facet_norm = 'kim_rnaseq_0717.4_genes.norm.r_boxplot_facet.txt'
file_for_graphing_facet_rlog = 'kim_rnaseq_0717.4_genes.rlog.r_boxplot_facet.txt'
file_for_line_facet_norm = 'kim_rnaseq_0717.4_genes.norm.r_line_facet.txt'
file_for_line_facet_rlog = 'kim_rnaseq_0717.4_genes.rlog.r_line_facet.txt'
# make_file_for_graphing_facetted(md_file, norm_counts, file_for_graphing_facet_norm, genes)
# make_file_for_graphing_facetted(md_file, r_log_counts, file_for_graphing_facet_rlog, genes)
make_file_for_linegraph_facetted(md_file, norm_counts, file_for_line_facet_norm, genes)
# make_file_for_linegraph_facetted(md_file, r_log_counts, file_for_line_facet_rlog, genes)
