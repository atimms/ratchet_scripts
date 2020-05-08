#!/usr/bin/env python
import sys
import subprocess
import os
import glob

##parameters
delim = '\t'
threads = '16'

##setup working directory where results will be
working_dir = '/data/atimms/kim_rnaseq_0717'
os.chdir(working_dir)

##programs
rnaseqc = '/home/atimms/programs/RNA-SeQC_v1.1.8.jar'
##ref files etc
genome_name = 'hg19'
star_ref_dir = '/data/atimms/references/star/'
fa_file = '/data/atimms/references/igenomes/' + genome_name + '/genome.fa'
gtf_file = '/data/atimms/references/igenomes/' + genome_name +  '/genes.gtf'
rrna_gtf = 'hg19_rRNA.gtf'
##general info
tissue_dict = {'H26857_Wver': 'Bulk', 'H26122_Wver': 'Bulk', 'H25174_Wver': 'Bulk', 'H26857_Whm': 'Bulk', 'H26712_B': 'Bulk', 'H26712_A': 'Bulk', 'H26938_Wver': 'Bulk', 'H26589_Wver': 'Bulk', 'H26374_W1': 'Bulk', 'H26360_W2': 'Bulk', 'H26566_Wver': 'Bulk', 'H26362_Wver': 'Bulk', 'H26360_Whm': 'Bulk', 'H26360_GrB_1': 'EGL', 'H26360_GrA_2': 'EGL', 'H26362_GrB_1': 'EGL', 'H26499_PKEGL1': 'EGL', 'H25174_EGL2': 'EGL', 'H26374_Gr1': 'EGL', 'H25174_EGL1': 'EGL', 'H26122_EGL1': 'EGL', 'H26374_GrB_1': 'EGL', 'H26360_Gr2': 'EGL', 'H26938_EGL1rep': 'EGL', 'H26566_EGL2': 'EGL', 'H26857_EGL1': 'EGL', 'H26566_EGL1': 'EGL', 'H26362_EGL1': 'EGL', 'H26122_EGL2': 'EGL', 'H26360_PkA_2': 'PCL', 'H26499_PK1_1rep': 'PCL', 'H26122_PK1_1': 'PCL', 'H26362_PK1_1': 'PCL', 'H26566_PK1_1': 'PCL', 'H26374_Pk1': 'PCL', 'H26360_Pk2': 'PCL', 'H26589_PK1_1': 'PCL', 'H26938_PK1_1rep': 'PCL', 'H25174_PK2_1': 'PCL', 'H26566_PK1_1rep': 'PCL', 'H26499_PK1_1': 'PCL', 'H26360_PkB_1': 'PCL', 'H25174_PK1_1': 'PCL', 'H26857_PK1_1': 'PCL', 'H26938_PK1_1': 'PCL', 'H26122_PK2_2': 'PCL', 'H26362_PkA_1': 'PCL', 'H26566_VZ': 'RL', 'H26566_RL': 'RL', 'H26362_RL': 'RL', 'H26857_RL': 'RL', 'H26857_VZ': 'RL', 'H26362_VZ': 'RL'}



##methods
def remove_non_coding_variants(infile, outfile):
	with open(infile, "r") as in_fh, open(outfile, "w") as out_fh:
		line_count, removed_count = 0, 0
		for line in in_fh:
			line_count += 1
			if line_count == 1:
				out_fh.write(line)
			else:
				line = line.split(delim)
				gene = line[0].strip()
				if gene.startswith(('MIR', 'SNOR')) or "-AS" in gene:
					print gene
					removed_count += 1
				else:
					out_fh.write(delim.join(line))
	print 'we have removed %s genes from %s total genes'%(removed_count, line_count -1)

def count_expressed_gene_per_tissue(count_file, normalized_counts_req):
	counts_dict = {}
	with open(count_file, "r") as in_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			line = line.rstrip().split(',')
			if line_count == 1:
				sample_index = [g.strip('"') for g in line[1:]]
				print sample_index
			else:
				gene = line[0].strip('"')
				sample_counts = line[1:]
				s_index = 0
				for sample_count in sample_counts:
					sample_name = sample_index[s_index]
					tissue = tissue_dict[sample_name]
					# print gene, sample_name, tissue, sample_count
					s_index += 1
					if gene in counts_dict:
						if tissue in counts_dict[gene]:
							counts_dict[gene][tissue].append(sample_count)
						else:
							counts_dict[gene][tissue] = [sample_count]
					else:
						counts_dict[gene] = {tissue:[sample_count]}
	##get gene counts
	final_count_dict = {}
	for g in counts_dict:
		print g
		print counts_dict[g]
		for tissue in counts_dict[g]:
			all_counts = counts_dict[g][tissue]
			f_all_counts = [float(i) for i in all_counts]
			average_all_counts = sum(f_all_counts) / len(f_all_counts)
			if tissue in final_count_dict:
				final_count_dict[tissue].append(average_all_counts)
			else:
				final_count_dict[tissue] = [average_all_counts]
	##now get numbers
	for t in final_count_dict:
		# print t, len(final_count_dict[t])
		total_genes = len(final_count_dict[t])
	for count_req in normalized_counts_req:
		for t in final_count_dict:
			passed_count = 0
			for av_count in final_count_dict[t]:
				if av_count >= count_req:
					passed_count += 1
			print 'for region %s out of %s genes we have %s genes that had an average normalized count of >= %s'%(t,total_genes,passed_count,count_req)

def add_if_in_genelist_deseq_results(infile, outfile, genes):
	delim = ','
	genelist = []
	##get list of genes
	with open(genes, "U") as gene_fh:
		for line in gene_fh:
			line = line.rstrip()
			genelist.append(line)
	print genelist
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
					print 'gene found:', gene
				else:
					out_fh.write(delim.join(line + ['not in genelist', '\n']))


def run_rnaseqc(sample_file, out_dir):
	rnaseqc_cmd = subprocess.Popen(['java', '-jar', rnaseqc, '-s', sample_file, '-t', gtf_file, '-r', fa_file, '-o', out_dir, '-rRNA', rrna_gtf])
	rnaseqc_cmd.wait()


def combine_rqc_rpkm_files(sample_names, outfile):
	gene_super_dict = {}
	for sample in sample_names:
		ind_gene_dict = {}
		infile = 'rnaseqc_results/' + sample + '/' + sample + '.metrics.tmp.txt.rpkm.gct'
		line_count = 0
		with open(infile, "r") as in_fh:
			for line in in_fh:
				line_count += 1
				if line_count > 3:
					line = line.rstrip().split(delim)
					gene = line[1]
					rpkm = float(line[2])
					if gene in ind_gene_dict:
						ind_gene_dict[gene].append(rpkm)
					else:
						ind_gene_dict[gene] = [rpkm]
		for g in ind_gene_dict:
			# print g, ind_gene_dict[g]
			ave_rpkm = sum(ind_gene_dict[g]) / len(ind_gene_dict[g])
			if g in gene_super_dict:
				gene_super_dict[g].append(str(ave_rpkm))
			else:
				gene_super_dict[g] = [str(ave_rpkm)]

	print 'dict contains %s genes'%len(gene_super_dict)


	with open(outfile, "w") as out_fh:
		out_fh.write(delim.join(['gene'] + sample_names + ['\n']))
		for g in gene_super_dict:
			if g != "":
				print g, gene_super_dict[g]
				out_fh.write(delim.join([g] + gene_super_dict[g] + ['\n']))


def get_samples_from_bs_rpkm(samples_wanted, infile, outfile):
	with open(infile, "r") as in_fh, open(outfile, "w") as out_fh:
		line_count = 0
		for line in in_fh:
			line = line.rstrip().split(delim)
			line_count += 1

			if line_count ==1:
				positions_wanted = [0]
				header = line
				for sample_wanted in samples_wanted:
					index_no = [i for i,x in enumerate(header) if x == sample_wanted]
					positions_wanted.extend(index_no)
				print positions_wanted, len(positions_wanted), len(samples_wanted)
				line_out = []
				for p in positions_wanted:
					line_out.append(line[p])
				out_fh.write(delim.join(line_out + ['\n']))
			else:
				line_out = []
				for p in positions_wanted:
					line_out.append(line[p])
				out_fh.write(delim.join(line_out + ['\n']))

def combine_ra_bs_rpkm_files(ra_rpkm, bs_rpkm, final_rpkm):
	rpkm_dict = {}
	header = []
	##add rnaaccess data
	with open(ra_rpkm, "r") as ra_fh:
		line_count = 0
		for line in ra_fh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count ==1:
				header.extend(line)
			else:
				gene = line[0]
				rpkms = line[1:]
				if gene in rpkm_dict:
					print 'gene %s already seen'%gene
				else:
					rpkm_dict[gene] = rpkms
	with open(bs_rpkm, "r") as bs_fh:
		line_count = 0
		for line in bs_fh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count ==1:
				header.extend(line[1:])
			else:
				gene = line[0]
				rpkms = line[1:]
				if gene in rpkm_dict:
					rpkm_dict[gene].extend(rpkms)
	with open(final_rpkm, "w") as out_fh:
		out_fh.write(delim.join(header + ['\n']))
		for g in rpkm_dict:
			print g, rpkm_dict[g], len(rpkm_dict[g])
			if len(rpkm_dict[g]) == 65:
				out_fh.write(delim.join([g] + rpkm_dict[g] + ['\n']))
			else:
				print 'gene %s had %s entries'%(g, len(rpkm_dict[g]))




##run methods

##remove non-coding genes
# count_files = ['kim_rnaseq_1016.star_fc.3_tissues.counts.txt', 'kim_rnaseq_1016.star_fc.4_tissues.counts.txt', 'kim_rnaseq_1016.star_fc.egl_pcl.counts.txt']
# for count_file in count_files:
# 	coding_count_file = count_file.rsplit('.', 1)[0] + '.coding_genes.txt'
# 	remove_non_coding_variants(count_file, coding_count_file)
##count 'expressed genes'
# norm_counts = 'kim_rnaseq_0717.norm_counts.4_tissues.csv'
# count_expressed_gene_per_tissue(norm_counts, [1,5,10,15,20])

##add if in genelist to deseq de results for making a volcano plot
# add_if_in_genelist_deseq_results('kim_rnaseq_0717.pcl_egl.de.csv', 'kim_rnaseq_0717.pcl_egl.de.asd161.csv', 'ASD_161gene_list.txt')
# add_if_in_genelist_deseq_results('kim_rnaseq_0717.pcl_egl.de.csv', 'kim_rnaseq_0717.pcl_egl.de.scz156.csv', 'SCZ_156genelist.txt')

##run rna-seqc on bam to get rpkm
##manually sorted added read group to bams and made sample file
sample_info = 'rnaseqc_sample_info.txt'
outdir = 'rnaseqc_results'
# run_rnaseqc(sample_info, outdir)

##combine rpkm with brainspan
rnaseqc_samples = ['H25174_EGL1', 'H25174-EGL1-tube3', 'H25174-EGL2-tube4', 'H25174-PK1-1-tube1', 'H25174-PK2-1-tube2', 'H25174-Wver-tube5', 'H26122-EGL1-tube33', 'H26122-EGL2-tube34', 'H26122-PK1-1-tube31', 'H26122-PK2-2-tube32', 'H26122-Wver-tube35', 'H26360_Gr2', 'H26360_GrA-2', 'H26360_GrB-1', 'H26360_Pk2', 'H26360_PkA-2', 'H26360_PkB-1', 'H26360_W2', 'H26360-Whm-tube21', 'H26362-EGL1-tube9', 'H26362_GrB-1', 'H26362-PK1-1-tube6', 'H26362_PkA-1', 'H26362-RL-tube7', 'H26362-VZ-tube8', 'H26362-Wver-tube10', 'H26374_Gr1', 'H26374_GrB-1', 'H26374_Pk1', 'H26374_W1', 'H26499-PK1-1rep-tube19', 'H26499-PK1-1-tube18', 'H26499-PKEGL1-tube20', 'H26566-EGL1-tube26', 'H26566-EGL2-tube27', 'H26566-PK1-1rep-tube23', 'H26566-PK1-1-tube22', 'H26566-RL-tube24', 'H26566-VZ-tube25', 'H26566-Wver-tube28', 'H26589-PK1-1-tube29', 'H26589-Wver-tube30', 'H26712_A', 'H26712_B', 'H26857-EGL1-tube13', 'H26857-PK1-1-tube11', 'H26857-RL-tube14', 'H26857-VZ-tube16', 'H26857-Whm-tube17', 'H26857-Wver-tube15', 'H26938-EGL1rep-tube38', 'H26938-PK1-1rep-tube37', 'H26938-PK1-1-tube36', 'H26938-Wver-tube39']
rnaseqc_rpkm_file = 'kim_rnaseq_0717.all_ra_samples.rpkm.txt'
brainspan_samples_wanted = ['d12960_CB', 'd13060_CBC', 'd12834_CB', 'd12888_CB', 'd12837_CBC', 'd12880_CBC', 'd12365_CBC', 'd12886_CBC', 'd12288_CBC', 'd12295_CBC', 'd263195015_CBC']
brainspan_rpkm_file = 'brainspan_rpkm_0916.txt'
brainspan_cb_rpkm_file = 'brainspan_cb_rpkm_0717.txt'
final_rpkm_file = 'kim_brainspan_rpkm_0717.txt'
# combine_rqc_rpkm_files(rnaseqc_samples, rnaseqc_rpkm_file)
# get_samples_from_bs_rpkm(brainspan_samples_wanted, brainspan_rpkm_file, brainspan_cb_rpkm_file)
combine_ra_bs_rpkm_files(rnaseqc_rpkm_file, brainspan_cb_rpkm_file, final_rpkm_file)

