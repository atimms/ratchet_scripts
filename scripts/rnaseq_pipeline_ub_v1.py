#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import shutil

'''
expects python module loaded i.e.
module load local_python/3.6.5 

'''


##parameters
delim = '\t'
star_threads = '16'

##programs
star = '/home/atimms/programs/STAR/bin/Linux_x86_64/STAR'
feature_counts = '/home/atimms/programs/subread-1.5.1-Linux-x86_64/bin/featureCounts'

##ref files etc
star_ref_dir = '/data/atimms/references/star/'
# fa_file = '/data/atimms/references/igenomes/' + genome_name + '/genome.fa'
gtf_file_path = ['/data/atimms/references/igenomes/', '/genes.gtf']
deseq_r_template = ''


##methods
def make_star_index_files(star_genome_dir, genome_fas, genome_gtf):
	star_index = subprocess.Popen([star, '--runMode', 'genomeGenerate', '--genomeDir', star_genome_dir, '--genomeFastaFiles', genome_fas, '--sjdbGTFfile', genome_gtf, '--runThreadN', star_threads])
	star_index.wait()


def combine_fqs(fqs_to_combine, fq_name):
	with open(fq_name, 'w') as cfq_fh:
		print("combining %s files named %s to make file %s"%(len(fqs_to_combine), fqs_to_combine, fq_name))
		cat_files = subprocess.Popen(['cat'] + fqs_to_combine, stdout=cfq_fh)
		cat_files.wait()

def star_align_paired_end(sample_name, fq_files, star_genome_dir):
	r1_fq = fq_files[0]
	r2_fq = fq_files[1]
	print(sample_name, r1_fq, r2_fq, star_genome_dir)
	star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', r1_fq, r2_fq, '--outSAMtype', 'BAM', 'Unsorted', '--readFilesCommand', 'zcat', '--outFileNamePrefix', sample_name + '.', '--runThreadN', star_threads])
	star_align.wait()

def star_align_single_end(sample_name, fq_file, star_genome_dir):
	print(sample_name, fq_file)
	star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', fq_file, '--outSAMtype', 'BAM', 'Unsorted', '--readFilesCommand', 'zcat', '--outFileNamePrefix', sample_name + '.', '--runThreadN', star_threads])
	star_align.wait()

def feature_count_paired_end(bam_suffix, genome_gtf, outfile, samples):
	bam_files = []
	for s in samples:
		bam = s + '.' + bam_suffix
		bam_files.append(bam)
	print(bam_files, len(bam_files))
	feat_count = subprocess.Popen([feature_counts,'-p', '-t', 'exon', '-g', 'gene_id', '-a', genome_gtf, '-o', outfile, '-T', '10', '-B'] + bam_files)
	feat_count.wait()

def feature_count_single_end(bam_suffix, genome_gtf, outfile, samples):
	bam_files = []
	for s in samples:
		bam = s + '.' + bam_suffix
		bam_files.append(bam)
	print(bam_files, len(bam_files))
	feat_count = subprocess.Popen([feature_counts,'-t', 'exon', '-g', 'gene_id', '-a', genome_gtf, '-o', outfile, '-T', '10'] + bam_files)
	feat_count.wait()

def format_feature_counts_file(infile, outfile):
	with open(outfile, "w") as outph, open(infile, "r") as inph:
		line_count = 0
		for line in inph:
			if line[0] != '#':
				line_count += 1
				line = line.strip('\n').split(delim)
				if line_count == 1:
					# print(line)
					samples = line[6:]
					samples = [s.split('.')[0] for s in samples]
					# samples = [sample_dict_cond[s][0] for s in samples]
					print(samples)
					header = ['gene'] + samples + ['\n']
					outph.write(delim.join(header))
				else:
					gene = line[0]
					if gene != '':
						lineout = [gene] + line[6:] + ['\n']
						# print lineout
						outph.write(delim.join(lineout))


def make_metadata_file_format_samplenames_order_as_counts(sample_dict, outfile, counts_file, header):
	with open(outfile, "w") as outph, open(counts_file, "r") as cou_ph:
		outph.write(delim.join(header + ['\n']))
		line_count = 0
		for line in cou_ph:
			line_count += 1
			if line_count == 2:
				print(line)
				line = line.strip('\n').split(delim)
				samples = line[6:]
				samples = [s.split('.')[0] for s in samples]
				for sample in samples:
					out_info = sample_dict[sample]
					print(sample, sample_dict[sample], out_info)
					outph.write(delim.join([sample] + out_info + ['\n']))


def call_all_rnaseq_methods(work_dir, info_file, genome, covariates):
	os.chdir(work_dir)
	project_name = info_file.split(".")[0]
	print('working on project:', project_name)
	##get info on project
	with open(info_file, "U") as infh:
		metadata_dict = {}
		line_count = 0
		sample_list = []
		for line in infh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				conditions = line[3:]
				print(conditions)
			else:
				sample = line[0]
				condition_by_sample = line[3:]
				metadata_dict[sample] = condition_by_sample
				sample_list.append(sample)
				fq1_info = line[1].split(',')
				fq2_info = line[2].split(',')
				# print fq1_info, len(fq1_info)
				# print fq2_info, len(fq2_info)
				##do we need to combine fqs
				if len(fq1_info) == 1:
					fq1 = fq1_info[0]
				elif len(fq1_info) >= 1:
					fq1 = sample + '.r1.fq.gz'
					combine_fqs(fq1_info, fq1)
				##single or paired sequencing
				if fq2_info == ['']:
					seq_type = 'single_end'
				else:
					seq_type = 'paired_end'
					##if same fqs for r1 and r2
					if len(fq1_info) == len(fq2_info):
						##combine files
						if len(fq2_info) == 1:
							fq2 = fq2_info[0]
						elif len(fq2_info) >= 1:
							fq2 = sample + '.r2.fq.gz'
							combine_fqs(fq2_info, fq2)
					else:
						print('issue with sample %s, there are different number of fq files for read 1 and 2'%sample)
				##run star
				star_index_dir = star_ref_dir + genome
				print('seq_type = ', seq_type)
				if seq_type == 'single_end':
					star_align_single_end(sample, fq1, star_index_dir)
				elif seq_type == 'paired_end':
					star_align_paired_end(sample, [fq1, fq2], star_index_dir)

	##run feature count on bams
	gtf_file = gtf_file_path[0] + genome + gtf_file_path[1]
	feature_count_results_file = project_name + '.feature_counts.all_samples.txt'
	deseq_feature_count_file = project_name + '.star_fc.counts.txt'
	star_bam_suffix = 'Aligned.out.bam'
	if seq_type == 'single_end':
		feature_count_single_end(star_bam_suffix, gtf_file, feature_count_results_file, sample_list)
	elif seq_type == 'paired_end':
		feature_count_paired_end(star_bam_suffix, gtf_file, feature_count_results_file, sample_list)
	format_feature_counts_file(feature_count_results_file, deseq_feature_count_file)
	##make metadata
	deseq_metadata_file = project_name + '.star_fc.metadata.txt'
	make_metadata_file_format_samplenames_order_as_counts(metadata_dict, deseq_metadata_file, feature_count_results_file, ['sample'] + conditions)

def count_and_make_files_for_deseq(work_dir, info_file, genome, covariates):
	os.chdir(work_dir)
	project_name = info_file.split(".")[0]
	print('working on project:', project_name)
	##get info on project
	with open(info_file, "U") as infh:
		metadata_dict = {}
		line_count = 0
		sample_list = []
		for line in infh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				conditions = line[3:]
				print(conditions)
			else:
				sample = line[0]
				condition_by_sample = line[3:]
				metadata_dict[sample] = condition_by_sample
				sample_list.append(sample)
				fq1_info = line[1].split(',')
				fq2_info = line[2].split(',')
				##single or paired sequencing
				if fq2_info == ['']:
					seq_type = 'single_end'
				else:
					seq_type = 'paired_end'
					##if same fqs for r1 and r2
					if len(fq1_info) == len(fq2_info):
						pass
					else:
						print('issue with sample %s, there are different number of fq files for read 1 and 2'%sample)


	##run feature count on bams
	gtf_file = gtf_file_path[0] + genome + gtf_file_path[1]
	feature_count_results_file = project_name + '.feature_counts.all_samples.txt'
	deseq_feature_count_file = project_name + '.star_fc.counts.txt'
	star_bam_suffix = 'Aligned.out.bam'
	if seq_type == 'single_end':
		feature_count_single_end(star_bam_suffix, gtf_file, feature_count_results_file, sample_list)
	elif seq_type == 'paired_end':
		feature_count_paired_end(star_bam_suffix, gtf_file, feature_count_results_file, sample_list)
	format_feature_counts_file(feature_count_results_file, deseq_feature_count_file)
	##make metadata
	deseq_metadata_file = project_name + '.star_fc.metadata.txt'
	make_metadata_file_format_samplenames_order_as_counts(metadata_dict, deseq_metadata_file, feature_count_results_file, ['sample'] + conditions)


def just_make_files_for_deseq(work_dir, info_file, genome, covariates):
	os.chdir(work_dir)
	project_name = info_file.split(".")[0]
	print('working on project:', project_name)
	##get info on project
	with open(info_file, "U") as infh:
		metadata_dict = {}
		line_count = 0
		sample_list = []
		for line in infh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				conditions = line[3:]
				print(conditions)
			else:
				sample = line[0]
				condition_by_sample = line[3:]
				metadata_dict[sample] = condition_by_sample
				sample_list.append(sample)
				fq1_info = line[1].split(',')
				fq2_info = line[2].split(',')
				# print fq1_info, len(fq1_info)
				# print fq2_info, len(fq2_info)
				##do we need to combine fqs
				if len(fq1_info) == 1:
					fq1 = fq1_info[0]
				elif len(fq1_info) >= 1:
					fq1 = sample + '.r1.fq.gz'
				##single or paired sequencing
				if fq2_info == ['']:
					seq_type = 'single_end'
				else:
					seq_type = 'paired_end'
					##if same fqs for r1 and r2
					if len(fq1_info) == len(fq2_info):
						##combine files
						if len(fq2_info) == 1:
							fq2 = fq2_info[0]
						elif len(fq2_info) >= 1:
							fq2 = sample + '.r2.fq.gz'
					else:
						print('issue with sample %s, there are different number of fq files for read 1 and 2'%sample)


	##run feature count on bams
	gtf_file = gtf_file_path[0] + genome + gtf_file_path[1]
	feature_count_results_file = project_name + '.feature_counts.all_samples.txt'
	deseq_feature_count_file = project_name + '.star_fc.counts.txt'
	star_bam_suffix = 'Aligned.out.bam'
	##make metadata
	deseq_metadata_file = project_name + '.star_fc.metadata.txt'
	make_metadata_file_format_samplenames_order_as_counts(metadata_dict, deseq_metadata_file, feature_count_results_file, ['sample'] + conditions)



##call methods

##kim cortex rnaseq
'''
# working_dir = '/data/atimms/kim_cortex_rnaseq_0417/chunk'
# conditions_to_use = []
##call all methods, so align, count and make filees for deseq
# call_all_rnaseq_methods(working_dir, 'chunk_all.txt', 'hg19')
# call_all_rnaseq_methods(working_dir, 'chunk_high.txt', 'hg19', ['dx', 'age'])
# call_all_rnaseq_methods(working_dir, 'chunk_low.txt', 'hg19', ['dx', 'age'])
# call_all_rnaseq_methods(working_dir, 'chunk_batch2.txt', 'hg19', ['dx', 'age'])
##just count and make files
# count_and_make_files_for_deseq(working_dir, 'chunk_high.txt', 'hg19', ['dx', 'age'])
# count_and_make_files_for_deseq(working_dir, 'chunk_all.txt', 'hg19', ['dx', 'age'])
# count_and_make_files_for_deseq(working_dir, 'chunk_batch1.txt', 'hg19', ['dx', 'age'])
working_dir = '/data/atimms/kim_cortex_rnaseq_0417/lcm'
call_all_rnaseq_methods(working_dir, 'lcm_all.txt', 'hg19', ['dx', 'age'])
'''


##eric cll rnaseq
'''
working_dir = '/data/atimms/eric_rnaseq_0417'
conditions_to_use = []
# call_all_rnaseq_methods(working_dir, 'cll_three_reps.txt', 'hg19', ['individual'])
# call_all_rnaseq_methods(working_dir, 'cll_individual.txt', 'hg19', ['individual'])
count_and_make_files_for_deseq(working_dir, 'cll_all_samples.txt', 'hg19', ['ind'])
'''

##glass rnaseq
'''
working_dir = '/data/atimms/glass_rnaseq_0717'
conditions_to_use = []
##call all methods, so align, count and make filees for deseq
call_all_rnaseq_methods(working_dir, 'glass_all.txt', 'hg19', ['sample', 'library'])
'''

##greb1l rnaseq
'''
working_dir = '/data/atimms/greb1l_rnaseq_0817'
conditions_to_use = []
##call all methods, so align, count and make filees for deseq
call_all_rnaseq_methods(working_dir, 'greb1l_all.txt', 'mm10', ['genotype'])

working_dir = '/data/atimms/dave_greb1l_rnaseq_1017'
conditions_to_use = []
##call all methods, so align, count and make filees for deseq
# call_all_rnaseq_methods(working_dir, 'greb1l_1017.txt', 'mm10', ['genotype', 'day', 'batch'])
# count_and_make_files_for_deseq(working_dir, 'greb1l_d4_1017.txt', 'mm10', ['day'])
# count_and_make_files_for_deseq(working_dir, 'greb1l_d0_1017.txt', 'mm10', ['day'])
count_and_make_files_for_deseq(working_dir, 'greb1l_d8_1017.txt', 'mm10', ['day'])
'''

##laure 1117
'''
working_dir = '/data/atimms/laura_rnaseq_1117'
conditions_to_use = []
##call all methods, so align, count and make filees for deseq
call_all_rnaseq_methods(working_dir, 'laura_1117_dys.txt', 'mm10', ['sample'])
'''

##kim 0218
working_dir = '/data/atimms/kim_rnaseq_0218'
conditions_to_use = []
##call all methods, so align, count and make filees for deseq
# call_all_rnaseq_methods(working_dir, 'kim_rnaseq_0218_all.txt', 'hg19', ['sample'])
##just use feature count and make the files
# count_and_make_files_for_deseq(working_dir, 'kim_rnaseq_0218_wep.txt', 'hg19', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'kim_rnaseq_0218_other.txt', 'hg19', ['sample'])
##just make the files after whole changed to bulk
# just_make_files_for_deseq(working_dir, 'kim_rnaseq_0218_wep.txt', 'hg19', ['sample'])
# just_make_files_for_deseq(working_dir, 'kim_rnaseq_0218_all.txt', 'hg19', ['sample'])

##now for the homebrew samples
##map etc
# call_all_rnaseq_methods(working_dir, 'kim_rnaseq_0218_homebrew.txt', 'hg19', ['sample'])
##then get DESeq files for comparison with all samples
# count_and_make_files_for_deseq(working_dir, 'kim_rnaseq_0218_homebrew_vs_all.txt', 'hg19', ['sample'])
##combined rl samples that are paired and compare to all egl pcl and bulk
# call_all_rnaseq_methods(working_dir, 'kim_rnaseq_0218_rl_combined.txt', 'hg19', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'kim_rnaseq_0218_rl_0518.txt', 'hg19', ['sample'])



##kim 0218 hipscs
working_dir = '/data/atimms/kim_hipscs_0318'
conditions_to_use = []
##run files
# call_all_rnaseq_methods(working_dir, 'kim_hipacs_0318_all.txt', 'hg19', ['sample'])
##get files for geo project
# call_all_rnaseq_methods(working_dir, 'GSE78474_0318.txt', 'hg19', ['sample'])
##so comapre hipscs with the geo data and some of kim 0218
# count_and_make_files_for_deseq(working_dir, 'kim_hipacs_0318_comparison.txt', 'hg19', ['sample'])

##acomy rnaseq 
working_dir = '/data/atimms/acomy_rnaseq_0618'
conditions_to_use = []
##call all methods, so align, count and make filees for deseq -- for the mouse experiment
# call_all_rnaseq_methods(working_dir, 'acomy_rnaseq_0618_mouse.txt', 'mm10', ['condition'])
##call all methods, so align, count and make filees for deseq -- for the acomy experiment
# call_all_rnaseq_methods(working_dir, 'acomy_rnaseq_0618_acomy.txt', 'acomy1', ['condition'])
##repeat the analysis changing a few parameters
# count_and_make_files_for_deseq(working_dir, 'acomy_rnaseq_0618_acomy_treatment.txt', 'acomy1', ['condition'])
# count_and_make_files_for_deseq(working_dir, 'acomy_rnaseq_0618_mouse_treatment.txt', 'mm10', ['condition'])
# count_and_make_files_for_deseq(working_dir, 'acomy_rnaseq_0618_acomy_reduced.txt', 'acomy1', ['condition'])


##acomy rnaseq 
working_dir = '/data/atimms/vishal_rnaseq_0618'
conditions_to_use = []
##make star index files for pig (only need to do once)
genome_to_index = 'Sscrofa11.1'
star_index_dir = star_ref_dir + genome_to_index
fa_file = '/data/atimms/references/igenomes/' + genome_to_index + '/genome.fa'
gtf_file = '/data/atimms/references/igenomes/' + genome_to_index + '/genes.gtf'
# make_star_index_files(star_index_dir, fa_file, gtf_file)

##call all methods, so align, count and make filees for deseq -- for the three different experiments
# call_all_rnaseq_methods(working_dir, 'vishal_0618_striatum.txt', 'Sscrofa11.1', ['condition'])
# call_all_rnaseq_methods(working_dir, 'vishal_0618_hippocampus.txt', 'Sscrofa11.1', ['condition'])
# call_all_rnaseq_methods(working_dir, 'vishal_0618_14samples.txt', 'Sscrofa11.1', ['condition'])
# call_all_rnaseq_methods(working_dir, 'test.txt', 'Sscrofa11.1', ['condition'])
# count_and_make_files_for_deseq(working_dir, 'vishal_0618_striatum_reduced.txt', 'Sscrofa11.1', ['condition'])
count_and_make_files_for_deseq(working_dir, 'vishal_0618_striatum_nogcsf.txt', 'Sscrofa11.1', ['condition'])











