#!/tools/BioBuilds-2015.04/bin/python
import sys
import subprocess
import os
import glob
# import dobyns_gemini_pipeline_v2

##parameters
delim = '\t'
threads = '16'

##programs
star = '/home/atimms/programs/STAR/bin/Linux_x86_64/STAR'
feature_counts = '/home/atimms/programs/subread-1.5.1-Linux-x86_64/bin/featureCounts'

##ref files etc
##hg19
genome_name = 'hg19'
star_index_dir = '/data/atimms/references/star/' + genome_name
fa_file = '/data/atimms/references/igenomes/' + genome_name + '/genome.fa'
gtf_file = '/data/atimms/references/igenomes/' + genome_name +'/genes.gtf'





def make_star_index_files(star_genome_dir, genome_fas, genome_gtf, threads_to_use):
	star_index = subprocess.Popen([star, '--runMode', 'genomeGenerate', '--genomeDir', star_genome_dir, '--genomeFastaFiles', genome_fas, '--sjdbGTFfile', genome_gtf, '--runThreadN', threads_to_use])
	star_index.wait()

def star_align_paired_end(sample_dict, star_genome_dir, threads_to_use):
	for sample_name in sample_dict:
		r1_fq = sample_dict[sample_name][0] + '_1.fastq.gz'
		r2_fq = sample_dict[sample_name][0] + '_2.fastq.gz'
		print sample_name, r1_fq, r2_fq
		star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', r1_fq, r2_fq, '--outSAMtype', 'BAM', 'Unsorted', '--readFilesCommand', 'zcat', '--outFileNamePrefix', sample_name + '.', '--runThreadN', threads_to_use])
		star_align.wait()

def star_align_single_end(sample_dict, star_genome_dir, threads_to_use):
	for sample_name in sample_dict:
		r1_fq = sample_dict[sample_name][0] + '_1.fastq.gz'
		print sample_name, r1_fq
		star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', r1_fq, '--outSAMtype', 'BAM', 'Unsorted', '--readFilesCommand', 'zcat', '--outFileNamePrefix', sample_name + '.', '--runThreadN', threads_to_use])
		star_align.wait()

def feature_count_paired_end_multiple_bams(bam_suffix, genome_gtf, outfile):
	bam_files = glob.glob('*' + bam_suffix)
	# print bam_files, len(bam_files)
	feat_count = subprocess.Popen([feature_counts,'-p', '-t', 'exon', '-g', 'gene_id', '-a', genome_gtf, '-o', outfile, '-T', '10', '-B'] + bam_files)
	feat_count.wait()

def feature_count_paired_end_multiple_bams_from_dict(bam_suffix, genome_gtf, outfile, sample_dict):
	bam_files = []
	for s in sample_dict:
		bam = s + '.' + bam_suffix
		bam_files.append(bam)
	# print bam_files, len(bam_files)
	feat_count = subprocess.Popen([feature_counts,'-p', '-t', 'exon', '-g', 'gene_id', '-a', genome_gtf, '-o', outfile, '-T', '10', '-B'] + bam_files)
	feat_count.wait()

def feature_count_single_end_multiple_bams_from_dict(bam_suffix, genome_gtf, outfile, sample_dict):
	bam_files = []
	for s in sample_dict:
		bam = s + '.' + bam_suffix
		bam_files.append(bam)
	# print bam_files, len(bam_files)
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
					print line
					samples = line[6:]
					# print samples
					samples = [s.split('.')[0] for s in samples]
					# print samples
					# samples = [sample_dict_cond[s][0] for s in samples]
					print samples
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
				print line
				line = line.strip('\n').split(delim)
				samples = line[6:]
				samples = [s.split('.')[0] for s in samples]
				for sample in samples:
					out_info = sample_dict[sample][1:]
					print sample, sample_dict[sample], out_info
					outph.write(delim.join([sample] + out_info + ['\n']))

def get_count_data_from_genelist(gene_list_file, counts_file, out_file):
	##make gene list
	gene_list = []
	with open(gene_list_file, "r") as gl_ph:
		for line in gl_ph:
			g = line.rstrip()
			gene_list.append(g)
	##print header and if gene 
	with open(out_file, "w") as outph, open(counts_file, "r") as cou_ph:
		line_count = 0
		for line in cou_ph:
			line_count += 1
			line = line.split(delim)
			if line_count == 2:
				outph.write(delim.join(line))
			else:
				gene = line[0]
				if gene in gene_list:
					outph.write(delim.join(line))

def combine_fq_file(samples):
	for sample in samples:
		print sample
		r1_fq = sample + '_1.fastq.gz'
		# r2_fq = sample + '_2.fastq.gz'
		r1_to_combine = glob.glob(sample + '*R1*fastq.gz')
		# r2_to_combine = glob.glob(sample + '*R2*fastq.gz')
		print r1_fq, r1_to_combine
		# print r2_fq, r2_to_combine
		with open(r1_fq, 'w') as r1_fh:
			cat_files = subprocess.Popen(['cat'] + r1_to_combine, stdout=r1_fh)
			cat_files.wait()
		# with open(r2_fq, 'w') as r2_fh:
		# 	cat_files = subprocess.Popen(['cat'] + r2_to_combine, stdout=r2_fh)
		# 	cat_files.wait()


##run methods

##make star index files for hg19 (only need to do once)
# make_star_index_files(star_index_dir, fa_file, gtf_file, threads)

##chunk of brain
##parameters
##setup working directory where results will be
working_dir = '/data/atimms/kim_cortex_rnaseq_0417/chunk'
os.chdir(working_dir)
project_name = 'kim_cortex_chunk_0417'

##so all samples with 4 pieces of info: fq name, cell_type, genotype and strain 
sample_dict_cond = {'KU1055': ['KU1055', 'KU1055', 'Control', '85', 'F', 'second'], 'KU1092': ['KU1092', 'KU1092', 'Control', '86', 'F', 'second'], 'KU1159-High': ['KU1159-High', 'KU1159', 'Control', '86', 'F', 'first'], 'KU1159-Low': ['KU1159-Low', 'KU1159', 'Control', '86', 'F', 'first'], 'KU1163': ['KU1163', 'KU1163', 'Control', '84', 'F', 'second'], 'KU1184': ['KU1184', 'KU1184', 'Control', '92', 'F', 'second'], 'KU1187': ['KU1187', 'KU1187', 'Control', '85', 'F', 'second'], 'T1952': ['T1952', 'T1952', 'Control', '81', 'M', 'second'], 'T2465-High': ['T2465-High', 'T2465', 'Control', '83', 'F', 'first'], 'T2465-Low': ['T2465-Low', 'T2465', 'Control', '83', 'F', 'first'], 'T3870': ['T3870', 'T3870', 'ET', '97', 'F', 'second'], 'T3891': ['T3891', 'T3891', 'ET', '87', 'F', 'second'], 'T3925': ['T3925', 'T3925', 'ET', '79', 'F', 'second'], 'T4175-High': ['T4175-High', 'T4175', 'ET', '87', 'F', 'first'], 'T4175-Low': ['T4175-Low', 'T4175', 'ET', '87', 'F', 'first'], 'T4197': ['T4197', 'T4197', 'ET', '90', 'F', 'second'], 'T4215-Hight': ['T4215-Hight', 'T4215', 'ET', '84', 'F', 'first'], 'T4215-Low': ['T4215-Low', 'T4215', 'ET', '84', 'F', 'first'], 'T4303': ['T4303', 'T4303', 'Control', '80', 'M', 'second'], 'T4388': ['T4388', 'T4388', 'Control', '88', 'F', 'second'], 'T4415': ['T4415', 'T4415', 'ET', '97', 'F', 'second'], 'T4444': ['T4444', 'T4444', 'ET', '83', 'M', 'second'], 'T4474': ['T4474', 'T4474', 'ET', '95', 'F', 'second'], 'T4482': ['T4482', 'T4482', 'ET', '91', 'M', 'second']}

star_bam_suffix = 'Aligned.out.bam'
feature_count_results_file = project_name + '.feature_counts.all_samples.txt'
deseq_feature_count_file = project_name + '.star_fc.counts.txt'
deseq_metadata_file = project_name + '.star_fc.metadata.txt'
metadata_header = ['sample', 'id', 'dx', 'age', 'sex', 'batch']

##combine fq files
# combine_fq_file(sample_dict_cond)

##star alignment
star_align_single_end(sample_dict_cond, star_index_dir, threads)

##featurecounts on all bams and then - set up to get correct bam i.e. first item in dict entry
feature_count_single_end_multiple_bams_from_dict(star_bam_suffix, gtf_file, feature_count_results_file, sample_dict_cond)
##format feature counts file
format_feature_counts_file(feature_count_results_file, deseq_feature_count_file)

##make metadata file
make_metadata_file_format_samplenames_order_as_counts(sample_dict_cond, deseq_metadata_file, feature_count_results_file, metadata_header)




