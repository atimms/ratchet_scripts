#!/tools/biobuilds-2016.11/bin/python
import sys
import subprocess
import os
import glob

##parameters
delim = '\t'
threads = '16'
##setup working directory where results will be
# working_dir = '/data/atimms/dave_ra_rnaseq_0217'
working_dir = '/data/atimms/dave_ra_rnaseq_seq_test_0717'
os.chdir(working_dir)


##programs
star = '/home/atimms/programs/STAR/bin/Linux_x86_64/STAR'
feature_counts = '/home/atimms/programs/subread-1.5.1-Linux-x86_64/bin/featureCounts'

##ref files etc
##hg19
genome_name = 'mm10'
star_index_dir = '/data/atimms/references/star/' + genome_name
fa_file = '/data/atimms/references/igenomes/' + genome_name + '/genome.fa'
gtf_file = '/data/atimms/references/igenomes/' + genome_name +'/genes.gtf'
star_bam_suffix = 'Aligned.out.bam'




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



##run methods

##make star index files for hg19 (only need to do once)
# make_star_index_files(star_index_dir, fa_file, gtf_file, threads)

'''
##GSE65697
##parameters
project_name = 'ra_rnaseq_GSE65697_0217'
feature_count_results_file = project_name + '.feature_counts.all_samples.txt'
deseq_feature_count_file = project_name + '.star_fc.counts.txt'
deseq_metadata_file = project_name + '.star_fc.metadata.txt'
metadata_header = ['gene', 'treatment']

##so all samples with 3 pieces of info: fq name, treatment type 
sample_dict_cond = {'mESC_LIF_1': ['SRR1792526', 'none'], 'mESC_LIF_2': ['SRR1792527', 'none'], 'mESC_LIF_3': ['SRR1792528', 'none'], 
		'mESC_RA_1': ['SRR1792529', 'RA'], 'mESC_RA_2': ['SRR1792530', 'RA'], 'mESC_RA_3': ['SRR1792531', 'RA']}

##star alignment
star_align_single_end(sample_dict_cond, star_index_dir, threads)

##featurecounts on all bams and then - set up to get correct bam i.e. first item in dict entry
feature_count_single_end_multiple_bams_from_dict(star_bam_suffix, gtf_file, feature_count_results_file, sample_dict_cond)
##format feature counts file
format_feature_counts_file(feature_count_results_file, deseq_feature_count_file)

##make metadata file
make_metadata_file_format_samplenames_order_as_counts(sample_dict_cond, deseq_metadata_file, feature_count_results_file, metadata_header)
'''

'''
##GSE71802
##parameters
project_name = 'ra_rnaseq_GSE71802_0217'
feature_count_results_file = project_name + '.feature_counts.all_samples.txt'
deseq_feature_count_file = project_name + '.star_fc.counts.txt'
deseq_metadata_file = project_name + '.star_fc.metadata.txt'
metadata_header = ['gene', 'ercc', 'treatment']

##so all samples with 4 pieces of info: fq name, ercc, treatment
sample_dict_cond = {'RAday0_with_mix1_1': ['SRR2147943', 'mix1', 'before'], 'RAday0_with_mix2_1': ['SRR2147944', 'mix2', 'before'], 
		'RAday0_with_mix3_1': ['SRR2147945', 'mix3', 'before'], 'RAday0_with_mix4_1': ['SRR2147946', 'mix4', 'before'], 
		'RAday0_with_mix1_2': ['SRR2147947', 'mix1', 'before'], 'RAday0_with_mix2_2': ['SRR2147948', 'mix2', 'before'], 
		'RAday0_with_mix3_2': ['SRR2147949', 'mix3', 'before'], 'RAday0_with_mix4_2': ['SRR2147950', 'mix4', 'before'], 
		'RAday4_with_mix1_1': ['SRR2147951', 'mix1', 'after'], 'RAday4_with_mix2_1': ['SRR2147952', 'mix2', 'after'], 
		'RAday4_with_mix3_1': ['SRR2147953', 'mix3', 'after'], 'RAday4_with_mix4_1': ['SRR2147954', 'mix4', 'after'], 
		'RAday4_with_mix1_2': ['SRR2147955', 'mix1', 'after'], 'RAday4_with_mix2_2': ['SRR2147956', 'mix2', 'after'], 
		'RAday4_with_mix3_2': ['SRR2147957', 'mix3', 'after'], 'RAday4_with_mix4_2': ['SRR2147958', 'mix4', 'after']}

##star alignment
# star_align_paired_end(sample_dict_cond, star_index_dir, threads)

##featurecounts on all bams and then - set up to get correct bam i.e. first item in dict entry
# feature_count_paired_end_multiple_bams_from_dict(star_bam_suffix, gtf_file, feature_count_results_file, sample_dict_cond)
##format feature counts file
# format_feature_counts_file(feature_count_results_file, deseq_feature_count_file)

##make metadata file
make_metadata_file_format_samplenames_order_as_counts(sample_dict_cond, deseq_metadata_file, feature_count_results_file, metadata_header)
'''

'''
##GSE75616
##parameters
project_name = 'ra_rnaseq_GSE75616_0217'
feature_count_results_file = project_name + '.feature_counts.all_samples.txt'
deseq_feature_count_file = project_name + '.star_fc.counts.txt'
deseq_metadata_file = project_name + '.star_fc.metadata.txt'
metadata_header = ['gene', 'cell_type']

##so all samples with 3 pieces of info: fq name, cell type
sample_dict_cond = {'ground_state_mESC_1': ['SRR2969537', 'ground_state_mESC'], 'ground_state_mESC_2': ['SRR2969538', 'ground_state_mESC'], 
		'ground_state_mESC_3': ['SRR2969539', 'ground_state_mESC'], 'multipotent_progenitor_1': ['SRR2969543', 'multipotent_progenitor'], 
		'multipotent_progenitor_2': ['SRR2969544', 'multipotent_progenitor'], 'multipotent_progenitor_3': ['SRR2969545', 'multipotent_progenitor'], 
		'neural_progenitor_1': ['SRR2969549', 'neural_progenitor'], 'neural_progenitor_2': ['SRR2969550', 'neural_progenitor'], 
		'neural_progenitor_3': ['SRR2969551', 'neural_progenitor']}

##star alignment
star_align_single_end(sample_dict_cond, star_index_dir, threads)

##featurecounts on all bams and then - set up to get correct bam i.e. first item in dict entry
feature_count_single_end_multiple_bams_from_dict(star_bam_suffix, gtf_file, feature_count_results_file, sample_dict_cond)
##format feature counts file
format_feature_counts_file(feature_count_results_file, deseq_feature_count_file)

##make metadata file
make_metadata_file_format_samplenames_order_as_counts(sample_dict_cond, deseq_metadata_file, feature_count_results_file, metadata_header)
'''

##GSE65697 --- analysis 0717
##cut down to 10 and 15 million
##10M
##parameters
project_name = 'GSE65697_10m_0717'
feature_count_results_file = project_name + '.feature_counts.all_samples.txt'
deseq_feature_count_file = project_name + '.star_fc.counts.txt'
deseq_metadata_file = project_name + '.star_fc.metadata.txt'
metadata_header = ['gene', 'treatment']

##so all samples with 3 pieces of info: fq name, treatment type 
sample_dict_cond = {'mESC_LIF_1_10m': ['SRR1792526_10m', 'none'], 'mESC_LIF_2_10m': ['SRR1792527_10m', 'none'], 'mESC_LIF_3_10m': ['SRR1792528_10m', 'none'], 
		'mESC_RA_1_10m': ['SRR1792529_10m', 'RA'], 'mESC_RA_2_10m': ['SRR1792530_10m', 'RA'], 'mESC_RA_3_10m': ['SRR1792531_10m', 'RA']}

##star alignment
# star_align_single_end(sample_dict_cond, star_index_dir, threads)

##featurecounts on all bams and then - set up to get correct bam i.e. first item in dict entry
feature_count_single_end_multiple_bams_from_dict(star_bam_suffix, gtf_file, feature_count_results_file, sample_dict_cond)
##format feature counts file
format_feature_counts_file(feature_count_results_file, deseq_feature_count_file)

##make metadata file
make_metadata_file_format_samplenames_order_as_counts(sample_dict_cond, deseq_metadata_file, feature_count_results_file, metadata_header)

##15M
##parameters
project_name = 'GSE65697_15m_0717'
feature_count_results_file = project_name + '.feature_counts.all_samples.txt'
deseq_feature_count_file = project_name + '.star_fc.counts.txt'
deseq_metadata_file = project_name + '.star_fc.metadata.txt'
metadata_header = ['gene', 'treatment']

##so all samples with 3 pieces of info: fq name, treatment type 
sample_dict_cond = {'mESC_LIF_1_15m': ['SRR1792526_15m', 'none'], 'mESC_LIF_2_15m': ['SRR1792527_15m', 'none'], 'mESC_LIF_3_15m': ['SRR1792528_15m', 'none'], 
		'mESC_RA_1_15m': ['SRR1792529_15m', 'RA'], 'mESC_RA_2_15m': ['SRR1792530_15m', 'RA'], 'mESC_RA_3_15m': ['SRR1792531_15m', 'RA']}

##star alignment
# star_align_single_end(sample_dict_cond, star_index_dir, threads)

##featurecounts on all bams and then - set up to get correct bam i.e. first item in dict entry
feature_count_single_end_multiple_bams_from_dict(star_bam_suffix, gtf_file, feature_count_results_file, sample_dict_cond)
##format feature counts file
format_feature_counts_file(feature_count_results_file, deseq_feature_count_file)

##make metadata file
make_metadata_file_format_samplenames_order_as_counts(sample_dict_cond, deseq_metadata_file, feature_count_results_file, metadata_header)

