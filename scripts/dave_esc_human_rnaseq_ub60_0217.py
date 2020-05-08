#!/tools/BioBuilds-2015.04/bin/python
import sys
import subprocess
import os
import glob
# import dobyns_gemini_pipeline_v2

##parameters
delim = '\t'
threads = '16'
##setup working directory where results will be
working_dir = '/data/atimms/dave_esc_0117/human_esc_bulk_rnaseq_0217'
os.chdir(working_dir)


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



##run methods

##make star index files for hg19 (only need to do once)
# make_star_index_files(star_index_dir, fa_file, gtf_file, threads)

##cell type experiment
##parameters
project_name = 'human_esc_rnaseq_cell_type_0217'

##so all samples with 4 pieces of info: fq name, cell_type, genotype and strain 
sample_dict_cond = {'H1_rep1': ['SRR2977621', 'H1'], 'H1_rep2': ['SRR2977622', 'H1'], 'H1_rep3': ['SRR2977623', 'H1'], 
		'H1_rep4': ['SRR2977624', 'H1'], 'H9_rep1': ['SRR2977625', 'H9'], 'H9_rep2': ['SRR2977626', 'H9'], 'H9_rep3': ['SRR2977627', 'H9'], 
		'DEC_rep1': ['SRR2977628', 'DEC'], 'DEC_rep2': ['SRR2977629', 'DEC'], 'EC_rep1': ['SRR2977630', 'EC'], 'EC_rep2': ['SRR2977631', 'EC'], 
		'EC_rep3': ['SRR2977632', 'EC'], 'HFF_rep1': ['SRR2977633', 'HFF'], 'HFF_rep2': ['SRR2977634', 'HFF'], 'HFF_rep3': ['SRR2977635', 'HFF'], 
		'NPC_rep1': ['SRR2977636', 'NPC'], 'NPC_rep2': ['SRR2977637', 'NPC'], 'TB_rep1': ['SRR2977638', 'TB'], 'TB_rep2': ['SRR2977639', 'TB']}

star_bam_suffix = 'Aligned.out.bam'
feature_count_results_file = project_name + '.feature_counts.all_samples.txt'
deseq_feature_count_file = project_name + '.star_fc.counts.txt'
deseq_metadata_file = project_name + '.star_fc.metadata.txt'
metadata_header = ['gene', 'cell_type']

##star alignment
# star_align_single_end(sample_dict_cond, star_index_dir, threads)

##featurecounts on all bams and then - set up to get correct bam i.e. first item in dict entry
# feature_count_single_end_multiple_bams_from_dict(star_bam_suffix, gtf_file, feature_count_results_file, sample_dict_cond)
##format feature counts file
# format_feature_counts_file(feature_count_results_file, deseq_feature_count_file)

##make metadata file
# make_metadata_file_format_samplenames_order_as_counts(sample_dict_cond, deseq_metadata_file, feature_count_results_file, metadata_header)


##time course experiment
##parameters
project_name = 'human_esc_rnaseq_time_course_0217'

##so all samples with 4 pieces of info: fq name, cell_type, genotype and strain 
sample_dict_cond = {'H9_12h_rep1': ['SRR2977640', '12h'], 'H9_12h_rep2': ['SRR2977641', '12h'], 'H9_12h_rep3': ['SRR2977642', '12h'], 
		'H9_24h_rep1': ['SRR2977643', '24h'], 'H9_24h_rep2': ['SRR2977644', '24h'], 'H9_24h_rep3': ['SRR2977645', '24h'], 
		'H9_36h_rep1': ['SRR2977646', '36h'], 'H9_36h_rep2': ['SRR2977647', '36h'], 'H9_36h_rep3': ['SRR2977648', '36h'], 
		'H9_72h_rep1': ['SRR2977649', '72h'], 'H9_72h_rep2': ['SRR2977650', '72h'], 'H9_72h_rep3': ['SRR2977651', '72h'], 
		'H9_96h_rep1': ['SRR2977652', '96h'], 'H9_96h_rep2': ['SRR2977653', '96h'], 'H9_96h_rep3': ['SRR2977654', '96h']}



star_bam_suffix = 'Aligned.out.bam'
feature_count_results_file = project_name + '.feature_counts.all_samples.txt'
deseq_feature_count_file = project_name + '.star_fc.counts.txt'
deseq_metadata_file = project_name + '.star_fc.metadata.txt'
metadata_header = ['gene', 'time_point']

extra_sample_dict_cond = {'H9_12h_rep1': ['SRR2977640', '12h', 'early'], 'H9_12h_rep2': ['SRR2977641', '12h', 'early'], 'H9_12h_rep3': ['SRR2977642', '12h', 'early'], 
		'H9_24h_rep1': ['SRR2977643', '24h', 'early'], 'H9_24h_rep2': ['SRR2977644', '24h', 'early'], 'H9_24h_rep3': ['SRR2977645', '24h', 'early'], 
		'H9_36h_rep1': ['SRR2977646', '36h', 'early'], 'H9_36h_rep2': ['SRR2977647', '36h', 'early'], 'H9_36h_rep3': ['SRR2977648', '36h', 'early'], 
		'H9_72h_rep1': ['SRR2977649', '72h', 'late'], 'H9_72h_rep2': ['SRR2977650', '72h', 'late'], 'H9_72h_rep3': ['SRR2977651', '72h', 'late'], 
		'H9_96h_rep1': ['SRR2977652', '96h', 'late'], 'H9_96h_rep2': ['SRR2977653', '96h', 'late'], 'H9_96h_rep3': ['SRR2977654', '96h', 'late']}

extra_deseq_metadata_file = project_name + '.star_fc.metadata.txt'
extra_metadata_header = ['gene', 'time_point', 'timing']
##star alignment
# star_align_single_end(sample_dict_cond, star_index_dir, threads)

##featurecounts on all bams and then 
# feature_count_single_end_multiple_bams_from_dict(star_bam_suffix, gtf_file, feature_count_results_file, sample_dict_cond)
##format feature counts file
# format_feature_counts_file(feature_count_results_file, deseq_feature_count_file)

##make metadata file
# make_metadata_file_format_samplenames_order_as_counts(sample_dict_cond, deseq_metadata_file, feature_count_results_file, metadata_header)
##add early/late category for time points
make_metadata_file_format_samplenames_order_as_counts(extra_sample_dict_cond, extra_deseq_metadata_file, feature_count_results_file, extra_metadata_header)





