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
working_dir = '/data/atimms/dave_esc_0117/esc_rnaseq_0117'
os.chdir(working_dir)


##programs
star = '/home/atimms/programs/STAR/bin/Linux_x86_64/STAR'
feature_counts = '/home/atimms/programs/subread-1.5.1-Linux-x86_64/bin/featureCounts'

##ref files etc
##GRCh38
# star_index_dir = '/data/atimms/references/star/GRCh38'
# gencode_fa_files = '/data/atimms/references/gencode/GRCh38/GRCh38.primary_assembly.genome.fa'
# genecode_gtf = '/data/atimms/references/gencode/GRCh38/gencode.v25.primary_assembly.annotation.gtf'
##hg19
genome_name = 'mm10'
star_index_dir = '/data/atimms/references/star/' + genome_name
fa_file = '/data/atimms/references/igenomes/' + genome_name + '/genome.fa'
gtf_file = '/data/atimms/references/igenomes/' + genome_name +'/genes.gtf'


##parameters
project_name = 'esc_rnaseq_0117'

##so all samples with 4 pieces of info: fq name, cell_type, genotype and strain 
sample_dict_cond = {'WT_ESC_1': ['SRR1047502', 'ESC', 'WT', '129/B6'], 'WT_ESC_2': ['SRR1047503', 'ESC', 'WT', '129/B6'], 
		'WT_iPSC_1': ['SRR1047504', 'iPSC', 'WT', '129/B6'], 'piPSC_1': ['SRR1047505', 'piPSC', 'WT', '129/B6'], 
		'piPSC_2': ['SRR1047506', 'piPSC', 'WT', '129/B6'], 'Nanog_null_ESC_1': ['SRR1047507', 'ESC', 'nanog_null', '129/B6'], 
		'Nanog_null_ESC_2': ['SRR1047508', 'ESC', 'nanog_null', '129/B6'], 'Nanog_null_iPSC_G2_1': ['SRR1047509', 'iPSC', 'nanog_null', '129/B6'], 
		'Nanog_null_iPSC_G2_2': ['SRR1047510', 'iPSC', 'nanog_null', '129/B6'], 'Nanog_null_iPSC_G5_1': ['SRR1047511', 'iPSC', 'nanog_null', '129/B6'], 
		'Nanog_null_iPSC_G5_2': ['SRR1047512', 'iPSC', 'nanog_null', '129/B6'], 'WT_ESC_3': ['SRR1047513', 'ESC', 'WT', 'c57/b6'], 
		'WT_ESC_4': ['SRR1047514', 'ESC', 'WT', 'c57/b6'], 'WT_iPSC_2': ['SRR1047515', 'iPSC', 'WT', 'c57/b6'], 
		'WT_iPSC_3': ['SRR1047516', 'iPSC', 'WT', 'c57/b6'], 'WT_MEFs_1': ['SRR1047517', 'MEF', 'WT', 'c57/b6'], 'WT_MEFs_2': ['SRR1047518', 'MEF', 'WT', 'c57/b6']}


star_bam_suffix = 'Aligned.out.bam'
feature_count_results_file = project_name + '.feature_counts.all_samples.txt'
deseq_feature_count_file = project_name + '.star_fc.counts.txt'
deseq_metadata_file = project_name + '.star_fc.metadata.txt'
genelist_file = 'kim_genelist_0916.txt'
counts_from_genelist_file = project_name + '.kim_genelist_0916.xls'
metadata_header = ['gene', 'cell_type', 'genotype', 'strain']


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

##star alignment
# star_align_paired_end(sample_dict_cond, star_index_dir, threads)

##featurecounts on all bams and then 
# feature_count_paired_end_multiple_bams(star_bam_suffix, gtf_file, feature_count_results_file)
##format feature counts file
# format_feature_counts_file(feature_count_results_file, deseq_feature_count_file)

##make metadata file
# make_metadata_file_format_samplenames(sample_dict_cond, deseq_metadata_file)
make_metadata_file_format_samplenames_order_as_counts(sample_dict_cond, deseq_metadata_file, feature_count_results_file, metadata_header)

# get_count_data_from_genelist(genelist_file, feature_count_results_file, counts_from_genelist_file)






