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
working_dir = '/data/atimms/kim_rnaseq_1016'
os.chdir(working_dir)
new_samples_fq_dir = '/data/atimms/kim_rnaseq_1016/new_fqs/'

##programs
star = '/home/atimms/programs/STAR/bin/Linux_x86_64/STAR'
feature_counts = '/home/atimms/programs/subread-1.5.1-Linux-x86_64/bin/featureCounts'

##ref files etc
##GRCh38
# star_index_dir = '/data/atimms/references/star/GRCh38'
# gencode_fa_files = '/data/atimms/references/gencode/GRCh38/GRCh38.primary_assembly.genome.fa'
# genecode_gtf = '/data/atimms/references/gencode/GRCh38/gencode.v25.primary_assembly.annotation.gtf'
##hg19
star_index_dir = '/data/atimms/references/star/hg19'
hg19_fa_files = '/data/atimms/references/igenomes/hg19/genome.fa'
hg19_gtf = '/data/atimms/references/igenomes/hg19/genes.gtf'


##parameters
project_name = 'kim_rnaseq_1016'
old_samples = ['H26360_Gr2', 'H26360_GrA-2', 'H26360_GrB-1', 'H26360_Pk2', 'H26360_PkA-2', 'H26360_PkB-1', 'H26360_W2', 
		'H26362_GrB-1', 'H26362_PkA-1', 'H26374_Gr1', 'H26374_GrB-1', 'H26374_Pk1', 'H26374_W1', 'H26712_A', 'H26712_B']
new_samples = ['13195-Wver-tube40', 'H25174-EGL1-tube3', 'H25174-EGL2-tube4', 'H25174-PK1-1-tube1', 'H25174-PK2-1-tube2', 'H25174-Wver-tube5', 'H26122-EGL1-tube33', 
		'H26122-EGL2-tube34', 'H26122-PK1-1-tube31', 'H26122-PK2-2-tube32', 'H26122-Wver-tube35', 'H26360-Whm-tube21', 'H26362-EGL1-tube9', 'H26362-PK1-1-tube6', 
		'H26362-RL-tube7', 'H26362-VZ-tube8', 'H26362-Wver-tube10', 'H26499-PK1-1rep-tube19', 'H26499-PK1-1-tube18', 'H26499-PKEGL1-tube20', 'H26566-EGL1-tube26', 
		'H26566-EGL2-tube27', 'H26566-PK1-1rep-tube23', 'H26566-PK1-1-tube22', 'H26566-RL-tube24', 'H26566-VZ-tube25', 'H26566-Wver-tube28', 'H26589-PK1-1-tube29', 
		'H26589-Wver-tube30', 'H26857-EGL1-tube13', 'H26857-PK1-1-tube11', 'H26857-RL-tube14', 'H26857-VZ-tube16', 'H26857-Whm-tube17', 'H26857-Wver-tube15', 'H26938-EGL1rep-tube38', 
		'H26938-PK1-1rep-tube37', 'H26938-PK1-1-tube36', 'H26938-Wver-tube39']
##so all samples with 3 pieces of info: new name, sample, tissue type - removed sample and updated groups
# sample_dict_cond = {'H26360_Gr2' :[ 'H26360_Gr2' , 'H26360' , 'GR' ], 'H26360_GrA-2' :[ 'H26360_GrA_2' , 'H26360' , 'GR' ], 'H26360_GrB-1' :[ 'H26360_GrB_1' , 'H26360' , 'GR' ], 
# 		'H26360_Pk2' :[ 'H26360_Pk2' , 'H26360' , 'PK' ], 'H26360_PkA-2' :[ 'H26360_PkA_2' , 'H26360' , 'PK' ], 'H26360_PkB-1' :[ 'H26360_PkB_1' , 'H26360' , 'PK' ], 
# 		'H26360_W2' :[ 'H26360_W2' , 'H26360' , 'Whole' ], 'H26362_GrB-1' :[ 'H26362_GrB_1' , 'H26362' , 'GR' ], 'H26362_PkA-1' :[ 'H26362_PkA_1' , 'H26362' , 'PK' ], 
# 		'H26374_Gr1' :[ 'H26374_Gr1' , 'H26374' , 'GR' ], 'H26374_GrB-1' :[ 'H26374_GrB_1' , 'H26374' , 'GR' ], 'H26374_Pk1' :[ 'H26374_Pk1' , 'H26374' , 'PK' ], 
# 		'H26374_W1' :[ 'H26374_W1' , 'H26374' , 'Whole' ], 'H26712_A' :[ 'H26712_A' , 'H26712' , 'Whole' ], 'H26712_B' :[ 'H26712_B' , 'H26712' , 'Whole' ], 
# 		'13195-Wver-tube40' :[ 'S13195_Wver' , 'S13195' , 'Wver' ], 'H25174-EGL1-tube3' :[ 'H25174_EGL1' , 'H25174' , 'EGL' ], 'H25174-EGL2-tube4' :[ 'H25174_EGL2' , 'H25174' , 'EGL' ], 
# 		'H25174-PK1-1-tube1' :[ 'H25174_PK1_1' , 'H25174' , 'PK' ], 'H25174-PK2-1-tube2' :[ 'H25174_PK2_1' , 'H25174' , 'PK' ], 'H25174-Wver-tube5' :[ 'H25174_Wver' , 'H25174' , 'Wver' ], 
# 		'H26122-EGL1-tube33' :[ 'H26122_EGL1' , 'H26122' , 'EGL' ], 'H26122-EGL2-tube34' :[ 'H26122_EGL2' , 'H26122' , 'EGL' ], 'H26122-PK1-1-tube31' :[ 'H26122_PK1_1' , 'H26122' , 'PK' ], 
# 		'H26122-PK2-2-tube32' :[ 'H26122_PK2_2' , 'H26122' , 'PK' ], 'H26122-Wver-tube35' :[ 'H26122_Wver' , 'H26122' , 'Wver' ], 'H26360-Whm-tube21' :[ 'H26360_Whm' , 'H26360' , 'Whm' ], 
# 		'H26362-EGL1-tube9' :[ 'H26362_EGL1' , 'H26362' , 'EGL' ], 'H26362-PK1-1-tube6' :[ 'H26362_PK1_1' , 'H26362' , 'PK' ], 'H26362-RL-tube7' :[ 'H26362_RL' , 'H26362' , 'RL' ], 
# 		'H26362-VZ-tube8' :[ 'H26362_VZ' , 'H26362' , 'VZ' ], 'H26362-Wver-tube10' :[ 'H26362_Wver' , 'H26362' , 'Wver' ], 'H26499-PK1-1rep-tube19' :[ 'H26499_PK1_1rep' , 'H26499' , 'PK' ], 
# 		'H26499-PK1-1-tube18' :[ 'H26499_PK1_1' , 'H26499' , 'PK' ], 'H26499-PKEGL1-tube20' :[ 'H26499_PKEGL1' , 'H26499' , 'PKEGL' ], 'H26566-EGL1-tube26' :[ 'H26566_EGL1' , 'H26566' , 'EGL' ], 
# 		'H26566-EGL2-tube27' :[ 'H26566_EGL2' , 'H26566' , 'EGL' ], 'H26566-PK1-1rep-tube23' :[ 'H26566_PK1_1rep' , 'H26566' , 'PK' ], 'H26566-PK1-1-tube22' :[ 'H26566_PK1_1' , 'H26566' , 'PK' ], 
# 		'H26566-RL-tube24' :[ 'H26566_RL' , 'H26566' , 'RL' ], 'H26566-VZ-tube25' :[ 'H26566_VZ' , 'H26566' , 'VZ' ], 'H26566-Wver-tube28' :[ 'H26566_Wver' , 'H26566' , 'Wver' ], 
# 		'H26589-PK1-1-tube29' :[ 'H26589_PK1_1' , 'H26589' , 'PK' ], 'H26589-Wver-tube30' :[ 'H26589_Wver' , 'H26589' , 'Wver' ], 'H26857-EGL1-tube13' :[ 'H26857_EGL1' , 'H26857' , 'EGL' ], 
# 		'H26857-PK1-1-tube11' :[ 'H26857_PK1_1' , 'H26857' , 'PK' ], 'H26857-RL-tube14' :[ 'H26857_RL' , 'H26857' , 'RL' ], 'H26857-VZ-tube16' :[ 'H26857_VZ' , 'H26857' , 'VZ' ], 
# 		'H26857-Whm-tube17' :[ 'H26857_Whm' , 'H26857' , 'Whm' ], 'H26857-Wver-tube15' :[ 'H26857_Wver' , 'H26857' , 'Wver' ], 'H26938-EGL1rep-tube38' :[ 'H26938_EGL1rep' , 'H26938' , 'EGL' ], 
# 		'H26938-PK1-1rep-tube37' :[ 'H26938_PK1_1rep' , 'H26938' , 'PK' ], 'H26938-PK1-1-tube36' :[ 'H26938_PK1_1' , 'H26938' , 'PK' ], 'H26938-Wver-tube39' :[ 'H26938_Wver' , 'H26938' , 'Wver' ]}
sample_dict_cond = {'H26360_Gr2' :[ 'H26360_Gr2' , 'H26360' , 'EGL' ], 'H26360_GrA-2' :[ 'H26360_GrA_2' , 'H26360' , 'EGL' ], 'H26360_GrB-1' :[ 'H26360_GrB_1' , 'H26360' , 'EGL' ], 
		'H26360_Pk2' :[ 'H26360_Pk2' , 'H26360' , 'PK' ], 'H26360_PkA-2' :[ 'H26360_PkA_2' , 'H26360' , 'PK' ], 'H26360_PkB-1' :[ 'H26360_PkB_1' , 'H26360' , 'PK' ], 
		'H26360_W2' :[ 'H26360_W2' , 'H26360' , 'Whole' ], 'H26362_GrB-1' :[ 'H26362_GrB_1' , 'H26362' , 'EGL' ], 'H26362_PkA-1' :[ 'H26362_PkA_1' , 'H26362' , 'PK' ], 
		'H26374_Gr1' :[ 'H26374_Gr1' , 'H26374' , 'EGL' ], 'H26374_GrB-1' :[ 'H26374_GrB_1' , 'H26374' , 'EGL' ], 'H26374_Pk1' :[ 'H26374_Pk1' , 'H26374' , 'PK' ], 
		'H26374_W1' :[ 'H26374_W1' , 'H26374' , 'Whole' ], 'H26712_A' :[ 'H26712_A' , 'H26712' , 'Whole' ], 'H26712_B' :[ 'H26712_B' , 'H26712' , 'Whole' ], 
		'H25174-EGL1-tube3' :[ 'H25174_EGL1' , 'H25174' , 'EGL' ], 'H25174-EGL2-tube4' :[ 'H25174_EGL2' , 'H25174' , 'EGL' ], 
		'H25174-PK1-1-tube1' :[ 'H25174_PK1_1' , 'H25174' , 'PK' ], 'H25174-PK2-1-tube2' :[ 'H25174_PK2_1' , 'H25174' , 'PK' ], 'H25174-Wver-tube5' :[ 'H25174_Wver' , 'H25174' , 'Whole' ], 
		'H26122-EGL1-tube33' :[ 'H26122_EGL1' , 'H26122' , 'EGL' ], 'H26122-EGL2-tube34' :[ 'H26122_EGL2' , 'H26122' , 'EGL' ], 'H26122-PK1-1-tube31' :[ 'H26122_PK1_1' , 'H26122' , 'PK' ], 
		'H26122-PK2-2-tube32' :[ 'H26122_PK2_2' , 'H26122' , 'PK' ], 'H26122-Wver-tube35' :[ 'H26122_Wver' , 'H26122' , 'Whole' ], 'H26360-Whm-tube21' :[ 'H26360_Whm' , 'H26360' , 'Whole' ], 
		'H26362-EGL1-tube9' :[ 'H26362_EGL1' , 'H26362' , 'EGL' ], 'H26362-PK1-1-tube6' :[ 'H26362_PK1_1' , 'H26362' , 'PK' ], 'H26362-RL-tube7' :[ 'H26362_RL' , 'H26362' , 'RL' ], 
		'H26362-VZ-tube8' :[ 'H26362_VZ' , 'H26362' , 'VZ' ], 'H26362-Wver-tube10' :[ 'H26362_Wver' , 'H26362' , 'Whole' ], 'H26499-PK1-1rep-tube19' :[ 'H26499_PK1_1rep' , 'H26499' , 'PK' ], 
		'H26499-PK1-1-tube18' :[ 'H26499_PK1_1' , 'H26499' , 'PK' ], 'H26499-PKEGL1-tube20' :[ 'H26499_PKEGL1' , 'H26499' , 'EGL' ], 'H26566-EGL1-tube26' :[ 'H26566_EGL1' , 'H26566' , 'EGL' ], 
		'H26566-EGL2-tube27' :[ 'H26566_EGL2' , 'H26566' , 'EGL' ], 'H26566-PK1-1rep-tube23' :[ 'H26566_PK1_1rep' , 'H26566' , 'PK' ], 'H26566-PK1-1-tube22' :[ 'H26566_PK1_1' , 'H26566' , 'PK' ], 
		'H26566-RL-tube24' :[ 'H26566_RL' , 'H26566' , 'RL' ], 'H26566-VZ-tube25' :[ 'H26566_VZ' , 'H26566' , 'VZ' ], 'H26566-Wver-tube28' :[ 'H26566_Wver' , 'H26566' , 'Whole' ], 
		'H26589-PK1-1-tube29' :[ 'H26589_PK1_1' , 'H26589' , 'PK' ], 'H26589-Wver-tube30' :[ 'H26589_Wver' , 'H26589' , 'Whole' ], 'H26857-EGL1-tube13' :[ 'H26857_EGL1' , 'H26857' , 'EGL' ], 
		'H26857-PK1-1-tube11' :[ 'H26857_PK1_1' , 'H26857' , 'PK' ], 'H26857-RL-tube14' :[ 'H26857_RL' , 'H26857' , 'RL' ], 'H26857-VZ-tube16' :[ 'H26857_VZ' , 'H26857' , 'VZ' ], 
		'H26857-Whm-tube17' :[ 'H26857_Whm' , 'H26857' , 'Whole' ], 'H26857-Wver-tube15' :[ 'H26857_Wver' , 'H26857' , 'Whole' ], 'H26938-EGL1rep-tube38' :[ 'H26938_EGL1rep' , 'H26938' , 'EGL' ], 
		'H26938-PK1-1rep-tube37' :[ 'H26938_PK1_1rep' , 'H26938' , 'PK' ], 'H26938-PK1-1-tube36' :[ 'H26938_PK1_1' , 'H26938' , 'PK' ], 'H26938-Wver-tube39' :[ 'H26938_Wver' , 'H26938' , 'Whole' ]}


all_samples = old_samples + new_samples
star_bam_suffix = 'Aligned.out.bam'
feature_count_results_file = project_name + '.feature_counts.all_samples.txt'
deseq_feature_count_file = project_name + '.star_fc.counts.txt'
deseq_metadata_file = project_name + '.star_fc.metadata.txt'
genelist_file = 'kim_genelist_0916.txt'
counts_from_genelist_file = project_name + '.kim_genelist_0916.xls'


def combine_fq_files(fq_dir, sample_names, final_dir):
	for sample_name in sample_names:
		print sample_name
		r1_fq = final_dir + '/' + sample_name + '.r1.fq.gz'
		r2_fq = final_dir + '/' + sample_name + '.r2.fq.gz'
		r1_to_combine = glob.glob(fq_dir + sample_name + '*R1*')
		r2_to_combine = glob.glob(fq_dir + sample_name + '*R2*')
		# print r1_fq, r1_to_combine
		with open(r1_fq, 'w') as r1_fh:
			cat_files = subprocess.Popen(['cat'] + r1_to_combine, stdout=r1_fh)
			cat_files.wait()
		with open(r2_fq, 'w') as r2_fh:
			cat_files = subprocess.Popen(['cat'] + r2_to_combine, stdout=r2_fh)
			cat_files.wait()

def make_star_index_files(star_genome_dir, genome_fas, genome_gtf, threads_to_use):
	star_index = subprocess.Popen([star, '--runMode', 'genomeGenerate', '--genomeDir', star_genome_dir, '--genomeFastaFiles', genome_fas, '--sjdbGTFfile', genome_gtf, '--runThreadN', threads_to_use])
	star_index.wait()

def star_align_paired_end(sample_names, star_genome_dir, threads_to_use, r1_info, r2_info):
	for sample_name in sample_names:
		r1_fq = glob.glob(sample_name + '*' + r1_info + '*')
		r2_fq = glob.glob(sample_name + '*' + r2_info + '*')
		print sample_name, r1_fq, r2_fq
		star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', r1_fq[0], r2_fq[0], '--outSAMtype', 'BAM', 'Unsorted', '--readFilesCommand', 'zcat', '--outFileNamePrefix', sample_name + '.', '--runThreadN', threads_to_use])
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
					# print line
					samples = line[6:]
					# print samples
					samples = [s.split('.')[0] for s in samples]
					# print samples
					samples = [sample_dict_cond[s][0] for s in samples]
					# print samples
					header = ['gene'] + samples + ['\n']
					outph.write(delim.join(header))
				else:
					gene = line[0]
					if gene != '':
						lineout = [gene] + line[6:] + ['\n']
						# print lineout
						outph.write(delim.join(lineout))


def make_metadata_file_format_samplenames_order_as_counts(sample_dict, outfile, counts_file):
	with open(outfile, "w") as outph, open(counts_file, "r") as cou_ph:
		outph.write(delim.join(['sample', 'individual', 'tissue', '\n']))
		line_count = 0
		for line in cou_ph:
			line_count += 1
			if line_count == 2:
				print line
				line = line.strip('\n').split(delim)
				samples = line[6:]
				samples = [s.split('.')[0] for s in samples]
				for sample in samples:
					print sample, sample_dict[sample]
					outph.write(delim.join(sample_dict[sample] + ['\n']))

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

##combine new fqs into 2 fastq files
# combine_fq_files(new_samples_fq_dir, new_samples, working_dir)


##make star index files for hg19 (only need to do once)
# make_star_index_files(star_index_dir, hg19_fa_files, hg19_gtf, threads)

##star alignment
# star_align_paired_end(new_samples, star_index_dir, threads, 'r1', 'r2')
# star_align_paired_end(old_samples, star_index_dir, threads, 'R1', 'R2')

##featurecounts on all bams and then 
# feature_count_paired_end_multiple_bams(star_bam_suffix, hg19_gtf, feature_count_results_file)
##format feature counts file
# format_feature_counts_file(feature_count_results_file, deseq_feature_count_file)

##make metadata file
# make_metadata_file_format_samplenames(sample_dict_cond, deseq_metadata_file)
##manually altered file so in the same order as the count file



##issue with sample 13195 (muscle contamination??) and updated dict tissue names, so repeat:
##featurecounts on all bams and then 
# feature_count_paired_end_multiple_bams_from_dict(star_bam_suffix, hg19_gtf, feature_count_results_file, sample_dict_cond)
##format feature counts file
# format_feature_counts_file(feature_count_results_file, deseq_feature_count_file)

##make metadata file
# make_metadata_file_format_samplenames_order_as_counts(sample_dict_cond, deseq_metadata_file, feature_count_results_file)

##analyze using R studio

##get_count_data_from_genelist

get_count_data_from_genelist(genelist_file, feature_count_results_file, counts_from_genelist_file)






