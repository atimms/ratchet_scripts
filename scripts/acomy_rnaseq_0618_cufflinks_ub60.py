#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import shutil

##parameters
delim = '\t'
threads = '16'

##setup working directory where results will be
working_dir = '/data/atimms/acomy_rnaseq_0618'
os.chdir(working_dir)

##reference files 
gtf_pre_post = ['/data/atimms/references/igenomes/', '/genes.gtf']
ref_pre_post = ['/data/atimms/references/igenomes/', '/genome']
ref_fa_pre_post = ['/data/atimms/references/igenomes/', '/genome.fa']

##programs
tophat2 = '/home/atimms/programs/tophat-2.1.1.Linux_x86_64/tophat2'
# picard_add_rg = '/Users/atimms/Desktop/ngs/programs/picard-tools-1.113/picard-tools-1.113/AddOrReplaceReadGroups.jar'
# picard_mkdup = '/Users/atimms/Desktop/ngs/programs/picard-tools-1.113/picard-tools-1.113/MarkDuplicates.jar'
# picard_reorder = '/Users/atimms/Desktop/ngs/programs/picard-tools-1.113/picard-tools-1.113/ReorderSam.jar'

##preferences
tophat_dir_suffix = '_tophat'
# cufflinks_dir_suffix = '_cufflinks'
# transcript_gtf_file = cufflinks_dir_suffix + '/transcripts.gtf'
# merged_gtf = './merged_asm/merged.gtf'
# pk_bams = './Pk1_tophat/accepted_hits.bam,./Pk2_tophat/accepted_hits.bam'
# gr_bams = './Gr1_tophat/accepted_hits.bam,./Gr2_tophat/accepted_hits.bam'
# w_bams = './W1_tophat/accepted_hits.bam,./W2_tophat/accepted_hits.bam'
# prefix = 'dwn_1015'
# cuffdiff_results = ['pk_gr', 'pk_w', 'gr_w']
# tophat_bam = 'accepted_hits.bam'
# htseq_suffix = '_htseq_counts.txt'





def run_tophat_paired_end(sample_name, fq1, fq2, genome, with_gtf):
	output_dir = sample_name + tophat_dir_suffix
	ref = ref_pre_post[0] + genome + ref_pre_post[1]
	##stranded
	# tophat_command = [tophat2, '-p', '6', '-G', gtf, '--library-type', 'fr-firststrand', '-o', output_dir, ref, fq1, fq2]
	if with_gtf == 'yes':
		gtf = gtf_pre_post[0] + genome + gtf_pre_post[1]
		##unstranded with gtf
		tophat_command = [tophat2, '-p', threads, '-G', gtf, '--library-type', 'fr-unstranded', '-o', output_dir, ref, fq1, fq2]
		tophat_se = subprocess.Popen(tophat_command)
		tophat_se.wait()
	elif with_gtf == 'no':
		##unstranded no gtf
		tophat_command = [tophat2, '-p', threads, '--library-type', 'fr-unstranded', '-o', output_dir, ref, fq1, fq2]
		tophat_se = subprocess.Popen(tophat_command)
		tophat_se.wait()

















##mouse data

##samples
sample_dict = {'ms1': ['ms1_USR17007944L_HM3YHBBXX_L6_1.fq.gz', 'ms1_USR17007944L_HM3YHBBXX_L6_2.fq.gz'], 
		'ms2': ['ms2_USR17007945L_HM3YHBBXX_L6_1.fq.gz', 'ms2_USR17007945L_HM3YHBBXX_L6_2.fq.gz'], 
		'ms3': ['ms3_USR17007946L_HM3YHBBXX_L7_1.fq.gz', 'ms3_USR17007946L_HM3YHBBXX_L7_2.fq.gz'], 
		'm2d1': ['m2d1_USR17007947L_HM3YHBBXX_L7_1.fq.gz', 'm2d1_USR17007947L_HM3YHBBXX_L7_2.fq.gz'], 
		'm2d2': ['m2d2_USR17007948L_HM3YHBBXX_L7_1.fq.gz', 'm2d2_USR17007948L_HM3YHBBXX_L7_2.fq.gz'], 
		'm2d3': ['m2d3_USR17007949L_HM3YHBBXX_L7_1.fq.gz', 'm2d3_USR17007949L_HM3YHBBXX_L7_2.fq.gz'],
		'm5d1': ['m5d1_USR17007950L_HM3YHBBXX_L7_1.fq.gz', 'm5d1_USR17007950L_HM3YHBBXX_L7_2.fq.gz'], 
		'm5d2': ['m5d2_USR17007951L_HM3YHBBXX_L7_1.fq.gz', 'm5d2_USR17007951L_HM3YHBBXX_L7_2.fq.gz'], 
		'm5d3': ['m5d3_USR17007952L_HM3YHBBXX_L7_1.fq.gz', 'm5d3_USR17007952L_HM3YHBBXX_L7_2.fq.gz']}
sample_condition_dict = {'H26360_Gr2': ['gr', 'H26360'], 'H26360_GrA-2': ['gr', 'H26360'], 'H26360_GrB-1': ['gr', 'H26360'], 'H26360_Pk2': ['pk', 'H26360'], 'H26360_PkA-2': ['pk', 'H26360'], 'H26360_PkB-1': ['pk', 'H26360'], 'H26360_W2': ['whole', 'H26360'], 'H26362_GrB-1': ['gr', 'H26362'], 'H26362_PkA-1': ['pk', 'H26362'], 'H26374_Gr1': ['gr', 'H26374'], 'H26374_GrB-1': ['gr', 'H26374'], 'H26374_Pk1': ['pk', 'H26374'], 'H26374_W1': ['whole', 'H26374'], 'H26712_A': ['whole', 'H26712'], 'H26712_B': ['whole', 'H26712']}
project_name = 'acomy_rnaseq_0618_cl_mouse'
genome_name = 'mm10'

##run methods
##tophat and cufflinks
# for sample in sample_dict:
	# run_tophat_paired_end(sample, sample_dict[sample][0], sample_dict[sample][1], genome_name, 'yes')
	# run_cufflinks(sample)



##mouse data

##samples
sample_dict = {'as1': ['as1_USR17007953L_HM3YHBBXX_L7_1.fq.gz','as1_USR17007953L_HM3YHBBXX_L7_2.fq.gz'],
		'as2': ['as2_USR17007954L_HM3YHBBXX_L8_1.fq.gz','as2_USR17007954L_HM3YHBBXX_L8_2.fq.gz'],
		'as3': ['as3_USR17007955L_HM3YHBBXX_L8_1.fq.gz','as3_USR17007955L_HM3YHBBXX_L8_2.fq.gz'],
		'a2d1': ['a2d1_USR17007956L_HM3YHBBXX_L8_1.fq.gz','a2d1_USR17007956L_HM3YHBBXX_L8_2.fq.gz'],
		'a2d2': ['a2d2_USR17007957L_HM3YHBBXX_L8_1.fq.gz','a2d2_USR17007957L_HM3YHBBXX_L8_2.fq.gz'],
		'a2d3': ['a2d3_USR17007958L_HM3YHBBXX_L8_1.fq.gz','a2d3_USR17007958L_HM3YHBBXX_L8_2.fq.gz'],
		'a5d1': ['a5d1_USR17007959L_HM3YHBBXX_L8_1.fq.gz','a5d1_USR17007959L_HM3YHBBXX_L8_2.fq.gz'],
		'a5d2': ['a5d2_USR17007960L_HM3YHBBXX_L8_1.fq.gz','a5d2_USR17007960L_HM3YHBBXX_L8_2.fq.gz'],
		'a5d3': ['a5d3_USR17007961L_HM3YHBBXX_L8_1.fq.gz','a5d3_USR17007961L_HM3YHBBXX_L8_2.fq.gz']}
		
sample_condition_dict = {'H26360_Gr2': ['gr', 'H26360'], 'H26360_GrA-2': ['gr', 'H26360'], 'H26360_GrB-1': ['gr', 'H26360'], 'H26360_Pk2': ['pk', 'H26360'], 'H26360_PkA-2': ['pk', 'H26360'], 'H26360_PkB-1': ['pk', 'H26360'], 'H26360_W2': ['whole', 'H26360'], 'H26362_GrB-1': ['gr', 'H26362'], 'H26362_PkA-1': ['pk', 'H26362'], 'H26374_Gr1': ['gr', 'H26374'], 'H26374_GrB-1': ['gr', 'H26374'], 'H26374_Pk1': ['pk', 'H26374'], 'H26374_W1': ['whole', 'H26374'], 'H26712_A': ['whole', 'H26712'], 'H26712_B': ['whole', 'H26712']}
project_name = 'acomy_rnaseq_0618_cl_acomy'
genome_name = 'acomy1'

##run methods
##tophat and cufflinks
for sample in sample_dict:
	run_tophat_paired_end(sample, sample_dict[sample][0], sample_dict[sample][1], genome_name, 'no')
	# run_cufflinks(sample)





