#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import re

##parameters
delim = '\t'
star_threads = '16'

##setup working directory where results will be
working_dir = '/data/atimms/cherry_nrl_splicing_1017'
os.chdir(working_dir)
##programs
star = '/home/atimms/programs/STAR/bin/Linux_x86_64/STAR'
feature_counts = '/home/atimms/programs/subread-1.5.1-Linux-x86_64/bin/featureCounts'


##ref files etc
star_ref_dir = '/data/atimms/references/star/'
# fa_file = '/data/atimms/references/igenomes/' + genome_name + '/genome.fa'
gtf_file_path = ['/data/atimms/references/igenomes/', '/genes.gtf']
deseq_r_template = ''


def combine_fqs(fqs_to_combine, fq_name):
	with open(fq_name, 'w') as cfq_fh:
		print "combining %s files named %s to make file %s"%(len(fqs_to_combine), fqs_to_combine, fq_name)
		cat_files = subprocess.Popen(['cat'] + fqs_to_combine, stdout=cfq_fh)
		cat_files.wait()

def star_align_paired_end(sample_name, fq_files, star_genome_dir):
	r1_fq = fq_files[0]
	r2_fq = fq_files[1]
	print sample_name, r1_fq, r2_fq, star_genome_dir
	star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', r1_fq, r2_fq, '--outSAMtype', 'BAM', 'Unsorted', '--readFilesCommand', 'zcat', '--outFileNamePrefix', sample_name + '.', '--runThreadN', star_threads])
	star_align.wait()

def sort_index_star_bam(sample_name, inbam_suffix, outbam_suffix):
	in_bam = sample_name + inbam_suffix
	out_bam = sample_name + outbam_suffix
	with open (out_bam, 'w') as out_fh:
		st_sort = subprocess.Popen(['samtools', 'sort','-O', 'bam', '-T', sample_name, in_bam], stdout=out_fh)
		st_sort.wait()
		st_index = subprocess.Popen(['samtools', 'index', out_bam])

##samtools view -bho t.bam Hu1_ret_whole_RNA-27964978.sorted.bam chr14:24538814-24593305
def get_region_from_bam(in_bam, out_bam, region):
	st_region = subprocess.Popen(['samtools', 'view', '-bho', out_bam, in_bam, region])
	st_region.wait()
	st_index = subprocess.Popen(['samtools', 'index', out_bam])
	st_index.wait()
##parameters to run


##run methods

##mouse data
'''
sample_dict = {'WT3': ['TJC__mm__Retina__WT3__90bp_L.fq.gz', 'TJC__mm__Retina__WT3__90bp_R.fq.gz'], 
'WT4': ['TJC__mm__Retina__WT4__90bp_L.fq.gz', 'TJC__mm__Retina__WT4__90bp_R.fq.gz']}
for sample in sample_dict:
	# star_align_paired_end(sample, sample_dict[sample], star_ref_dir + 'mm10')
	sort_index_star_bam(sample, '.Aligned.out.bam', '.sorted.bam')
'''

##human data
sample_dict = {'Hu1_ret_nuclei_RNA-27963984': [['Hu1-ret-nuclei-RNA_S4_L001_R1_001.fastq.gz', 'Hu1-ret-nuclei-RNA_S4_L002_R1_001.fastq.gz', 'Hu1-ret-nuclei-RNA_S4_L003_R1_001.fastq.gz', 'Hu1-ret-nuclei-RNA_S4_L004_R1_001.fastq.gz'], 
			['Hu1-ret-nuclei-RNA_S4_L001_R2_001.fastq.gz', 'Hu1-ret-nuclei-RNA_S4_L002_R2_001.fastq.gz', 'Hu1-ret-nuclei-RNA_S4_L003_R2_001.fastq.gz', 'Hu1-ret-nuclei-RNA_S4_L004_R2_001.fastq.gz']], 
			'Hu1_ret_whole_RNA-27964978': [['Hu1-ret-whole-RNA_S1_L001_R1_001.fastq.gz', 'Hu1-ret-whole-RNA_S1_L002_R1_001.fastq.gz', 'Hu1-ret-whole-RNA_S1_L003_R1_001.fastq.gz', 'Hu1-ret-whole-RNA_S1_L004_R1_001.fastq.gz'],
			['Hu1-ret-whole-RNA_S1_L001_R2_001.fastq.gz', 'Hu1-ret-whole-RNA_S1_L002_R2_001.fastq.gz', 'Hu1-ret-whole-RNA_S1_L003_R2_001.fastq.gz', 'Hu1-ret-whole-RNA_S1_L004_R2_001.fastq.gz']],
			'Hu2_ret_whole_RNA-27995973': [['Hu2-ret-whole-RNA_S2_L001_R1_001.fastq.gz', 'Hu2-ret-whole-RNA_S2_L002_R1_001.fastq.gz', 'Hu2-ret-whole-RNA_S2_L003_R1_001.fastq.gz', 'Hu2-ret-whole-RNA_S2_L004_R1_001.fastq.gz'],
			['Hu2-ret-whole-RNA_S2_L001_R2_001.fastq.gz', 'Hu2-ret-whole-RNA_S2_L002_R2_001.fastq.gz', 'Hu2-ret-whole-RNA_S2_L003_R2_001.fastq.gz', 'Hu2-ret-whole-RNA_S2_L004_R2_001.fastq.gz']]}
for sample in sample_dict:
	# read_1 = sample + '.r1.fq.gz'
	# read_2 = sample + '.r2.fq.gz'
	# combine_fqs(sample_dict[sample][0], read_1)
	# combine_fqs(sample_dict[sample][1], read_2)
	# star_align_paired_end(sample, [read_1, read_2], star_ref_dir + 'hg19')
	# sort_index_star_bam(sample, '.Aligned.out.bam', '.sorted.bam')
	get_region_from_bam(sample + '.sorted.bam', sample + '.sorted.nrl.bam', 'chr14:24538814-24593305')











