#!/usr/bin/env python
import sys
import subprocess
import os
import glob

##parameters
delim = '\t'
coverage_required = [5,10,20]


##programs
bedtools = '/home/atimms/programs/bedtools2-master/bin/bedtools' ##version 2.25



##methods
def get_ref_files_for_genome(genome):
	##files depending on genome 
	if genome == 'mm10':
		fasta = '/data/atimms/references/mm10/mm10.fa'
		exome_bed = '/data/atimms/references/mm10/mm10_refGene_coding_exons.bed'
	else:
		print 'genome %s not supported'%genome
	return fasta, exome_bed



##get exome coverage for bam files
def calculate_exome_coverage(samples, bam_suffix, genome, prefix):
	##parameters and header
	cov_head = []
	for cr in coverage_required:
		cr1 = 'total bp >=' + str(cr)
		cr2 = 'percentage covered >=' + str(cr)
		cov_head.extend([cr1, cr2])
	header = delim.join(['bam file', 'target size', 'total bp', 'coverage'] + cov_head + ['\n'])
	##open results file and write header
	exome_bed = get_ref_files_for_genome(genome)[1]
	output_file = prefix + '.exome_coverage_data.txt'
	print output_file
	with open(output_file, "w") as outfh:
		outfh.write(header)
		##go through bam files and get coverage histogram
		for sample in samples:
			bam = sample + bam_suffix
			print 'calculating exome coverge for bam file', bam
			##get temp coverage files
			coverage_histogram = sample + '.hist.temp'
			hist_fh = open(coverage_histogram, 'w')
			# bedtools_cov = subprocess.Popen([bedtools, 'coverage', '-a', exome_bed, '-b', bam, '-hist', '-sorted'], stdout=hist_fh)
			# bedtools_cov.wait()
			bedtools_cov = subprocess.Popen([bedtools, 'coverage', '-a', exome_bed, '-b', bam, '-hist', '-sorted'], stdout=subprocess.PIPE)
			grep_cmd = subprocess.Popen(['grep', '^all'], stdin=bedtools_cov.stdout, stdout=hist_fh)
			grep_cmd.wait()
			hist_fh.close()
			##calculate numbers from histogram
			with open(coverage_histogram, "r") as cov_temp:
				seq_total = 0
				for line in cov_temp:
					line = line.strip('\n').split(delim)
					target_size = float(line[3])
					if line[1] != '0':
						seq_total += int(line[1]) *  int(line[2])
			cov_stats = []
			for cr in coverage_required:
				coverage_bp = 0
				with open(coverage_histogram, "r") as cov_temp:
					for line in cov_temp:
						line = line.strip('\n').split(delim)
						target_size = float(line[3])
						if int(line[1]) >= cr:
							coverage_bp += int(line[2])
					cov_stats.extend([str(coverage_bp), str((coverage_bp / target_size) * 100)])
			line_out = delim.join([bam, str(target_size), str(seq_total), str(seq_total / target_size)] + cov_stats + ['\n'])
			outfh.write(line_out)


def calculate_genome_coverage(samples, bam_suffix, genome, prefix):
	##parameters and header
	cov_head = []
	for cr in coverage_required:
		cr1 = 'total bp >=' + str(cr)
		cr2 = 'percentage covered >=' + str(cr)
		cov_head.extend([cr1, cr2])
	header = delim.join(['bam file', 'target size', 'total bp', 'coverage'] + cov_head + ['\n'])
	##open results file and write header
	fasta = get_ref_files_for_genome(genome)[0]
	output_file = prefix + '.genome_coverage_data.txt'
	print output_file
	with open(output_file, "w") as outfh:
		outfh.write(header)
		##go through bam files and get coverage histogram
		for sample in samples:
			bam = sample + bam_suffix
			print 'calculating genome coverge for bam file', bam
			##get temp coverage files
			coverage_histogram = sample + '.hist.temp'
			hist_fh = open(coverage_histogram, 'w')
			bedtools_cov = subprocess.Popen([bedtools, 'genomecov', '-ibam', bam, '-g', fasta], stdout=subprocess.PIPE)
			grep_cmd = subprocess.Popen(['grep', '^genome'], stdin=bedtools_cov.stdout, stdout=hist_fh)
			grep_cmd.wait()
			hist_fh.close()
			##calculate numbers from histogram
			with open(coverage_histogram, "r") as cov_temp:
				seq_total = 0
				for line in cov_temp:
					line = line.strip('\n').split(delim)
					target_size = float(line[3])
					if line[1] != '0':
						seq_total += int(line[1]) *  int(line[2])
			cov_stats = []
			for cr in coverage_required:
				coverage_bp = 0
				with open(coverage_histogram, "r") as cov_temp:
					for line in cov_temp:
						line = line.strip('\n').split(delim)
						target_size = float(line[3])
						if int(line[1]) >= cr:
							coverage_bp += int(line[2])
					cov_stats.extend([str(coverage_bp), str((coverage_bp / target_size) * 100)])
			line_out = delim.join([bam, str(target_size), str(seq_total), str(seq_total / target_size)] + cov_stats + ['\n'])
			outfh.write(line_out)

def samtools_flagstat(samples, bam_suffix):
	for sample in samples:
		output_file = sample + '.samtools_flagstat.txt'
		bam = sample + bam_suffix
		with open(output_file, "w") as outfh:
			st_flag = subprocess.Popen(['samtools', 'flagstat', bam], stdout=outfh)
			st_flag.wait()


##run methods
samplenames = ['J321_756_5', 'K402_0014_6', 'K402_0035_1', 'K416_0006_5', 'K416_0007_9', 'K416_0035_2', 'K417_0024_8', 'K421_0060_3', 'K421_0260_8', 'K433_0180_5', 'K434_0036_7', 'K434_0068_7', 'K434_0070_8']
# samplenames = ['J321_756_5', 'K402_0014_6']
bamfile_suffix = '.bwa_mkdup.bam'
genome_name = 'mm10'
final_prefix = 'test'

##exome cov
calculate_exome_coverage(samplenames, bamfile_suffix, genome_name, final_prefix)
##genome cov
calculate_genome_coverage(samplenames, bamfile_suffix, genome_name, final_prefix)
##samtools flagstats
samtools_flagstat(samplenames, bamfile_suffix)

