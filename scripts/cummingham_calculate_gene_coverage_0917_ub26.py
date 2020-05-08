#!/usr/bin/env python
import sys
import subprocess
import os
import glob

##parameters
delim = '\t'
coverage_required = [5,10,20]


##programs
bedtools = '/home/atimms/programs/bedtools2-master/bin/bedtools' 





##get exome coverage for bam files
def calculate_bed_target_coverage(samples, bam_suffix, bed_file, prefix):
	##parameters and header
	cov_head = []
	for cr in coverage_required:
		cr1 = 'total bp >=' + str(cr)
		cr2 = 'percentage covered >=' + str(cr)
		cov_head.extend([cr1, cr2])
	header = delim.join(['bam file', 'target size', 'total bp', 'coverage'] + cov_head + ['\n'])
	##open results file and write header
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
			# bedtools_cov = subprocess.Popen([bedtools, 'coverage', '-a',bed_file , '-b', bam, '-hist', '-sorted'], stdout=subprocess.PIPE)
			bedtools_cov = subprocess.Popen([bedtools, 'coverage', '-a',bed_file , '-b', bam, '-hist'], stdout=subprocess.PIPE)
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

##run methods
# samplenames = ['163502', '163505', '163508', '163509', '163525', '163548', '163557', '163559', '163561', '163562', '163563', '163582', '163594', '163598', '163603', '163616', '163628', '163646', '163655', '163658', '163662', '163663', '163700', '95493', '95521', '95599', '95682']
samplenames = ['95599']
bamfile_suffix = '.accepted_hits.merged.nodups.realigned.recal.bam'
# bamfile_suffix = '.sorted.bam'
genenames = ['FGFR1', 'FGFR2', 'FGFR3', 'TWIST1', 'TCF12', 'EFNB1']
# genenames = ['TWIST1']
for genename in genenames:
	final_prefix = 'cs_rnaseq_0917.' + genename
	bedfile = genename + '.bed'
	calculate_bed_target_coverage(samplenames, bamfile_suffix, bedfile, final_prefix)





