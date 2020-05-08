#!/tools/BioBuilds-2015.04/bin/python
import sys
import subprocess
import os
import glob

##parameters
delim = '\t'


##file names etc
ref_dir = '/data/atimms/references/'
fasta = ref_dir + 'human_g1k_v37.fasta'
indels_mills = ref_dir + 'Mills_and_1000G_gold_standard.indels.b37.vcf'
indels_1000g = ref_dir + '1000G_phase1.indels.b37.vcf'
rtc_intervals = ref_dir + 'Mills_1000G.intervals'
dbsnp = ref_dir + 'dbsnp_138.b37.vcf'
exome_capture_bed = '/data/atimms/references/dobyns_exome.in_any_padded_target.1015.bed'
final_bam_suffix = '.bwa_gatk.bam'
bamlist = 'bams.list'
exome_bed_for_coverage = 'hg19_ccds_merged_nochr.bed'


##programs
samtools = '/home/atimms/programs/samtools-1.2/samtools'
bedtools = '/home/atimms/programs/bedtools2-master/bin/bedtools'
picard = '/home/atimms/programs/picard-tools-1.139/picard.jar'
bwa = '/home/atimms/programs/bwa-0.7.12/bwa'
gatk = '/home/atimms/programs/GenomeAnalysisTK.jar'
bcftools = '/home/atimms/programs/bcftools-1.2/bcftools'
freebayes = '/home/atimms/programs/freebayes/bin/freebayes'
bgzip = '/home/atimms/programs/samtools-bcftools-htslib-1.0_x64-linux/bin/bgzip'
plink = '/home/atimms/programs/plink'
scalpel_discovery = '/home/atimms/programs/scalpel-0.5.2/scalpel-discovery'
scalpel_export = '/home/atimms/programs/scalpel-0.5.2/scalpel-export'

##get exome coverage for bam files
def calculate_exome_coverage(samples, bam_suffix, exome_bed, prefix):
	##parameters and header
	coverage_required = [1,5,10,20,50]
	cov_head = []
	for cr in coverage_required:
		cr1 = 'total bp >=' + str(cr)
		cr2 = 'percentage covered >=' + str(cr)
		cov_head.extend([cr1, cr2])
	header = delim.join(['bam file', 'target size', 'total bp', 'coverage'] + cov_head + ['\n'])
	##open results file and write header
	output_file = prefix + '.coverage_data.txt'
	print output_file
	with open(output_file, "w") as outfh:
		outfh.write(header)
		##go through bam files and get coverage histogram
		for sample in samples:
			bam = sample + bam_suffix
			print 'calculating coverge for bam file', bam
			##get temp coverage files
			coverage_histogram = sample + '.hist.temp'
			hist_fh = open(coverage_histogram, 'w')
			# bedtools_cov = subprocess.Popen(['bedtools', 'coverage', '-a', exome_bed, '-b', bam, '-hist'], stdout=hist_fh)
			# bedtools_cov.wait()
			bedtools_cov = subprocess.Popen(['bedtools', 'coverage', '-a', exome_bed, '-b', bam, '-hist', '-sorted'], stdout=subprocess.PIPE)
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

# for prefix in ['LR16-412']:
# for prefix in ['LR14-088']:
for prefix in ['LR16-306', 'LR14-088f']:

	sample_dict = [prefix]
	calculate_exome_coverage(sample_dict, final_bam_suffix, exome_bed_for_coverage, prefix)

