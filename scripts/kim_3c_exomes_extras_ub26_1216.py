#!/usr/bin/env python
import sys
import subprocess
import os
import glob

##parameters
delim = '\t'
working_dir = '/data/atimms/3c_exomes_1216'
os.chdir(working_dir)

##files etc
ref_dir = '/data/atimms/references/'
fasta = ref_dir + 'human_g1k_v37.fasta'

##programs
gatk = '/home/atimms/programs/gatk_3.6/GenomeAnalysisTK.jar'
tabix = '/home/atimms/programs/htslib-1.3/tabix'
bcftools = '/home/atimms/programs/bcftools-1.3/bcftools'
bgzip = '/home/atimms/programs/htslib-1.3/bgzip'

def merge_vcf_files(vcfs, out_vcf):
	vcfs_to_combine_cmd = []
	for vcf in vcfs:
		tabix_vcf = subprocess.Popen([tabix, vcf])
		tabix_vcf.wait()
		vcfs_to_combine_cmd.extend(['--variant', vcf])
	combine_var = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineVariants', '-R', fasta, '-nt', '15'] + vcfs_to_combine_cmd + ['-o', out_vcf, '-genotypeMergeOptions', 'UNSORTED', '-filteredAreUncalled'])
	combine_var.wait()

def filter_vcf(in_vcf, regions_bed, out_vcf):
	bgzip_vcf = subprocess.Popen([bgzip, in_vcf])
	bgzip_vcf.wait()
	tabix_vcf = subprocess.Popen([tabix, '-p', 'vcf', in_vcf + '.gz'])
	tabix_vcf.wait()
	bcftools_filter = subprocess.Popen([bcftools, 'view', '-R', regions_bed, '-o', out_vcf, '-O', 'z', in_vcf + '.gz'])
	bcftools_filter.wait()
	tabix_vcf = subprocess.Popen([tabix, '-p', 'vcf', out_vcf])
	tabix_vcf.wait()

def make_list_of_bams(bam_files, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for bam in bam_files:
			outfh.write(bam + '\n')

def run_gatk_ug(pos_vcf, out_vcf, bams, bamlist, regions_bed):
	make_list_of_bams(bams, bamlist)
	gatk_ug = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'UnifiedGenotyper', '-R', fasta, '-gt_mode', 'GENOTYPE_GIVEN_ALLELES', '-alleles', pos_vcf, '-L', regions_bed, '-out_mode', 'EMIT_ALL_SITES', '-I', bamlist, '-o', out_vcf])
	gatk_ug.wait()



##run methods

vcfs_to_combine = glob.glob('*ntersected_vcfs/0002.vcf.gz')
combined_vcf = 'temp_combined.vcf'
regions_file = 'regions_of_interest.bed'
regions_vcf = 'temp_regions_of_interest.vcf.gz'
genotyped_vcf = '3c_genotyped_vars.122116.vcf'
bamfiles = glob.glob('*bwa_gatk.bam')
file_of_bams = 'bam.list'

##combine all 3c vars into one vcf
# merge_vcf_files(vcfs_to_combine, combined_vcf)
##filter that vcf for variants kim is interested in
filter_vcf(combined_vcf, regions_file, regions_vcf)
run_gatk_ug(regions_vcf, genotyped_vcf, bamfiles, file_of_bams, regions_file)