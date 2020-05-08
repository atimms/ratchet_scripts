#!/usr/bin/env python
import sys
import subprocess
import os
import glob

##parameters
delim = '\t'


##file names etc
ref_dir = '/data/atimms/references/'
fasta = ref_dir + 'human_g1k_v37.fasta'
exome_capture_bed = '/data/atimms/references/dobyns_exome.in_any_padded_target.1015.bed'
final_bam_suffix = '.bwa_gatk.bam'
bamlist = 'bams.list'
mh_ref = '/home/atimms/programs/MosaicHunter/resources/'
mh_rpts = mh_ref + 'all_repeats.b37.bed'
mh_common_site = mh_ref + 'WES_Agilent_50M.error_prone.b37.bed'
mh_dbsnp = mh_ref + 'dbsnp_137.b37.tsv'


##programs
mosaic_hunter = '/home/atimms/programs/MosaicHunter/build/mosaichunter.jar'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
gatk = '/home/atimms/programs/gatk_3.6/GenomeAnalysisTK.jar'

#annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/data/atimms/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,dbnsfp30a,popfreq_all_20150413,avsnp147']
av_operation = ['-operation', 'g,r,r,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput']

def variant_calling_paired_mutect2(bams, name_prefix):
	tumor_bam = bams[0]
	normal_bam = bams[1]
	vcf_temp1 = name_prefix + '.temp_1.vcf'
	vcf_temp2 = name_prefix + '.temp_2.vcf.gz'
	final_vcf = name_prefix + '.mutect.vcf.gz'
	##run mutect caller
	call_mutect = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'MuTect2', '-R', fasta, '-I:tumor', tumor_bam, '-I:normal', normal_bam, '-L', exome_capture_bed, '-o', vcf_temp1])
	call_mutect.wait()
	bgzip_run = subprocess.Popen([bgzip, vcf_temp1])
	bgzip_run.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines, and left normalize indels
	bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-O', 'z', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
	bcf_index.wait()

##annotate vcf file
def table_annovar_vcf(vcf, project_prefix):
		out_prefix = project_prefix
		command = [table_annovar] + av_buildver + [vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', out_prefix]
		annovar = subprocess.Popen(command)
		annovar.wait()


def variant_calling_trio_mosaichunter(pro_name, bams, pro_sex, pro_cov):
	##clean bam files
	clean_bams = []
	for bam in bams:
		sample_name = bam.split('.')[0]
		clean_bam = sample_name + '.temp_clean.bam'
		clean_bams.append(clean_bam)
		##add coverage info
		##samtools view -h -f 0x2 input.bam | perl -ne 'print if (/^@/||(/NM:i:(\d+)/&&$1<=4))' | samtools view -Sb - >cleaner.bam
		get_paired = subprocess.Popen(['samtools', 'view', '-h', '-f', '0x2', bam], stdout=subprocess.PIPE)
		rm_mismatched = subprocess.Popen(['perl', '-ne', 'print if (/^@/||(/NM:i:(\d+)/&&$1<=4))'], stdin=get_paired.stdout, stdout=subprocess.PIPE)
		st_view = subprocess.Popen(['samtools', 'view','-b', '-o', clean_bam], stdin=rm_mismatched.stdout)
		st_view.wait()
		bgzip_run = subprocess.Popen(['samtools', 'index', clean_bam])
		bgzip_run.wait()
	print 'clean_bams:', clean_bams
	# clean_bams = ['/data/atimms/dobyns_mosaic_test_0317/LR16-173_avm_lesion.temp_clean.bam', '/data/atimms/dobyns_mosaic_test_0317/LR16-173f.temp_clean.bam', '/data/atimms/dobyns_mosaic_test_0317/LR16-173m.temp_clean.bam']
	##run mosaic hunter in trio 'mode'
	outdir = pro_name + '_mh_trio_results'
	run_mh_trio = subprocess.Popen(['java', '-jar', mosaic_hunter, 'exome', '-P', 'input_file=' + clean_bams[0], '-P', 'mosaic_filter.father_bam_file=' + 
			clean_bams[1], '-P', 'mosaic_filter.mother_bam_file=' + clean_bams[2], '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=trio',
			 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
			 '-P', 'output_dir=' + outdir, '-P', 'syscall_filter.depth=' + pro_cov])
	run_mh_trio.wait()
	for b in clean_bams:
		os.remove(b)
		os.remove(b + '.bai')
		print 'deleting', b, b + '.bai'


def annotate_mh_output(s_name):
	mh_var_file = s_name + '_mh_trio_results/final.passed.tsv'
	avinput = s_name + '.avinput'
	##convert mh output to an avinput file
	with open(mh_var_file, "r") as infh, open(avinput, "w") as outfh:
		for line in infh:
			line = line.strip('\n').split(delim)
			chrom = line[0]
			pos = line[1]
			ref = line[2]
			if ref == line[6]:
				alt = line[8]
				alleles = line[7] + ',' + line[9]
			else:
				alt = line[6]
				alleles = line[9] + ',' + line[7]
				# print line
				# print chrom, pos, ref, alt, alleles
			post_prob = line[29]
			# print line
			print chrom, pos, ref, alt, alleles, post_prob
			outfh.write(delim.join([chrom, pos, pos, ref, alt, alleles, post_prob, '\n']))
	##run annovar
	command = [table_annovar] + av_buildver + [avinput] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', s_name]
	annovar = subprocess.Popen(command)
	annovar.wait()

def annotate_mutect_output(s_name):
	vcf = name_prefix + '.mutect.vcf.gz'
	##run annovar
	command = [table_annovar] + av_buildver + [vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', s_name]
	annovar = subprocess.Popen(command)
	annovar.wait()

def run_mosaic_variant_calling(work_dir, ped_name, bam_files, proband_sex, proband_coverage, analysis_type):
	os.chdir(work_dir)
	if analysis_type == 'paired':
		variant_calling_paired_mutect2(bam_files, ped_name)
		annotate_mutect_output(ped_name)
		
	if analysis_type == 'trio':
		sample_name = bam_files[0].split('.')[0]
		variant_calling_trio_mosaichunter(sample_name, bam_files, proband_sex, proband_coverage)
		annotate_mh_output(sample_name)

##run methods


##paired analysis with mutect
working_dir = '/data/atimms/dobyns_mosaic_test_0317'
# run_mosaic_variant_calling(working_dir, 'LR16-173_paired', ['LR16-173_avm_lesion.bwa_gatk.bam', 'LR16-173_saliva.bwa_gatk.bam'], '', '', 'paired')
run_mosaic_variant_calling(working_dir, 'LR13-356_paired', ['LR13-356_avm_lesion.bwa_gatk.bam', 'LR13-356_saliva.bwa_gatk.bam'], '', '', 'paired')


