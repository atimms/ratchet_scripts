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
av_options = ['-otherinfo', '-remove', '-nastring', '.','-arg', '-splicing 10 ,,,,,']
av_options_vcf = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput','-arg', '-splicing 10 ,,,,,']

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
	outfiles = []
	##clean bam files
	clean_bams = []
	for bam in bams:
		sample_name = bam.split('.')[0]
		clean_bam = sample_name + '.temp_clean.bam'
		clean_bams.append(clean_bam)
		##clean the bam
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
	##default
	outdir = pro_name + '_mh_trio_results_default'
	run_mh_trio = subprocess.Popen(['java', '-jar', mosaic_hunter, 'exome', '-P', 'input_file=' + clean_bams[0], '-P', 'mosaic_filter.father_bam_file=' + 
			clean_bams[1], '-P', 'mosaic_filter.mother_bam_file=' + clean_bams[2], '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=trio',
			 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
			 '-P', 'output_dir=' + outdir, '-P', 'syscall_filter.depth=' + pro_cov])
	run_mh_trio.wait()
	outfiles.append(outdir + '/final.passed.tsv')
	##lower cov
	outdir = pro_name + '_mh_trio_results_lowcov'
	run_mh_trio = subprocess.Popen(['java', '-jar', mosaic_hunter, 'exome', '-P', 'input_file=' + clean_bams[0], '-P', 'mosaic_filter.father_bam_file=' + 
			clean_bams[1], '-P', 'mosaic_filter.mother_bam_file=' + clean_bams[2], '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=trio',
			 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
			 '-P', 'output_dir=' + outdir, '-P', 'syscall_filter.depth=' + pro_cov, '-P', 'depth_filter.min_depth=10'])
	run_mh_trio.wait()
	outfiles.append(outdir + '/final.passed.tsv')
	##lower aaf
	outdir = pro_name + '_mh_trio_results_lowaaf'
	run_mh_trio = subprocess.Popen(['java', '-jar', mosaic_hunter, 'exome', '-P', 'input_file=' + clean_bams[0], '-P', 'mosaic_filter.father_bam_file=' + 
			clean_bams[1], '-P', 'mosaic_filter.mother_bam_file=' + clean_bams[2], '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=trio',
			 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
			 '-P', 'output_dir=' + outdir, '-P', 'syscall_filter.depth=' + pro_cov, '-P', 'min_minor_allele_number=2', '-P', 'min_minor_allele_percentage=1', 
			 '-P', 'base_number_filter.min_minor_allele_number=2', '-P', 'base_number_filter.min_minor_allele_percentage=1'])
	run_mh_trio.wait()
	outfiles.append(outdir + '/final.passed.tsv')
	##no filter
	outdir = pro_name + '_mh_trio_results_nofilter'
	run_mh_trio = subprocess.Popen(['java', '-jar', mosaic_hunter, 'exome', '-P', 'input_file=' + clean_bams[0], '-P', 'mosaic_filter.father_bam_file=' + 
			clean_bams[1], '-P', 'mosaic_filter.mother_bam_file=' + clean_bams[2], '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=trio',
			 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
			 '-P', 'output_dir=' + outdir, '-P', 'syscall_filter.depth=' + pro_cov, '-P', 'min_p_value=0', '-P', ' strand_bias_filter.p_value_cutoff=0', 
			 '-P', 'within_read_position_filter.p_value_cutoff=0'])
	run_mh_trio.wait()
	outfiles.append(outdir + '/final.passed.tsv')
	##all 3
	outdir = pro_name + '_mh_trio_results_all3'
	run_mh_trio = subprocess.Popen(['java', '-jar', mosaic_hunter, 'exome', '-P', 'input_file=' + clean_bams[0], '-P', 'mosaic_filter.father_bam_file=' + 
			clean_bams[1], '-P', 'mosaic_filter.mother_bam_file=' + clean_bams[2], '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=trio',
			 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
			 '-P', 'output_dir=' + outdir, '-P', 'syscall_filter.depth=' + pro_cov, '-P', 'depth_filter.min_depth=10', '-P', 'min_minor_allele_number=2', 
			 '-P', 'min_minor_allele_percentage=1', '-P', 'min_p_value=0', '-P', 'base_number_filter.min_minor_allele_number=2', '-P', 
			 'base_number_filter.min_minor_allele_percentage=1', '-P', 'strand_bias_filter.p_value_cutoff=0', '-P', 'within_read_position_filter.p_value_cutoff=0'])
	run_mh_trio.wait()
	outfiles.append(outdir + '/final.passed.tsv')
	##delete clean bams
	for b in clean_bams:
		os.remove(b)
		os.remove(b + '.bai')
		print 'deleting', b, b + '.bai'
	return outfiles

def reformat_annovar(multianno, annotated):
	with open(multianno, "r") as infh, open(annotated, "w") as outfh:
		line_count = 0
		for line in infh:
			line_count += 1
			line = line.strip('\n').split(delim)
			if line_count == 1:
				header = line[:66] + ['ref/alt', 'aaf', 'prior probability', '\n']
				outfh.write(delim.join(header))
			else:
				location = line[5]
				if location == 'exonic' or location == 'splicing':
					max_maf = line[46]
					if max_maf == '.':
						max_maf = 0
					if float(max_maf) <= 0.01:
						alleles = line[66].split(',')
						aaf = float(alleles[1])/ (int(alleles[1]) + int(alleles[0]))
						print location, alleles, aaf
						line_out = line[:67] + [str(aaf), line[67], '\n']
						outfh.write(delim.join(line_out))

def annotate_mh_output(s_name, files_to_annotate):
	for mh_var_file in files_to_annotate:
		print mh_var_file
		analysis = mh_var_file.split('/')[0].split('_')[-1]
		print analysis
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
		#filter annovar output
		multianno = s_name + '.hg19_multianno.txt'
		annotated = s_name + '.mosaichunter_trio.' + analysis +  '.xls'
		reformat_annovar(multianno, annotated)


def annotate_mutect_output(s_name):
	vcf = name_prefix + '.mutect.vcf.gz'
	##run annovar
	command = [table_annovar] + av_buildver + [vcf] + av_ref_dir + av_protocol + av_operation + av_options_vcf + ['-out', s_name]
	annovar = subprocess.Popen(command)
	annovar.wait()

def run_mosaic_variant_calling(work_dir, ped_name, bam_files, proband_sex, proband_coverage, analysis_type):
	os.chdir(work_dir)
	if analysis_type == 'paired':
		variant_calling_paired_mutect2(bam_files, ped_name)
		annotate_mutect_output(ped_name)
	if analysis_type == 'trio':
		sample_name = bam_files[0].split('.')[0]
		files_to_annotate = variant_calling_trio_mosaichunter(sample_name, bam_files, proband_sex, proband_coverage)
		annotate_mh_output(sample_name, files_to_annotate)

##run methods


##paired analysis with mutect
# working_dir = '/data/atimms/dobyns_mosaic_test_0317'
# run_mosaic_variant_calling(working_dir, 'LR16-173_paired', ['LR16-173_avm_lesion.bwa_gatk.bam', 'LR16-173_saliva.bwa_gatk.bam'], '', '', 'paired')
# run_mosaic_variant_calling(working_dir, 'LR13-356_paired', ['LR13-356_avm_lesion.bwa_gatk.bam', 'LR13-356_saliva.bwa_gatk.bam'], '', '', 'paired')
working_dir = '/data/atimms/novartis_fcd_exomes_0317'
# run_mosaic_variant_calling(working_dir, 'LR12-243_paired', ['LR12-243_brain.bwa_gatk.bam', 'LR12-243_gl.bwa_gatk.bam'], 'M', '297', 'paired')
# run_mosaic_variant_calling(working_dir, 'LR12-245_paired', ['LR12-245_brain.bwa_gatk.bam', 'LR12-245_gl.bwa_gatk.bam'], 'F', '301', 'paired')
# run_mosaic_variant_calling(working_dir, 'LR12-250_paired', ['LR12-250_brain.bwa_gatk.bam', 'LR12-250_gl.bwa_gatk.bam'], 'F', '291', 'paired')
# run_mosaic_variant_calling(working_dir, 'LR12-255_paired', ['LR12-255_brain.bwa_gatk.bam', 'LR12-255_gl.bwa_gatk.bam'], 'M', '233', 'paired')
# run_mosaic_variant_calling(working_dir, 'LR12-259_paired', ['LR12-259_brain.bwa_gatk.bam', 'LR12-259_gl.bwa_gatk.bam'], 'F', '272', 'paired')
run_mosaic_variant_calling(working_dir, 'LR12-260_paired', ['LR12-260_brain.bwa_gatk.bam', 'LR12-260_gl.bwa_gatk.bam'], 'F', '344', 'paired')
# run_mosaic_variant_calling(working_dir, 'LR12-269_paired', ['LR12-269_brain.bwa_gatk.bam', 'LR12-269_gl.bwa_gatk.bam'], 'M', '203', 'paired')
##mom not realated so don't analyze this way.
# run_mosaic_variant_calling(working_dir, 'LR13-354', ['LR13-354_brain.bwa_gatk.bam', 'LR13-354f.bwa_gatk.bam', 'LR13-354m.bwa_gatk.bam'], 'F', '', 'trio')
##trio analysis with mosaic humnter
##what the fields are
##run_mosaic_variant_calling(working_dir, 'name', ['pro_bam', 'dad_bam', 'mom_bam'], 'sex', 'pro_cov', 'analysis')
##empty one
##run_mosaic_variant_calling(working_dir, '', ['', '', ''], '', '', 'trio')

##jimmys files
# working_dir = '/data/atimms/dobyns_mosaic_test_0317'
# run_mosaic_variant_calling(working_dir, 'LR16-173', ['LR16-173_avm_lesion.bwa_gatk.bam', 'LR16-173f.bwa_gatk.bam', 'LR16-173m.bwa_gatk.bam'], 'F', '320', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR13-356', ['LR13-356_avm_lesion.bwa_gatk.bam', 'LR13-356f.bwa_gatk.bam', 'LR13-356m.bwa_gatk.bam'], 'M', '385', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR16-030', ['LR16-030.bwa_gatk.bam', 'LR16-030f.bwa_gatk.bam', 'LR16-030m.bwa_gatk.bam'], 'M', '360', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR16-031', ['LR16-031.bwa_gatk.bam', 'LR16-031f.bwa_gatk.bam', 'LR16-031m.bwa_gatk.bam'], 'M', '980', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR16-172', ['LR16-172.bwa_gatk.bam', 'LR16-172f.bwa_gatk.bam', 'LR16-172m.bwa_gatk.bam'], 'M', '355', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR14-252', ['LR14-252_frozen_scalp.bwa_gatk.bam', 'LR14-252f.bwa_gatk.bam', 'LR14-252m.bwa_gatk.bam'], 'F', '475', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR14-252', ['LR14-252_skin.bwa_gatk.bam', 'LR14-252f.bwa_gatk.bam', 'LR14-252m.bwa_gatk.bam'], 'F', '190', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR16-065', ['LR16-065_eye.bwa_gatk.bam', 'LR16-065f.bwa_gatk.bam', 'LR16-065m.bwa_gatk.bam'], 'M', '490', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR16-065', ['LR16-065_saliva.bwa_gatk.bam', 'LR16-065f.bwa_gatk.bam', 'LR16-065m.bwa_gatk.bam'], 'M', '160', 'trio')
##kims peds
# working_dir = '/data/atimms/dobyns_mosaic_test_0317'
# run_mosaic_variant_calling(working_dir, 'LR04-350', ['LR04-350.bwa_gatk.bam', 'LR04-350f.bwa_gatk.bam', 'LR04-350m.bwa_gatk.bam'], 'M', '72', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR08-108', ['LR08-108.bwa_gatk.bam', 'LR08-108f.bwa_gatk.bam', 'LR08-108m.bwa_gatk.bam'], 'M', '109', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR08-337', ['LR08-337.bwa_gatk.bam', 'LR08-337f.bwa_gatk.bam', 'LR08-337m.bwa_gatk.bam'], 'M', '81', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR09-203', ['LR09-203.bwa_gatk.bam', 'LR09-203f.bwa_gatk.bam', 'LR09-203m.bwa_gatk.bam'], 'M', '69', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR09-387', ['LR09-387.bwa_gatk.bam', 'LR09-387f.bwa_gatk.bam', 'LR09-387m.bwa_gatk.bam'], 'M', '191', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR14-048', ['LR14-048.bwa_gatk.bam', 'LR14-048f.bwa_gatk.bam', 'LR14-048m.bwa_gatk.bam'], 'M', '132', 'trio')
##ghayda's data --- need average cov
# working_dir = '/data/atimms/novartis_fcd_exomes_0317'
# run_mosaic_variant_calling(working_dir, 'LR12-243', ['LR12-243_brain.bwa_gatk.bam', 'LR12-243f.bwa_gatk.bam', 'LR12-243m.bwa_gatk.bam'], 'M', '297', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR12-245', ['LR12-245_brain.bwa_gatk.bam', 'LR12-245f.bwa_gatk.bam', 'LR12-245m.bwa_gatk.bam'], 'F', '301', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR12-250', ['LR12-250_brain.bwa_gatk.bam', 'LR12-250f.bwa_gatk.bam', 'LR12-250m.bwa_gatk.bam'], 'F', '291', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR12-255', ['LR12-255_brain.bwa_gatk.bam', 'LR12-255f.bwa_gatk.bam', 'LR12-255m.bwa_gatk.bam'], 'M', '233', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR12-259', ['LR12-259_brain.bwa_gatk.bam', 'LR12-259f.bwa_gatk.bam', 'LR12-259m.bwa_gatk.bam'], 'F', '272', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR12-260', ['LR12-260_brain.bwa_gatk.bam', 'LR12-260f.bwa_gatk.bam', 'LR12-260m.bwa_gatk.bam'], 'F', '344', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR12-269', ['LR12-269_brain.bwa_gatk.bam', 'LR12-269f.bwa_gatk.bam', 'LR12-269m.bwa_gatk.bam'], 'M', '203', 'trio')
##mom not realated so don't analyze this way.
# run_mosaic_variant_calling(working_dir, 'LR13-354', ['LR13-354_brain.bwa_gatk.bam', 'LR13-354f.bwa_gatk.bam', 'LR13-354m.bwa_gatk.bam'], 'F', '', 'trio')





