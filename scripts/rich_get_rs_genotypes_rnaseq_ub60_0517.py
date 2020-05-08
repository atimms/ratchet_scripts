#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import filtering_annotated

##parameters
delim = '\t'
genome_name = 'hg19'
fasta = '/data/atimms/references/igenomes/' + genome_name + '/genome.fa'

#annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/data/atimms/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'avsnp147']
av_operation = ['-operation', 'f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput']
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
dbsnp_file = '/data/atimms/references/annovar/hg19/hg19_avsnp147.txt'


##methods
def table_annovar_vcf(vcf, project_prefix):
	out_prefix = project_prefix
	command = [table_annovar] + av_buildver + [vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', out_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()

def annotate_vcf_file(vcf, project_prefix, samples):
	table_annovar_vcf(vcf, project_prefix)


def get_genotypes_from_vcf(vcf, rs_list, outfile):
	with open(vcf, "r") as vcf_fh, open(outfile, "w") as out_fh:
		line_count = 0
		for line in vcf_fh:
			line_count += 1
			if line[:2] != '##':
				
				if line[0] == "#":
					line = line.rstrip().split(delim)
					samples = line[9:]
					head_out = ['rs_number', 'chr', 'pos', 'ref', 'alt'] + samples + ['\n']
					out_fh.write(delim.join(head_out))
				else:
					line = line.rstrip().split(delim)
					# info = line[7].split(';')
					info = line[7]
					chrom_etc = line[:2] + line[3:5]
					genotypes = line[9:]
					for rs in rs_list:
						if rs + ';' in info:
							print rs, chrom_etc
							line_out = [rs] + chrom_etc 
							for geno in genotypes:
								geno = geno.split(':')[0]
								line_out.append(geno)
							out_fh.write(delim.join(line_out + ['\n']))

def make_bed_from_rs_numbers(rs_list, dbsnp, bed_file):
	with open(dbsnp, "r") as in_fh, open(bed_file, "w") as out_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			# print line
			rs = line[5]
			# print rs, rs_numbers
			if rs in rs_list:
				chrom = 'chr' + line[0]
				start = int(line[1]) - 4
				rest = line[2:] + ['\n']
				print rs, line, chrom, start, rest
				out_fh.write(delim.join([chrom, str(start)] + rest))

def make_list_of_bams(sample_dict, bam_suffix, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for sample in sample_dict:
			bam = sample + bam_suffix
			outfh.write(bam + '\n')

def  index_bam_file_for_freebayes(samples, bam_suffix):
	bam_files = []
	for sample in samples:
		bam = sample + bam_suffix
		bam_files.append(bam)
		##check bai file exists and make if it isn't there
		if os.path.isfile(bam + '.bai'):
			print 'bam %s alreaded indexed'%bam
		else:
			print 'indexing bam file:', bam
			st_index = subprocess.Popen(['samtools', 'index', bam])
			st_index.wait()

def genotype_vars_with_freebayes(bam_file, final_vcf, bedfile):
	vcf_temp1 = 'temp_fb1.vcf'
	vcf_temp2 = 'temp_fb2.vcf'
	with open(vcf_temp1, 'w') as vcf_fh:
		freebayes = subprocess.Popen(['freebayes', '-f', fasta, '-L', bam_file, '-t', bedfile, '--haplotype-length', '0', '--min-alternate-count', '1', '--min-alternate-fraction', '0', '--pooled-continuous', '--report-monomorphic'], stdout=vcf_fh)
		freebayes.wait()
	bgzip = subprocess.Popen(['bgzip', vcf_temp1])
	bgzip.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	os.remove( vcf_temp1 + '.gz')


def make_genotype_file(bedfile, vcffile, outfile):
	##make list from bed file
	bed_var_dict = {}
	with open(bedfile, "r") as inb_fh, open(outfile, "w") as out_fh:
		for line_bed in inb_fh:
			line_bed = line_bed.rstrip().split(delim)
			bed_var = '_'.join([line_bed[0]] + line_bed[2:5])
			rs_number = line_bed[5]
			print bed_var
			bed_var_dict[bed_var] = rs_number
	with open(vcffile, "r") as vcf_fh, open(outfile, "w") as out_fh:
		for line in vcf_fh:
			if line[:2] != '##':
				if line[0] == "#":
					line = line.rstrip().split(delim)
					samples = line[9:]
					head_out = ['rs_number', 'chr', 'pos', 'ref', 'alt'] + samples + ['\n']
					out_fh.write(delim.join(head_out))
				else:
					line = line.rstrip().split(delim)
					# info = line[7].split(';')
					info = line[7]
					chrom_etc = line[:2] + line[3:5]
					genotypes = line[9:]
					vcf_var = '_'.join(line[:2] + line[3:5])
					print vcf_var
					if vcf_var in bed_var_dict:
						rs = bed_var_dict[vcf_var]
						print 'found it', vcf_var,rs
						line_out = [rs] + chrom_etc 
						for geno in genotypes:
							geno = geno.split(':')[0]
							if geno == '0/0':
								geno = 'hom_wt'
							elif geno == '0/1' or geno == '1/0':
								geno = 'het'
							elif geno == '1/1':
								geno = 'hom_mut'
							elif geno == '.':
								geno = 'not_covered'
							line_out.append(geno)
						out_fh.write(delim.join(line_out + ['\n']))


##bri batch
'''
working_dir = '/data/atimms/rich_rnaseq/asthma_rnaseq_0217'
os.chdir(working_dir)
project_name = 'asthma_bri_rs_numbers'
bri_samples = ['N353', 'T112', 'T147', 'T148', 'T149', 'T150', 'T151', 'T152', 'T154', 'T155', 'T241', 'T352', 'T354', 'T355', 
		'T356', 'T359', 'T361', 'T362', 'T363', 'T365', 'T366', 'T368', 'T370', 'T372', 'T373']
rs_numbers = ['rs1990760', 'rs112363590', 'rs3205166', 'rs55717844', 'rs10813831', 'rs112437427', 'rs1131665', 'rs1061505', 
		'rs1061502', 'rs113323265', 'rs112815033', 'rs7972447', 'rs1138106178', 'rs17857295', 'rs7269320', 'rs113404862', 
		'rs3775290', 'rs6830345', 'rs111488413', 'rs113258886', 'rs113275119', 'rs112666655', 'rs864058', 'rs5741881', 
		'rs7486100', 'rs3746660', 'rs3775291', 'rs3853839', 'rs2292151', 'rs3747414', 'rs2274863']

##deleted multianno file - so re make with just dbsnp
##these don't work as not all info in vcf file
# annotate_vcf_file('asthma_bri_0217.intersected_vcfs/0002.vcf', project_name, bri_samples)
# get_genotypes_from_vcf(project_name + '.hg19_multianno.vcf', rs_numbers, project_name + '.genotypes.xls')
##so get rs info and genotype using samtools
make_bed_from_rs_numbers(rs_numbers, dbsnp_file, project_name+'.bed')
make_list_of_bams(bri_samples, '.star.gatk.bam', 'bam.list')
# index_bam_file_for_freebayes(bri_samples, '.star.gatk.bam')
genotype_vars_with_freebayes('bam.list', project_name + '.vcf', project_name+'.bed')
make_genotype_file(project_name+'.bed', project_name + '.vcf', project_name + '.genotypes.xls')
'''

##amgen_batch
working_dir = '/data/atimms/rich_rnaseq/amgen_asthma_rnaseq_0317'
os.chdir(working_dir)
project_name = 'asthma_amgen_rs_numbers'
amgen_samples = ['T130', 'T219', 'T132', 'T133', 'T134', 'T136', 'T137', 'N135', 'T303', 'T321', 'T320', 'T323', 'T325', 'T324', 'T104', 'T131', 'T108', 'T222', 'T223', 'T329', 'T128', 'T145', 'T122', 'T127', 'T126', 'T125', 'T124', 'T313', 'T143', 'T144', 'T140', 'T312', 'T216', 'T314', 'T319', 'T332', 'T117', 'T114', 'T115', 'T111', 'T119']
rs_numbers = ['rs1990760', 'rs112363590', 'rs3205166', 'rs55717844', 'rs10813831', 'rs112437427', 'rs1131665', 'rs1061505', 
		'rs1061502', 'rs113323265', 'rs112815033', 'rs7972447', 'rs1138106178', 'rs17857295', 'rs7269320', 'rs113404862', 
		'rs3775290', 'rs6830345', 'rs111488413', 'rs113258886', 'rs113275119', 'rs112666655', 'rs864058', 'rs5741881', 
		'rs7486100', 'rs3746660', 'rs3775291', 'rs3853839', 'rs2292151', 'rs3747414', 'rs2274863']

##deleted multianno file - so re make with just dbsnp
##these don't work as not all info in vcf file
# annotate_vcf_file('asthma_bri_0217.intersected_vcfs/0002.vcf', project_name, bri_samples)
# get_genotypes_from_vcf(project_name + '.hg19_multianno.vcf', rs_numbers, project_name + '.genotypes.xls')
##so get rs info and genotype using samtools
# make_bed_from_rs_numbers(rs_numbers, dbsnp_file, project_name+'.bed')
# make_list_of_bams(amgen_samples, '.star.gatk.bam', 'bam.list')
# index_bam_file_for_freebayes(amgen_samples, '.star.gatk.bam')
# genotype_vars_with_freebayes('bam.list', project_name + '.vcf', project_name+'.bed')
make_genotype_file(project_name+'.bed', project_name + '.vcf', project_name + '.genotypes.xls')
