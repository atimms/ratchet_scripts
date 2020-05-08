#!/usr/bin/env python
import os
import subprocess
import glob
import shutil

'''
module load biobuilds/2017.05
'''


##set input variables and parameters
delim = '\t'

##programs and ref files
vt = '/data/atimms/gemini/vt/vt'
gemini = '/data/atimms/gemini/bin/gemini'
fasta = '/data/atimms/references/human_g1k_v37.fasta'
# fasta = '/data/atimms/references/grch37.fa'
snpeff_jar = '/data/atimms/gemini/snpEff/snpEff.jar'
bcftools_12 = '/home/atimms/programs/bcftools-1.2/bcftools'

##annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/data/atimms/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,dbnsfp35a,mcap13,revel,gnomad211_exome,gnomad211_genome,cosmic88_coding,cosmic88_noncoding,rmsk,genomicSuperDups']
av_operation = ['-operation', 'g,f,f,f,f,f,f,f,r,r']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput', '-arg', '-splicing 10 ,,,,,,,,,']
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'

##methods
def annotate_with_annovar(vcf, out_prefix):
	##run_table_annovar(avinput, av_prefix):
	command = [table_annovar] + av_buildver + [vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', out_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()


def format_clinvar_file(infile, outfile):
	clinvar_labels = ['CLNDISDB', 'CLNDN', 'CLNHGVS', 'CLNREVSTAT', 'CLNSIG', 'CLNVC', 'CLNVCSO']
	dropped_line_count, dropped_clinvar_count = 0,0
	with open(infile, 'r') as in_fh, open(outfile, 'w') as out_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			line = line.rstrip().split(delim)
			if line_count == 1:
				header = line[:120] + clinvar_labels
				out_fh.write(delim.join(header) + '\n')
			else:
				# print(line)
				# print(len(line))
				##remove few rows with data missing
				if len(line) == 131:
					line_out = line[:120]
					clinvar_data = line[130].split(';')
					clinvar_labels_plus = [i + '=' for i in clinvar_labels]
					c_counts = 0
					##make sure we have all clivar data
					for c in clinvar_data:
						if c.startswith(tuple(clinvar_labels_plus)):
							c_counts += 1
					# print(len(c_counts), len(clinvar_labels))
					# print(clinvar_data)
					if c_counts == len(clinvar_labels):
						for cd in clinvar_data:
							if cd.startswith(tuple(clinvar_labels_plus)):
								line_out.append(cd.split('=')[1])
						out_fh.write(delim.join(line_out) + '\n')
					else:
						dropped_clinvar_count += 1

				else:
					dropped_line_count += 1
					# print(line)
					# print(len(line))
	print(dropped_line_count, dropped_clinvar_count)

def filter_vars_gene_clinsig(infile, gene_infiles, clinsigs_wanted):
	for gene_infile in gene_infiles:
		gene_list = []
		with open(gene_infile, 'r') as gi_fh:
			for line in gi_fh:
				gene = line.rstrip()
				gene_list.append(gene)
		print(gene_infile, len(gene_list))
		outfile = infile.split('.')[0] + '.' + gene_infile.split('.')[0] + '.annotated.xls'
		with open(infile, 'r') as in_fh, open(outfile, 'w') as out_fh:
			line_count = 0
			for line in in_fh:
				line_count += 1
				if line_count == 1:
					out_fh.write(line)
				else:
					line = line.split(delim)
					gene = line[6]
					# print(line)
					clinsig = line[124]
					if gene in gene_list and clinsig in clinsigs_wanted:
						out_fh.write(delim.join(line))


##run all other methods
def annotate_clinvar_vcf(working_directory, out_prefix, in_vcf, gene_files):
	os.chdir(working_directory)
	##annoate
	# annotate_with_annovar(in_vcf, out_prefix)
	##change header and format
	multianno = out_prefix + '.hg19_multianno.txt'
	annotated =  out_prefix + '.annotated.xls'
	# format_clinvar_file(multianno, annotated)
	##filter by status/gene lists
	##only keep these labels
	clinsigs_to_keep = ['Benign', 'Likely_benign', 'Benign/Likely_benign', 'Pathogenic', 
		'Likely_pathogenic', 'Pathogenic/Likely_pathogenic', 'Uncertain_significance']
	filter_vars_gene_clinsig(annotated, gene_files, clinsigs_to_keep)

##run methods
clinvar_vcf = 'clinvar_20190731.vcf.gz'
project = clinvar_vcf.split('.')[0]
w_dir = '/data/atimms/gm_exome_project_0719'
gene_lists = ['pirozzi_genes.txt', 'dbdb_genes.txt']

annotate_clinvar_vcf(w_dir, project, clinvar_vcf, gene_lists)


