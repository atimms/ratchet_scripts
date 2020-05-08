#!/usr/bin/python
import sys
import subprocess

'''
To do:

Dependancies:
must have samtools, bwa, tabix and this file in PATH

update path to mipgen programs

update parameters below 
- don't touch anything below the line

make file with list of genenames i.e. example.txt

move to directory with file

if designing from gene list in text file
type:  migen_pipeline.py <file with list of genenames> 'genelist'

if designing from bed file
type:  migen_pipeline.py <bedfile> 'bed'

if designing from bed file, for non-unique regions
type:  migen_pipeline.py <bedfile> 'bed_not_unique'

if checking list of genes
type:  migen_pipeline.py <file with list of genenames> 'gene_check'

all results files will be in that directory and will start with the 
same name as the original file i.e. example
'''

##parameters to update

##programs
extract_coding_exons = '/home/atimms/programs/MIPGEN-master/tools/extract_coding_gene_exons.sh'
mipgen = '/home/atimms/programs/MIPGEN-master/mipgen'

##files
refGene = '/data/atimms/references/refGene_noStrangeChr.txt'
#fa file must be indexed
hg19_fasta = '/data/atimms/references/human_g1k_v37.fasta'
#snpfile must be compressed with bgzip and indexed with tabix i.e. >bgzip <file> >tabix -p vcf <file>
snpfile = '/data/atimms/references/common_all.vcf.gz'

##mipgen instructions
min_size = 162
max_size = 162
#tag_sizes = '4,4'
tag_sizes = '5,0'
##can be 'std' or 'ff10'
design_parameters = 'ff10'
##minimum length of extension and ligation arms (default is 16 and 18) - recommended >=20 (especially for low gc)
ext_min_length = 16
lig_min_length = 18


########################################################################################

##To do............
#auto check files are correct and indexed
#add more parameters for mipgen
#summary data for mipgen section


message = '''
###########################################
# pipeline for using mipgen to make MIPs: #
###########################################

Usage:

Dependancies:
must have samtools, bwa, tabix and this file in PATH

update path to mipgen programs

update parameters in migen_pipeline.py above the line 

make file with list of genenames i.e. example.txt

move to directory with file

if designing from gene list in text file
type:  migen_pipeline.py <file with list of genenames> 'genelist'

if designing from bed file
type:  migen_pipeline.py <bedfile> 'bed'

if designing from bed file, for non-unique regions
type:  migen_pipeline.py <bedfile> 'bed_not_unique'

if checking list of genes
type:  migen_pipeline.py <file with list of genenames> 'gene_check'

all results files will be in that directory and will start with the 
same name as the original file i.e. example


##################
# summary of run #
##################

'''

##parameters
genelist = sys.argv[1]
analysis = sys.argv[2]
delim = '\t'
bed_suffix = '_coding_exons.bed'

def make_refGene_into_list_of_genes(refGene):
	gene_list = []
	with open(refGene, "r") as refgene:
		for line in refgene:
			line = line.strip('\n').split(delim)
			gene = line[12]
			gene_list.append(gene)
	return gene_list

def are_genenames_in_bed(gene_list):
	refgene_genes = make_refGene_into_list_of_genes(refGene)
	not_in_refgene = []
	with open(gene_list, "U") as file:
		for gene in file:
			gene = gene.rstrip()
			if gene not in refgene_genes:
				not_in_refgene.append(gene)
	return not_in_refgene

def convert_line_end_to_unix(gene_list):
	outfile = gene_list.split('.')[0] + '.temp'
	with open(gene_list, "U") as gfile, open(outfile, "w") as final:
		for g in gfile:
			final.write(g.strip() + '\n')

def extract_coding_exons_from_genelist(gene_list):
	outfile = gene_list.split('.')[0] + bed_suffix
	command = extract_coding_exons + ' ' + gene_list + ' ' + refGene + ' > ' + outfile 
	subprocess.call(command, shell=True)

def run_mipgen(gene_list):
	name = gene_list.split('.')[0]
	bed = name + bed_suffix
	if design_parameters == 'std':
		command = mipgen + ' -regions_to_scan ' + bed + ' -project_name ' + name + ' -min_capture_size ' + str(min_size) + ' -max_capture_size ' + str(max_size) + ' -bwa_genome_index ' + hg19_fasta + ' -snp_file ' + snpfile + ' -tag_sizes ' + tag_sizes + ' -ext_min_length ' + str(ext_min_length) + ' -lig_min_length ' + str(lig_min_length)
	elif design_parameters == 'ff10':
		command = mipgen + ' -regions_to_scan ' + bed + ' -project_name ' + name + ' -min_capture_size ' + str(min_size) + ' -max_capture_size ' + str(max_size) + ' -bwa_genome_index ' + hg19_fasta + ' -snp_file ' + snpfile + ' -tag_sizes ' + tag_sizes + ' -score_method mixed' + ' -feature_flank 10' + ' -ext_min_length ' + str(ext_min_length) + ' -lig_min_length ' + str(lig_min_length)
	else:
		print 'design parameter choice not recognised, please change'
	subprocess.call(command, shell=True)

def delete_intermediate_files():
	command = 'rm -f *.sai *.fq *sam *.all_mips.txt *.temp *.fa'
	subprocess.call(command, shell=True)

def call_and_summarize_results(genes):
	convert_line_end_to_unix(genes)
	converted_genes = genes.split('.')[0] + '.temp'
	extract_coding_exons_from_genelist(converted_genes)
	missing_list = are_genenames_in_bed(converted_genes)
	missing_genes = ' '.join(missing_list)
	run_mipgen(converted_genes)
	delete_intermediate_files()
	print message
	if len(missing_list) >= 1:
		print "the folowing genes weren't found in the reference file:", missing_genes + '\n'
	else:
		print "all genes in genelist were found in reference file" + '\n'
	exons = genes.split('.')[0] + bed_suffix
	print "the file '%s' was created, which contains all coding exons from the genes in file '%s'"% (exons, genes) + '\n'

def run_mipgen_not_unique(bed):
	name = bed.split('.')[0]
	if design_parameters == 'std':
		command = mipgen + ' -regions_to_scan ' + bed + ' -project_name ' + name + ' -min_capture_size ' + str(min_size) + ' -max_capture_size ' + str(max_size) + ' -bwa_genome_index ' + hg19_fasta + ' -snp_file ' + snpfile + ' -tag_sizes ' + tag_sizes + ' -check_copy_number off'
	elif design_parameters == 'ff10':
		command = mipgen + ' -regions_to_scan ' + bed + ' -project_name ' + name + ' -min_capture_size ' + str(min_size) + ' -max_capture_size ' + str(max_size) + ' -bwa_genome_index ' + hg19_fasta + ' -snp_file ' + snpfile + ' -tag_sizes ' + tag_sizes + ' -score_method mixed' + ' -feature_flank 10' + ' -check_copy_number off'
	else:
		print 'design parameter choice not recognised, please change'
	subprocess.call(command, shell=True)
	delete_intermediate_files()

def run_mipgen_bed(bed):
	name = bed.split('.')[0]
	if design_parameters == 'std':
		command = mipgen + ' -regions_to_scan ' + bed + ' -project_name ' + name + ' -min_capture_size ' + str(min_size) + ' -max_capture_size ' + str(max_size) + ' -bwa_genome_index ' + hg19_fasta + ' -snp_file ' + snpfile + ' -tag_sizes ' + tag_sizes + ' -ext_min_length ' + str(ext_min_length) + ' -lig_min_length ' + str(lig_min_length)
	elif design_parameters == 'ff10':
		command = mipgen + ' -regions_to_scan ' + bed + ' -project_name ' + name + ' -min_capture_size ' + str(min_size) + ' -max_capture_size ' + str(max_size) + ' -bwa_genome_index ' + hg19_fasta + ' -snp_file ' + snpfile + ' -tag_sizes ' + tag_sizes + ' -score_method mixed' + ' -feature_flank 10' + ' -ext_min_length ' + str(ext_min_length) + ' -lig_min_length ' + str(lig_min_length)
	else:
		print 'design parameter choice not recognised, please change'
	subprocess.call(command, shell=True)
	delete_intermediate_files()

def check_if_genes_in_refgene(genes):
	convert_line_end_to_unix(genes)
	converted_genes = genes.split('.')[0] + '.temp'
	missing_list = are_genenames_in_bed(converted_genes)
	missing_genes = ' '.join(missing_list)
	print message
	if len(missing_list) >= 1:
		print "the folowing genes weren't found in the reference file:", missing_genes + '\n'
	else:
		print "all genes in genelist were found in reference file" + '\n'

##call methods depending on type of input
if analysis == 'genelist':	
	call_and_summarize_results(genelist)
elif analysis == 'bed_not_unique':
	run_mipgen_not_unique(genelist)
elif analysis == 'bed':
	run_mipgen_bed(genelist)

elif analysis == 'gene_check':
	check_if_genes_in_refgene(genelist)
