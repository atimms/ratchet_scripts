#!/usr/bin/env python
import sys
import subprocess
import os
import glob


##parameters
delim = '\t'

##setup working directory where results will be
working_dir = '/data/atimms/serge_virus_id_0217/'
clark_prog_dir = '/home/atimms/programs/CLARKSCV1.2.3/'



##programs etc
clark_db_dir = '/data/atimms/clark_databases'
read1_files = working_dir + 'read1_files.txt'
read2_files = working_dir + 'read2_files.txt'
results_files = working_dir + 'results_files.txt'
sample_names = ['DNA_Cortex', 'DNA_Stem', 'RNA_Cortex', 'RNA_Stem']
seqtk = '/home/atimms/programs/seqtk/seqtk'
bwa = '/home/atimms/programs/bwa-0.7.12/bwa'

##methods
def remove_human_reads(clark_results_file):
	print 'for file:', clark_results_file
	non_human_read_ids = clark_results_file.split('.')[0] + '.temp.read_ids.txt'
	print 'writing non human read ids to ', non_human_read_ids
	with open (clark_results_file, 'r') as in_fh, open (non_human_read_ids, 'w') as out_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			if line_count > 1:
				line = line.rstrip().split(',')
				t_id = line[2]
				##if not human
				if t_id != '9606':
					read_id = line[0]
					out_fh.write(read_id + '\n')
	##extract reads using seqtk
	##seqtk subseq in.fq name.lst > out.fq
	r1_fq = clark_results_file.split('.')[0] + '.R1.fq'
	r2_fq = clark_results_file.split('.')[0] + '.R2.fq'
	nonh_r1_fq = clark_results_file.split('.')[0] + '.non_human.R1.fq'
	nonh_r2_fq = clark_results_file.split('.')[0] + '.non_human.R2.fq'
	with open(nonh_r1_fq, 'w') as r1_fh:
		run_seqtk_r1  = subprocess.Popen([seqtk, 'subseq', r1_fq, non_human_read_ids], stdout=r1_fh)
		run_seqtk_r1.wait()
	with open(nonh_r2_fq, 'w') as r2_fh:
		run_seqtk_r2  = subprocess.Popen([seqtk, 'subseq', r2_fq, non_human_read_ids], stdout=r2_fh)
		run_seqtk_r2.wait()
	return nonh_r1_fq, nonh_r2_fq

def map_against_genomes(fq_files, fa_file):
	sample_name = fq_files[0].split('.')[0]
	r1_bam = sample_name + '.non_human.clark_default.7genomes.r1.bam'
	r2_bam = sample_name + '.non_human.clark_default.7genomes.r2.bam'
	with open(r1_bam, 'w') as bam_fh:
		bwa_r1 = subprocess.Popen(['bwa', 'mem', '-t', '10', fa_file, fq_files[0]], stdout=subprocess.PIPE)
		sam_bam_r1 = subprocess.Popen(['samtools', 'view', '-u', '-F', '4', '-'], stdin=bwa_r1.stdout, stdout=subprocess.PIPE)
		st_sort_r1 = subprocess.Popen(['samtools', 'sort','-O', 'bam', '-T', sample_name, '-'], stdin=sam_bam_r1.stdout, stdout=bam_fh)
		st_sort_r1.wait()
		st_r1_index = subprocess.Popen(['samtools', 'index', r1_bam])
		st_r1_index.wait()
	with open(r2_bam, 'w') as bam_fh:
		bwa_r2 = subprocess.Popen(['bwa', 'mem', '-t', '10', fa_file, fq_files[1]], stdout=subprocess.PIPE)
		sam_bam_r2 = subprocess.Popen(['samtools', 'view', '-u', '-F', '4', '-'], stdin=bwa_r2.stdout, stdout=subprocess.PIPE)
		st_sort_r2 = subprocess.Popen(['samtools', 'sort','-O', 'bam', '-T', sample_name, '-'], stdin=sam_bam_r2.stdout, stdout=bam_fh)
		st_sort_r2.wait()
		st_r2_index = subprocess.Popen(['samtools', 'index', r2_bam])
		st_r2_index.wait()



def get_read_ids_from_clark_results(s_dict, clark_results_file):
	print 'for file:', clark_results_file
	##get list of all tax ids we want
	tax_ids = []
	for s in s_dict:
		tax_id = s_dict[s][1]
		tax_ids.append(tax_id)
	print "we have %s tax ids we're interested in"%len(tax_ids)
	print "which are:", tax_ids
	##make dict of read ids
	read_id_dict = {}
	with open (clark_results_file, 'r') as in_fh:
		for line in in_fh:
			line = line.rstrip().split(',')
			t_id = line[2]
			if t_id in tax_ids:
				read_id = line[0]
				if t_id in read_id_dict:
					read_id_dict[t_id].append(read_id)
				else:
					read_id_dict[t_id] = [read_id]
	##get infor from read dict
	for t in read_id_dict:
		print t, len(read_id_dict[t])
		with open (read_file, 'w') as _fh:
			pass



##run methods
##make databases the first time set_targets is run
##human, virus and bacteria
# os.chdir(clark_prog_dir)
# set_target = subprocess.Popen(['set_targets.sh', clark_db_dir, 'bacteria', 'viruses', 'human'])
# set_target.wait()

##run against one sample -- test threads i.e. -n and s version i.e. --spaced and full mode i.e. -m 0 and kmer length i.e. -k 20
# run_clark = subprocess.Popen(['classify_metagenome.sh', '-P', working_dir + 'test_100k.R1.fq.gz', working_dir + 'test_100k.R2.fq.gz', '-R', working_dir + 'test_100k.results.csv', '--gzipped'])
# run_clark.wait()
##get results for the one sample
# with open(working_dir + 'test_100k.summary.csv', "w") as outfh:
# 	est_ab = subprocess.Popen(['./estimate_abundance.sh', '-F', working_dir + 'test_100k.results.csv.csv', '-D', clark_db_dir], stdout=outfh)
# 	est_ab.wait()

##run against all samples, on at a time -- test threads i.e. -n and s version i.e. --spaced and full mode i.e. -m 0 and kmer length i.e. -k 20
##have to not use gzipped files as temp space is too small

##default
'''
for sample_name in sample_names:
	os.chdir(clark_prog_dir)
	read1 = working_dir + sample_name + '.R1.fq'
	read2 = working_dir + sample_name + '.R2.fq'
	results_file = working_dir + sample_name + '.default_results'
	summary_file = working_dir + sample_name + '.default_summary.csv'
	##set targey -
	set_target = subprocess.Popen(['set_targets.sh', clark_db_dir, 'bacteria', 'viruses', 'human'])
	set_target.wait()
	##align sequences
	run_multi_clark = subprocess.Popen(['classify_metagenome.sh', '-P', read1, read2, '-R', results_file])
	run_multi_clark.wait()
	##get summary of those results
	with open(summary_file, "w") as outfh:
		est_ab = subprocess.Popen(['./estimate_abundance.sh', '-F', results_file + '.csv', '-D', clark_db_dir], stdout=outfh)
		est_ab.wait()
'''
##spaced mode
'''
for sample_name in sample_names:
	os.chdir(clark_prog_dir)
	read1 = working_dir + sample_name + '.R1.fq'
	read2 = working_dir + sample_name + '.R2.fq'
	results_file = working_dir + sample_name + '.spaced_results'
	summary_file = working_dir + sample_name + '.spaced_summary.csv'
	##set target -
	set_target = subprocess.Popen(['set_targets.sh', clark_db_dir, 'bacteria', 'viruses', 'human'])
	set_target.wait()
	##build set of spaced kmers
	build_spaced = subprocess.Popen(['./buildSpacedDB.sh'])
	build_spaced.wait()
	##align sequences
	run_multi_clark = subprocess.Popen(['classify_metagenome.sh', '-P', read1, read2, '--spaced', '-R', results_file])
	run_multi_clark.wait()
	##get summary of those results
	with open(summary_file, "w") as outfh:
		est_ab = subprocess.Popen(['./estimate_abundance.sh', '-F', results_file + '.csv', '-D', clark_db_dir], stdout=outfh)
		est_ab.wait()
'''

##spaced mode and full
'''
for sample_name in sample_names:
	os.chdir(clark_prog_dir)
	read1 = working_dir + sample_name + '.R1.fq'
	read2 = working_dir + sample_name + '.R2.fq'
	results_file = working_dir + sample_name + '.spaced_full_results'
	summary_file = working_dir + sample_name + '.spaced_full_summary.csv'
	##set target -
	set_target = subprocess.Popen(['set_targets.sh', clark_db_dir, 'bacteria', 'viruses', 'human'])
	set_target.wait()
	##build set of spaced kmers
	build_spaced = subprocess.Popen(['./buildSpacedDB.sh'])
	build_spaced.wait()
	##align sequences
	run_multi_clark = subprocess.Popen(['classify_metagenome.sh', '-P', read1, read2, '--spaced', '-m', '0', '-R', results_file])
	run_multi_clark.wait()
	##get summary of those results
	with open(summary_file, "w") as outfh:
		est_ab = subprocess.Popen(['./estimate_abundance.sh', '-F', results_file + '.csv', '-D', clark_db_dir], stdout=outfh)
		est_ab.wait()
'''

##vector i.e custom
##make Cusom dir within clark_db_dir and move UniVec.fa (downloaded from ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec) into it


##map reads corresponding to certain species --- not used
##dict of species serge in interesetd in in format serge name: [clark name, clark id]
# species_reads_to_get = {'mycobacterium tuberculosis': ['mycobacterium tuberculosis', '1773'], 'HHV-7 ': ['Human betaherpesvirus 7', '10372'], 
# 'HCV': ['Hepatitis C virus', '11103'], 'HHV-6B ': ['Human betaherpesvirus 6B', '32604'], 'Campylobacter jejuni': ['Campylobacter jejuni', '197'], 
# 'Staphylococcus haemolyticus': ['Staphylococcus haemolyticus', '1283'], 'human enterovirus ': ['Enterovirus H', '310907']}
# clark_results_files =['DNA_Cortex.default_results.csv', 'RNA_Cortex.default_results.csv', 'DNA_Stem.default_results.csv', 'RNA_Stem.default_results.csv']
# clark_results_files =['DNA_Cortex.default_results.csv']
# for c_file in clark_results_files:
# 	get_read_ids_from_clark_results(species_reads_to_get, c_file)


##deplete human reads and then map against 7 genomes
##manually get all genome and put in one fa file and then index with bwa index 
ref_fa = 'interesting_sequences.fa'
clark_results_files =['DNA_Cortex.default_results.csv', 'RNA_Cortex.default_results.csv', 'DNA_Stem.default_results.csv', 'RNA_Stem.default_results.csv']
# clark_results_files =['DNA_Cortex.default_results.csv']

for c_file in clark_results_files:
	read1, read2 = remove_human_reads(c_file)
	map_against_genomes([read1, read2], ref_fa)











