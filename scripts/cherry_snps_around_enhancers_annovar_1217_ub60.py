#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import pybedtools as pbt 


'''
module load local_python/3.6.5 
'''


##parameters
delim = '\t'
working_dir = '/data/atimms/cherry_snps_around_enhancers_1217'
os.chdir(working_dir)

##programs
gatk = '/home/atimms/programs/GenomeAnalysisTK.jar'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
ann_var='/home/atimms/programs/annovar/annotate_variation.pl'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
bigWigToBedGraph = '/home/atimms/programs/bigWigToBedGraph'
bigWigToWig = '/home/atimms/programs/bigWigToWig'

##annovar folders and parameters
av_genome = 'hg38'
annovar_ref_dir = '/data/atimms/references/annovar/' + av_genome
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/data/atimms/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'gnomad_genome', 'kaviar_20150923,1000g2015aug']
av_operation = ['-operation', 'f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.']



##methods
def download_annovar_dbs(db_name):
	#$ann_var -buildver hg19 -downdb clinvar_20170130 $av_ref_hg19 -webfrom annovar
	dl_ann = subprocess.Popen([ann_var, '-buildver', av_genome,  '-downdb', db_name, annovar_ref_dir, '-webfrom', 'annovar'])
	dl_ann.wait()

def filter_db_file(freq_wanted, gen_db, snp_suffix, indel_suffix, all_suffix):
	snp_outfile = gen_db + snp_suffix
	all_outfile = gen_db + all_suffix
	indel_outfile = gen_db + indel_suffix
	if gen_db == 'gnomad':
		infile = '/data/atimms/references/annovar/hg38/hg38_gnomad_genome.txt'
		with open(infile, "r") as in_fh, open(snp_outfile, "w") as snp_fh, open(all_outfile, "w") as all_fh, open(indel_outfile, "w") as ind_fh:
			lc, spc, apc, ipc = 0,0,0,0
			for line in in_fh:
				lc += 1
				if lc >1:
					line = line.split(delim)
					all_af = line[5]
					ref = line[3]
					alt = line[4]
					##ignore dots
					if all_af != '.':
						all_af = float(all_af)
						if all_af >= freq_wanted:
							apc += 1
							all_fh.write(delim.join(line[:5]) + '\n')
							if len(ref) == 1 and len(alt) == 1 and alt != '-' and ref != '_':
								spc += 1
								snp_fh.write(delim.join(line[:5]) + '\n')
							else:
								ipc += 1
								ind_fh.write(delim.join(line[:5]) + '\n')					
	if gen_db == '1000g':
		infile = '/data/atimms/references/annovar/hg38/hg38_ALL.sites.2015_08.txt'
		with open(infile, "r") as in_fh, open(snp_outfile, "w") as snp_fh, open(all_outfile, "w") as all_fh, open(indel_outfile, "w") as ind_fh:
			lc, spc, apc, ipc = 0,0,0,0
			for line in in_fh:
				lc += 1
				if lc >1:
					line = line.split(delim)
					all_af = line[4]
					ref = line[2]
					alt = line[3]
					##ignore dots
					if all_af != '.':
						all_af = float(all_af)
						chrom = line[0]
						start = str(int(line[1]) -1)
						end = line[1]
						r_a = line[2:4]
						line_out = [chrom, start, end] + r_a
						if all_af >= freq_wanted:
							apc += 1
							all_fh.write(delim.join(line_out) + '\n')
							if len(ref) == 1 and len(alt) == 1:
								spc += 1
								snp_fh.write(delim.join(line_out) + '\n')
							else:
								ipc += 1
								ind_fh.write(delim.join(line_out) + '\n')
	if gen_db == 'kaviar':
		infile = '/data/atimms/references/annovar/hg38/hg38_kaviar_20150923.txt'
		with open(infile, "r") as in_fh, open(snp_outfile, "w") as snp_fh, open(all_outfile, "w") as all_fh, open(indel_outfile, "w") as ind_fh:
			lc, spc, apc, ipc = 0,0,0,0
			for line in in_fh:
				lc += 1
				if lc >1:
					line = line.split(delim)
					all_af = line[5]
					ref = line[3]
					alt = line[4]
					##ignore dots
					if all_af != '.':
						all_af = float(all_af)
						if all_af >= freq_wanted:
							apc += 1
							all_fh.write(delim.join(line[:5]) + '\n')
							if len(ref) == 1 and len(alt) == 1 and alt != '-' and ref != '_':
								spc += 1
								snp_fh.write(delim.join(line[:5]) + '\n')
							else:
								ipc += 1
								ind_fh.write(delim.join(line[:5]) + '\n')	
	# print('%s total variants, of which %s SNPs and %s SNPs and indels and %s indels had an allele frequency of > %s'%(lc, spc, apc, ipc, freq_wanted))


def process_bed_file(ebed, gen_db, snp_suffix, indel_suffix, all_suffix, freq_wanted):
	snp_infile = gen_db + snp_suffix
	all_infile = gen_db + all_suffix
	indel_infile = gen_db + indel_suffix
	temp_bed = 'temp.bed'
	windows_bed = 'windows.bed'
	snp_counts_bed = gen_db + '.' + str(freq_wanted) + '.snp_counts.bed'
	all_counts_bed = gen_db + '.' + str(freq_wanted) + '.all_counts.bed'
	indel_counts_bed = gen_db + '.' + str(freq_wanted) + '.indel_counts.bed'
	####make bed file of the region +- 1000bp
	with open(ebed, "r") as ebed_fh, open(temp_bed, "w") as tbed_fh:
		rc = 0
		for line in ebed_fh:
			rc += 1
			line = line.split(delim)
			chrom = line[0][3:]
			start = int(line[1]) - 1000
			end = int(line[2]) + 1024
			tbed_fh.write(delim.join([chrom, str(start), str(end), str(rc)]) + '\n')
	##make windows bed
	with open(windows_bed, "w") as wbed_fh: 
		##bedtools makewindows -b input.bed -n 3 -i srcwinnum
		bt_mw = subprocess.Popen(['bedtools', 'makewindows', '-b', temp_bed, '-w', '50', '-s', '25', '-i', 'srcwinnum'], stdout=wbed_fh)
		bt_mw.wait()
	##get count data
	with open(snp_counts_bed, "w") as scbed_fh: 
		bt_mw = subprocess.Popen(['bedtools', 'intersect', '-a', windows_bed, '-b', snp_infile, '-wa', '-c'], stdout=scbed_fh)
		bt_mw.wait()
	with open(all_counts_bed, "w") as acbed_fh: 
		bt_mw = subprocess.Popen(['bedtools', 'intersect', '-a', windows_bed, '-b', all_infile, '-wa', '-c'], stdout=acbed_fh)
		bt_mw.wait()
	with open(indel_counts_bed, "w") as icbed_fh: 
		bt_mw = subprocess.Popen(['bedtools', 'intersect', '-a', windows_bed, '-b', indel_infile, '-wa', '-c'], stdout=icbed_fh)
		bt_mw.wait()


def make_counts_file(infile, outfile):
	w_dict = {}
	with open(infile, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			window = int(line[3].split('_')[1])
			count = int(line[4])
			if window in w_dict:
				w_dict[window] += count
			else:
				w_dict[window] = count
	with open(outfile, "w") as out_fh:
		out_fh.write('window' + delim + 'count' + '\n')
		for w in sorted(w_dict.iterkeys()):
			out_fh.write(str(w) + delim + str(w_dict[w]) + '\n')

def combined_count_files(snp_infile, indel_infile, all_infile, outfile):
	c_dict = {}
	with open(snp_infile, "r") as snp_fh:
		lc = 0
		for line in snp_fh:
			lc += 1
			line = line.rstrip().split(delim)
			if lc >1:
				w = line[0]
				c = line[1]
				if int(w) <= 80:
					c_dict[w] = [c]
	with open(indel_infile, "r") as indel_fh:
		lc = 0
		for line in indel_fh:
			lc += 1
			line = line.rstrip().split(delim)
			if lc >1:
				w = line[0]
				c = line[1]
				if int(w) <= 80:
					c_dict[w].append(c)
	with open(all_infile, "r") as all_fh:
		lc = 0
		for line in all_fh:
			lc += 1
			line = line.rstrip().split(delim)
			if lc >1:
				w = line[0]
				c = line[1]
				if int(w) <= 80:
					c_dict[w].append(c)
	with open(outfile, "w") as out_fh:
		out_fh.write('window' + delim + 'count' + '\n')
		for w in sorted(w_dict.iterkeys()):
			counts = [str(i) for i in w_dict[w]]
			out_fh.write(delim.join([w] + counts) + '\n')


def combined_count_files_for_r(snp_infile, indel_infile, all_infile, outfile):
	with open(outfile, "w") as out_fh:
		out_fh.write('window' + delim + 'variant_type' + delim + 'count' + '\n')
		with open(snp_infile, "r") as snp_fh:
			lc = 0
			for line in snp_fh:
				lc += 1
				line = line.rstrip().split(delim)
				if lc >1:
					w = line[0]
					c = line[1]
					if int(w) <= 80:
						out_fh.write(delim.join([w, 'snp', c]) + '\n')
		with open(indel_infile, "r") as indel_fh:
			lc = 0
			for line in indel_fh:
				lc += 1
				line = line.rstrip().split(delim)
				if lc >1:
					w = line[0]
					c = line[1]
					if int(w) <= 80:
						out_fh.write(delim.join([w, 'indel', c]) + '\n')
		with open(all_infile, "r") as all_fh:
			lc = 0
			for line in all_fh:
				lc += 1
				line = line.rstrip().split(delim)
				if lc >1:
					w = line[0]
					c = line[1]
					if int(w) <= 80:
						out_fh.write(delim.join([w, 'snp_indel', c]) + '\n')



def combine_count_bed(gen_db, outfile_suffix, freq_wanted):
	snp_counts_bed = gen_db + '.' + str(freq_wanted) + '.snp_counts.bed'
	all_counts_bed = gen_db + '.' +  str(freq_wanted) + '.all_counts.bed'
	indel_counts_bed = gen_db + '.' +  str(freq_wanted) + '.indel_counts.bed'
	snps_outfile = gen_db + '.snp_counts' + outfile_suffix
	all_outfile = gen_db + '.all_counts' + outfile_suffix
	indel_outfile = gen_db + '.indel_counts' + outfile_suffix
	comb_outfile = gen_db + '.combined_counts' + outfile_suffix
	# make_counts_file(snp_counts_bed, snps_outfile)
	# make_counts_file(all_counts_bed, all_outfile)
	# make_counts_file(indel_counts_bed, indel_outfile)
	# combined_count_files(snps_outfile, indel_outfile, all_outfile, comb_outfile)
	combined_count_files_for_r(snps_outfile, indel_outfile, all_outfile, comb_outfile)

def process_bed_file_gerp(ebed, gerp_bedfile, out_file):
	temp_bed = 'temp.bed'
	temp2_bed = 'temp2.bed'
	windows_bed = 'windows.bed'
	window_gerp_dict = {}
	# '''
	####make bed file of the region +- 1000bp
	with open(ebed, "r") as ebed_fh, open(temp_bed, "w") as tbed_fh:
		rc = 0
		for line in ebed_fh:
			rc += 1
			line = line.split(delim)
			chrom = line[0][3:]
			start = int(line[1]) - 1000
			end = int(line[2]) + 1024
			tbed_fh.write(delim.join([chrom, str(start), str(end), str(rc)]) + '\n')
	##make windows bed
	with open(windows_bed, "w") as wbed_fh: 
		##bedtools makewindows -b input.bed -n 3 -i srcwinnum
		bt_mw = subprocess.Popen(['bedtools', 'makewindows', '-b', temp_bed, '-w', '50', '-s', '25', '-i', 'srcwinnum'], stdout=wbed_fh)
		bt_mw.wait()
	
	##get count data
	with open(temp2_bed, "w") as t2bed_fh: 
		bt_mw = subprocess.Popen(['bedtools', 'intersect', '-a', windows_bed, '-b', gerp_bedfile, '-wao'], stdout=t2bed_fh)
		bt_mw.wait()
	# '''
	##take temp bed and get average gerp per window
	with open(temp2_bed, "r") as t22bed_fh: 
		for line in t22bed_fh:
			# print(line)
			line = line.rstrip().split(delim)
			info = '_'.join(line[:4])
			rs = line[7]
			##removes windows with no data
			if rs != '.':
				if info in window_gerp_dict:
					window_gerp_dict[info].append(float(rs))
				else:
					window_gerp_dict[info] = [float(rs)]
	with open(out_file, "w") as out_fh: 
		for i in window_gerp_dict:
			average_rs = sum(window_gerp_dict[i]) / len(window_gerp_dict[i])
			print(i, window_gerp_dict[i], average_rs)
			i_out = i.split('_')
			out_fh.write(delim.join(i_out + [str(average_rs)]) + '\n')

def process_bed_file_phastcon(ebed, gerp_bedfile, out_file):
	temp_bed = 'temp.bed'
	temp1_bed = 'temp1.bed'
	temp2_bed = 'temp2.bed'
	# temp_pc_bed = 'temp_pc.bed'
	windows_bed = 'windows.bed'
	window_gerp_dict = {}
	# '''
	####make bed file of the region +- 1000bp
	with open(ebed, "r") as ebed_fh, open(temp_bed, "w") as tbed_fh:
		rc = 0
		for line in ebed_fh:
			rc += 1
			line = line.split(delim)
			chrom = line[0][3:]
			start = int(line[1]) - 1000
			end = int(line[2]) + 1024
			tbed_fh.write(delim.join([chrom, str(start), str(end), str(rc)]) + '\n')

	##make windows bed
	with open(temp1_bed, "w") as t1bed_fh: 
		##bedtools makewindows -b input.bed -n 3 -i srcwinnum
		bt_mw = subprocess.Popen(['bedtools', 'makewindows', '-b', temp_bed, '-w', '50', '-s', '25', '-i', 'srcwinnum'], stdout=t1bed_fh)
		bt_mw.wait()
	##sort bed file
	with open(windows_bed, "w") as wbed_fh: 
		bt_sort = subprocess.Popen(['bedtools', 'sort', '-i', temp1_bed], stdout=wbed_fh)
		bt_sort.wait()
	'''
	##sort phastcon file #sorted manually
	with open(temp_pc_bed, "w") as tpcbed_fh: 
		bt_sort = subprocess.Popen(['bedtools', 'sort', '-i', gerp_bedfile], stdout=tpcbed_fh)
		bt_sort.wait()
	'''

	##get count data
	with open(temp2_bed, "w") as t2bed_fh: 
		bt_mw = subprocess.Popen(['bedtools', 'intersect', '-a', windows_bed, '-b', gerp_bedfile, '-wao', '-sorted'], stdout=t2bed_fh)
		bt_mw.wait()
	
	##take temp bed and get average gerp per window
	with open(temp2_bed, "r") as t22bed_fh: 
		for line in t22bed_fh:
			# print(line)
			line = line.rstrip().split(delim)
			info = '_'.join(line[:4])
			rs = line[7]
			##removes windows with no data
			if rs != '.':
				if info in window_gerp_dict:
					window_gerp_dict[info].append(float(rs))
				else:
					window_gerp_dict[info] = [float(rs)]
	with open(out_file, "w") as out_fh: 
		for i in window_gerp_dict:
			average_rs = sum(window_gerp_dict[i]) / len(window_gerp_dict[i])
			# print(i, window_gerp_dict[i], average_rs)
			i_out = i.split('_')
			out_fh.write(delim.join(i_out + [str(average_rs)]) + '\n')

def convert_gerp_to_bed(in_file, out_file):
	check_dict = {}
	with open(in_file, "r") as in_fh, open(out_file, "w") as out_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.rstrip().split(delim)
			if lc >1:
				# print(line)
				chrom = line[0]
				pos = line[2]
				rs = line[5]
				cp = chrom + '_' + pos
				start = str(int(pos) - 1)
				if cp in check_dict:
					if rs != check_dict[cp]:
						print('issue with pos:', cp)
				else:
					out_fh.write(delim.join([chrom, start, pos, rs]) + '\n')
					check_dict[cp] = rs
					# print(chrom, start, pos, rs)

def combine_windows_gerp_bed(infile, outfile):
	w_dict = {}
	with open(infile, "r") as in_fh:
		for line in in_fh:
			print(line)
			line = line.rstrip().split(delim)
			window = int(line[4])
			rs = float(line[5])
			if window in w_dict:
				w_dict[window].append(rs)
			else:
				w_dict[window] = [rs]
	with open(outfile, "w") as out_fh:
		out_fh.write('window' + delim + 'count' + '\n')
		for w in sorted(w_dict.keys()):
			rs_ave = sum(w_dict[w]) / len(w_dict[w])
			print(w, w_dict[w], rs_ave)
			out_fh.write(str(w) + delim + str(rs_ave) + '\n')


def convert_bw_to_bedgraph(in_bw, out_bed):
	bt_mw = subprocess.Popen([bigWigToBedGraph, in_bw, out_bed])
	bt_mw.wait()


##run methods
# enhancer_bed = 'Hu-ret-ATAC-H3K27ac-dedup-GREAT-diseasegenes-ATACsummit.bed'

##download genome databases
'''
dbs_to_download = ['gnomad_genome', 'kaviar_20150923', '1000g2015aug']
dbs_to_download = ['hrcr1']
dbs_to_download = ['hrcr1']
for db in dbs_to_download:
	download_annovar_dbs(db)
'''


'''
##do the thing --parameters
freq_req_list = [0.01, 0.05, 0.001]
genome_dbs = ['1000g', 'gnomad', 'kaviar']
# enhancer_bed = 'test.bed'

for freq_req in freq_req_list:
	genome_db_snps_suffix = '.' + str(freq_req) + 'maf.snps.bed'
	genome_db_snps_indels_suffix = '.' + str(freq_req) + 'maf.snps_indels.bed'
	genome_db_indels_suffix = '.' + str(freq_req) + 'maf.indels.bed'
	counts_file_suffix = '.' + str(freq_req) + 'maf.counts.txt'
	for genome_db in genome_dbs:
		##filter annoatation dbs for query
		# filter_db_file(freq_req, genome_db, genome_db_snps_suffix, genome_db_indels_suffix, genome_db_snps_indels_suffix)
		##process regions from bed and get windows and count data
		# process_bed_file(enhancer_bed, genome_db, genome_db_snps_suffix, genome_db_indels_suffix, genome_db_snps_indels_suffix, freq_req)
		##combine counts from all enhancers
		combine_count_bed(genome_db, counts_file_suffix, freq_req)
'''

##lets do the same for gerp score 
'''
##parameters
enhancer_bed = 'Hu-ret-ATAC-H3K27ac-dedup-GREAT-diseasegenes-ATACsummit_hg19.bed'
# gerp_file = 'test.txt' ##testing
gerp_file = 'hg19_ljb26_gerp++.txt'
gerp_bed = 'hg19_gerp.bed'
gerp_results = 'gerp_window_scores.bed'
gerp_combined_results = 'gerp_window_averages.txt'

##have gerp file from annovar reference
convert_gerp_to_bed(gerp_file, gerp_bed)

##process regions from bed and get windows and count data
process_bed_file_gerp(enhancer_bed, gerp_bed, gerp_results)
##combine counts from all enhancers
combine_windows_gerp_bed(gerp_results, gerp_combined_results)
'''

##now phast con scores

##7way
##parameters
enhancer_bed = 'Hu-ret-ATAC-H3K27ac-dedup-GREAT-diseasegenes-ATACsummit.bed'
phastcon_prefix = 'hg38.phastCons7way'
phastcon_bw = phastcon_prefix + '.bw'
phastcon_bed = phastcon_prefix + '.bed'
# phastcon_bed_nochr = phastcon_prefix + '.no_chr_sorted.bed'
phastcon_bed_nochr = phastcon_prefix + '.no_chr.bed'
phastcon_results = phastcon_prefix + '.results.bed'
phastcon_comb_results = phastcon_prefix + '.comb_results.txt'

##make bed file from bw
# convert_bw_to_bedgraph(phastcon_bw, phastcon_bed)
##then maunually: removed 'chr' >sed 's/^chr//' hg38.phastCons7way.bed > temp_pc7.bed 
##process regions from bed and get windows and count data
# process_bed_file_phastcon(enhancer_bed, phastcon_bed_nochr, phastcon_results)
##combine counts from all enhancers
# combine_windows_gerp_bed(phastcon_results, phastcon_comb_results)


##100way
##parameters
enhancer_bed = 'Hu-ret-ATAC-H3K27ac-dedup-GREAT-diseasegenes-ATACsummit.bed'
phastcon_prefix = 'hg38.phastCons100way'
phastcon_bw = phastcon_prefix + '.bw'
phastcon_bed = phastcon_prefix + '.bed'
phastcon_bed_nochr = phastcon_prefix + '.no_chr.bed'
phastcon_results = phastcon_prefix + '.results.bed'
phastcon_comb_results = phastcon_prefix + '.comb_results.txt'

##make bed file from bw
# convert_bw_to_bedgraph(phastcon_bw, phastcon_bed)
##then maunually: removed 'chr' >sed 's/^chr//' hg38.phastCons7way.bed > temp_pc7.bed 
##process regions from bed and get windows and count data
# process_bed_file_phastcon(enhancer_bed, phastcon_bed_nochr, phastcon_results)
##combine counts from all enhancers
# combine_windows_gerp_bed(phastcon_results, phastcon_comb_results)


##100way phylop
##parameters
enhancer_bed = 'Hu-ret-ATAC-H3K27ac-dedup-GREAT-diseasegenes-ATACsummit.bed'
phastcon_prefix = 'hg38.phyloP100way'
phastcon_bw = phastcon_prefix + '.bw'
phastcon_bed = phastcon_prefix + '.bed'
phastcon_bed_nochr = phastcon_prefix + '.no_chr.bed'
phastcon_results = phastcon_prefix + '.results.bed'
phastcon_comb_results = phastcon_prefix + '.comb_results.txt'

##make bed file from bw
# convert_bw_to_bedgraph(phastcon_bw, phastcon_bed)
##then maunually: removed 'chr' >sed 's/^chr//' hg38.phastCons7way.bed > temp_pc7.bed 
##process regions from bed and get windows and count data
process_bed_file_phastcon(enhancer_bed, phastcon_bed_nochr, phastcon_results)
##combine counts from all enhancers
combine_windows_gerp_bed(phastcon_results, phastcon_comb_results)



