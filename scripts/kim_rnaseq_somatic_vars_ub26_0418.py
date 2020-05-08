#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import filtering_annotated
import shutil
# import dobyns_gemini_pipeline_v2

##parameters
delim = '\t'
threads = '16'

##programs
star = '/home/atimms/programs/STAR/bin/Linux_x86_64/STAR'
feature_counts = '/home/atimms/programs/subread-1.5.1-Linux-x86_64/bin/featureCounts'
picard = '/home/atimms/programs/picard-tools-1.139/picard.jar'
gatk = '/home/atimms/programs/GenomeAnalysisTK.jar'
freebayes = '/home/atimms/programs/freebayes/bin/freebayes'
bcftools_12 = '/home/atimms/programs/bcftools-1.2/bcftools'
vt = '/home/atimms/programs/vt/vt'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
pisces = '/opt/Pisces_5.1.6.54/Pisces.exe'

##ref files and file names
genome_name = 'hg19'
star_index_dir = '/data/atimms/references/star/' + genome_name
fa_file = '/data/atimms/references/igenomes/' + genome_name + '/genome.fa'
gtf_file = '/data/atimms/references/igenomes/' + genome_name +'/genes.gtf'
star_bam_suffix = '.Aligned.out.bam'
proccessed_bam_suffix = '.star.gatk.bam'
ref_dir = '/data/atimms/references/'
indels_mills = ref_dir + 'Mills_and_1000G_gold_standard.indels.hg19.sites.vcf'
dbsnp = ref_dir + 'dbsnp_138.hg19.vcf'
# indels_1000g = ref_dir + '1000G_phase1.indels.b37.vcf'
pisces_fastq_dir = '/data/atimms/references/illumina_hg19/'
bamlist = 'bams.list'

#annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/data/atimms/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,dbnsfp30a,popfreq_all_20150413,avsnp147,generic,generic,generic,generic', '-genericdbfile']
av_operation = ['-operation', 'g,r,r,f,f,f,f,f,f,f']
# av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,dbnsfp30a,popfreq_all_20150413,avsnp147']
# av_operation = ['-operation', 'g,r,r,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput']

##definitions for determining var function
exonic_definitions = ['exonic', 'splicing']
non_protein_changing_definitions =['synonymous SNV', 'synonymous_SNV']
no_record_definition = ['', '.']
##list of columns for both methods (all start at 1):
##for filtering variants: refgene function, frequency dbs, refgene exonic function, refgene geneame, quality and coverage
filtering_cols = [15, 56, 18, 16, [7,8]]
frequencies = [0.01]

def make_star_index_files(star_genome_dir, genome_fas, genome_gtf, threads_to_use):
	star_index = subprocess.Popen([star, '--runMode', 'genomeGenerate', '--genomeDir', star_genome_dir, '--genomeFastaFiles', genome_fas, '--sjdbGTFfile', genome_gtf, '--runThreadN', threads_to_use])
	star_index.wait()

def star_align_paired_end_2_pass(sample_list, star_genome_dir, threads_to_use, fq_prefixes):
	for sample_name in sample_list:
		r1_fq = sample_name + fq_prefixes[0]
		r2_fq = sample_name + fq_prefixes[1]
		print sample_name, r1_fq, r2_fq
		star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', r1_fq, r2_fq, '--outSAMtype', 'BAM', 'Unsorted', '--outSAMmapqUnique', '60', '--readFilesCommand', 'zcat', '--twopassMode', 'Basic', '--outFileNamePrefix', sample_name + '.', '--runThreadN', threads_to_use])
		star_align.wait()



def proccess_bam_files(sample_names):
	for sample_name in sample_names:
		in_bam = sample_name + star_bam_suffix
		arg_bam = sample_name + '.with_rg.bam'
		mkdup_bam = sample_name + '.mkdup.bam'
		reorder_bam = sample_name + '.reorder.bam'
		split_bam = sample_name + '.split.bam'
		final_bam = sample_name + proccessed_bam_suffix
		picard_arg = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'AddOrReplaceReadGroups', 'INPUT=' + in_bam, 'OUTPUT=' + arg_bam, 'SO=coordinate', 'RGID=' + sample_name, 'RGLB=' + sample_name, 'RGPL=Illumina', 'RGPU=machine1', 'RGSM=' + sample_name ])
		picard_arg.wait()	
		picard_md = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=SILENT', 'INPUT=' + arg_bam, 'OUTPUT=' + mkdup_bam,'CREATE_INDEX=true', 'METRICS_FILE=' + sample_name + '.metrics'])
		picard_md.wait()
		picard_rs = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'ReorderSam', 'INPUT=' + mkdup_bam, 'OUTPUT=' + reorder_bam, 'CREATE_INDEX=true', 'REFERENCE=' + fa_file])
		picard_rs.wait()
		gatk_snt = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'SplitNCigarReads',  '-R', fa_file, '-I', reorder_bam, '-o', split_bam, '-U', 'ALLOW_N_CIGAR_READS'])
		gatk_snt.wait()
		gatk_br = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'BaseRecalibrator', '-nct', '6', '-R', fa_file, '-I', split_bam, '-knownSites', indels_mills, '-knownSites', dbsnp, '-o', sample_name + '.recal_data.table'])
		gatk_br.wait()
		gatk_pr = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'PrintReads', '-nct', '6', '-R', fa_file, '-I', split_bam, '-BQSR', sample_name + '.recal_data.table', '-o', final_bam])
		gatk_pr.wait()
		##remove intermeidate files
		files_to_go = glob.glob('*with_rg.bam') + glob.glob('*mkdup.bam') + glob.glob('*mkdup.bai') + glob.glob('*metrics') + glob.glob('*reorder.bam')+ glob.glob('*reorder.bai') + glob.glob('*split.bam') + glob.glob('*split.bai')
		print 'removing files:'
		for f in files_to_go:
			os.remove(f)
			print f


def variant_calling_pisces(ped_name, samples, bam_suffix):
	out_dir = ped_name + '_pisces'
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
	print bam_files
	##run picses on all bam
	# run_pisces = subprocess.Popen(['mono', 'Pisces', '-G', pisces_fastq_dir, '-BamPaths', ','.join(bam_files), '-MinVF', '0.01', '-i', exome_capture_bed, '-t', '16', '-OutFolder', out_dir])
	##pisces now loaded by module so can just call Pisces and it adds the mono
	
	# run_pisces = subprocess.Popen(['Pisces', '-G', pisces_fastq_dir, '-BamPaths', ','.join(bam_files), '-MinVF', '0.01', '-i', exome_capture_bed, '-t', '16', '-OutFolder', out_dir])
	run_pisces = subprocess.Popen(['mono', pisces, '-G', pisces_fastq_dir, '-BamPaths', ','.join(bam_files), '-MinVF', '0.01', '-t', '16', '-OutFolder', out_dir])
	run_pisces.wait()


def combine_vcf_files(vcfs, out_vcf):
	gz_vcfs = []
	for vcf in vcfs:
		bgzip_vcf = subprocess.Popen(['bgzip', vcf])
		bgzip_vcf.wait()
		bcf_index = subprocess.Popen([bcftools_12, 'index', vcf + '.gz'])
		bcf_index.wait()
		gz_vcfs.append(vcf + '.gz')
	bcf_merge = subprocess.Popen([bcftools_12, 'merge'] + gz_vcfs + ['-O', 'v', '-o', out_vcf, '-m', 'none'])
	bcf_merge.wait()

def combine_vcf_files_gatk(vcfs, out_vcf):
	vcfs_with_v = []
	for vcf in vcfs:
		vcf_with_v = ['-V', vcf]
		vcfs_with_v.extend(vcf_with_v)
	print vcfs_with_v
	# combine_var = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineVariants', '-R',fa_file, '-nt', '15', '--variant', vcfs[0],'--variant', vcfs[1],'-o', out_vcf, '-genotypeMergeOptions', 'UNIQUIFY'])
	combine_var = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineVariants', '-R',fa_file, '-nt', '15'] + vcfs_with_v + ['-o', out_vcf, '-genotypeMergeOptions', 'REQUIRE_UNIQUE'])
	combine_var.wait()

def intesect_two_vcf_files(vcf1, vcf2, output_dir):
	for vcf in [vcf1, vcf2]:
		bgzip_vcf = subprocess.Popen(['bgzip', vcf])
		bgzip_vcf.wait()
		bcf_index = subprocess.Popen(['bcftools', 'index', vcf + '.gz'])
		bcf_index.wait()
	bcf_isec = subprocess.Popen(['bcftools', 'isec', vcf1 + '.gz', vcf2 + '.gz', '-p', output_dir ])
	bcf_isec.wait()

##convert vcf file to individual avinput file
def convert_to_annovar(vcf, project_prefix):
	av_prefix = project_prefix
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', av_prefix])
	con_ann.wait()

##annotate vcf file
def table_annovar_vcf(vcf, project_prefix, av_files):
	out_prefix = project_prefix
	command = [table_annovar] + av_buildver + [vcf] + av_ref_dir + av_protocol + [','.join(av_files)] + av_operation + av_options + ['-out', out_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()


	
def format_avinput(project_prefix):
	multianno = project_prefix + '.hg19_multianno.txt'
	outfile = project_prefix + '.annotated.txt'
	head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 'genomicSuperDups', 'SIFT_score', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score', 'Polyphen2_HVAR_pred', 'LRT_score', 'LRT_pred', 'MutationTaster_score', 'MutationTaster_pred', 'MutationAssessor_score', 'MutationAssessor_pred', 'FATHMM_score', 'FATHMM_pred', 'PROVEAN_score', 'PROVEAN_pred', 'VEST3_score', 'CADD_raw', 'CADD_phred', 'DANN_score', 'fathmm-MKL_coding_score', 'fathmm-MKL_coding_pred', 'MetaSVM_score', 'MetaSVM_pred', 'MetaLR_score', 'MetaLR_pred', 'integrated_fitCons_score', 'integrated_confidence_value', 'GERP++_RS', 'phyloP7way_vertebrate', 'phyloP20way_mammalian', 'phastCons7way_vertebrate', 'phastCons20way_mammalian', 'SiPhy_29way_logOdds', 'PopFreqMax', '1000G_ALL', '1000G_AFR', '1000G_AMR', '1000G_EAS', '1000G_EUR', '1000G_SAS', 'ExAC_ALL', 'ExAC_AFR', 'ExAC_AMR', 'ExAC_EAS', 'ExAC_FIN', 'ExAC_NFE', 'ExAC_OTH', 'ExAC_SAS', 'ESP6500siv2_ALL', 'ESP6500siv2_AA', 'ESP6500siv2_EA', 'CG46', 'avsnp147', 'na', 'na', 'coverage', 'chr', 'pos', 'id', 'ref2', 'alt2', 'qual', 'filter', 'info', 'format', 'H25578.adrenal.1435', 'H25578.brain.1436', 'H25578.heart.1437', 'H25578.kidney.1438', 'H25578.liver.1439', 'H25578.lung.1440', 'H25578.maternal.1433', 'H25578.placenta.1434', 'H25578.skin.1441']
	head_out = delim.join(head + ['\n'])
	with open(multianno, "r") as multi_fh, open(outfile, "w") as out_fh:
		out_fh.write(head_out)
		line_count = 0
		for line in multi_fh:
			line_count += 1
			if line_count >1:
				out_fh.write(line)

def split_info_field(info_list):
	indices = [i for i, s in enumerate(info_list) if 'ANNOVAR_DATE' in s]
	#print indices
	i_count = 0
	final_list = []
	for info in info_list:
		#print info
		if i_count > indices[0] and info != 'ALLELE_END':
			info2 = info.split('=')[1]
			#print info2
			final_list.append(info2)
		i_count += 1
	return final_list

def format_avinput(description, sample_list, av_files, av_suffix):
	for sample in sample_list:
		names = []
		for av in av_files:
			name = av.split('.')[1]
			names.append(name)
		head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Zygosity', 'Qual', 'Coverage', 'Filter', 'Pos', 'Ref2', 'Alt2', 'Format', 'Info', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 'genomicSuperDups', 'SIFT_score', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score', 'Polyphen2_HVAR_pred', 'LRT_score', 'LRT_pred', 'MutationTaster_score', 'MutationTaster_pred', 'MutationAssessor_score', 'MutationAssessor_pred', 'FATHMM_score', 'FATHMM_pred', 'PROVEAN_score', 'PROVEAN_pred', 'VEST3_score', 'CADD_raw', 'CADD_phred', 'DANN_score', 'fathmm-MKL_coding_score', 'fathmm-MKL_coding_pred', 'MetaSVM_score', 'MetaSVM_pred', 'MetaLR_score', 'MetaLR_pred', 'integrated_fitCons_score', 'integrated_confidence_value', 'GERP++_RS', 'phyloP7way_vertebrate', 'phyloP20way_mammalian', 'phastCons7way_vertebrate', 'phastCons20way_mammalian', 'SiPhy_29way_logOdds', 'PopFreqMax', '1000G_ALL', '1000G_AFR', '1000G_AMR', '1000G_EAS', '1000G_EUR', '1000G_SAS', 'ExAC_ALL', 'ExAC_AFR', 'ExAC_AMR', 'ExAC_EAS', 'ExAC_FIN', 'ExAC_NFE', 'ExAC_OTH', 'ExAC_SAS', 'ESP6500siv2_ALL', 'ESP6500siv2_AA', 'ESP6500siv2_EA', 'CG46', 'avsnp147'] + names
		head_out = delim.join(head + ['\n'])
		avinput = description + '.' + sample + av_suffix
		outfile = sample + '.annotated.txt'
		with open(avinput, "r") as av, open(outfile, "w") as final:
			final.write(head_out)
			for line in av:
				line = line.strip('\n').split(delim)
				stuff = line[0:8] + [line[14]] + [line[9]] + line[11:13]
				info = line[15].split(';')
				info_list = split_info_field(info)
				other_stuff = line[16:]
				line_out = delim.join(stuff + other_stuff + info_list +['\n'])
				final.write(line_out)

def move_to_av_folder(project_prefix):
	avinputs = glob.glob('temp.' + project_prefix + '*put')
	# print avinputs
	for put in avinputs:
		cp_vcf = subprocess.Popen(['cp', put, str(av_ref_dir[0])])
		cp_vcf.wait()
	return avinputs

##calls all annovar methods
def annotate_vcf_file(vcf, project_prefix, samples, av_suffix):
	convert_to_annovar(vcf, 'temp')
	avs = move_to_av_folder(project_prefix)
	table_annovar_vcf(vcf, project_prefix, avs)
	post_annovar_vcf = project_prefix + '.' + av_genome + '_multianno.vcf'
	convert_to_annovar(post_annovar_vcf, project_prefix)
	format_avinput(project_prefix, samples, avs, av_suffix)


def format_annovar_multi(multianno, samplenames, av_files):
	print multianno
	print samplenames
	print av_files
	outfile = multianno.split('.')[0] + '.picses_gatk.annotated.xls'
	fifty_cov_out = multianno.split('.')[0] + '.picses_gatk.min_cov50.annotated.xls'
	gt_names = [s + ' genotype' for s in samplenames]
	cov_names = [s + ' coverage' for s in samplenames]
	aaf_names = [s + ' aaf' for s in samplenames]
	print gt_names
	with open(multianno, "r") as multi, open(outfile, "w") as final, open(fifty_cov_out, "w") as final_50:
		line_count = 0
		for line in multi:
			line_count += 1
			line = line.strip('\n').split(delim)
			if line_count == 1:
				head_out = line[:66] + gt_names + cov_names + aaf_names + ['min coverage', '\n']
				final.write(delim.join(head_out))
				final_50.write(delim.join(head_out))
			else:
				genotypes = line[82:86]
				coverages, cov_ints, aafs = [], [], []
				for g in genotypes:
					g = g.split(':')
					if g[0] == './.':
						coverages.append('na')
						cov_ints.append(0)
						aafs.append('no data')
					else:
						print g, line[:5]
						
						coverage = g[2]
						if coverage == '.':
							coverage = '0'
						coverages.append(coverage)
						cov_int = int(coverage)
						cov_ints.append(cov_int)
						alleles = g[1]
						# print alleles, coverage
						if alleles == '.':
							aaf = 'no info'
						else:
							alleles = alleles.split(',')
							if len(alleles) != 2:
								aaf = 'multiallelic'
							else:
								aaf = str(float(alleles[1]) / cov_int)
							aafs.append(aaf)
				min_cov = min(cov_ints)
				print coverages, aafs, min_cov
				line_out = line[:70] + coverages + aafs + [str(min_cov), '\n']
				final.write(delim.join(line_out))
				if min_cov >= 50:
					final_50.write(delim.join(line_out))




def annotate_vcf_file_gatk(vcf, project_prefix, samples):
	convert_to_annovar(vcf, 'temp')
	avs = move_to_av_folder(project_prefix)
	table_annovar_vcf(vcf, project_prefix, avs)
	post_annovar_txt = project_prefix + '.' + av_genome + '_multianno.txt'
	format_annovar_multi(post_annovar_txt, samples, avs)

##simple list maker -- give string and dictionary to get length from
def make_list(string, dict):
	l = []
	for i in range(len(dict)):
		l.append(string)
	return l

def make_list_of_bams(bam_files, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for bam in bam_files:
			outfh.write(bam + '\n')

def genotype_vars_with_gatk(bamlist, combined_vcf, gatk_norm_vcf):
	vcf_temp1 = 'temp_21.vcf'
	vcf_temp2 = 'temp_22.vcf.gz'
	##remove old version
	if os.path.isfile(vcf_temp1 + '.gz'):
		os.remove(vcf_temp1 + '.gz')
	if os.path.isfile(vcf_temp2):
		os.remove(vcf_temp2)
	gatk_ug = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'UnifiedGenotyper', '-R', fa_file, '-gt_mode', 'GENOTYPE_GIVEN_ALLELES', '-alleles', combined_vcf, '-dcov', '2500', '-out_mode', 'EMIT_ALL_SITES', '-I', bamlist, '-o', vcf_temp1])
	gatk_ug.wait()
	##split multi-allelic variants calls in separate lines, and left normalize indels
	bgzip_run = subprocess.Popen(['bgzip', vcf_temp1])
	bgzip_run.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-O', 'z', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fa_file, '-O', 'v', '-o', gatk_norm_vcf, vcf_temp2])
	bcf_norm2.wait()

def filter_ann_txt_files(samples, cols_to_filter, freq_req):
	##filter rnaseq data
	for sample in samples:
		for freq in freq_req:
			##only 'rare' (<=1%)
			filtering_annotated.filter(working_dir, "and", sample + '.annotated.txt', "1.temp", [cols_to_filter[1]], ['<='], [freq])
			##q>=30 and coverage >=5
			filtering_annotated.filter(working_dir, "and", "1.temp", "2.temp", cols_to_filter[4], ['>=', '>='], [30,5])
			##exonic_variants in refGene
			filtering_annotated.filter(working_dir, "or", "2.temp", "3.temp", [cols_to_filter[0], cols_to_filter[0]], ['==','=='], [exonic_definitions[0],exonic_definitions[1]])
			##get all protein changing
			filtering_annotated.filter(working_dir, "and", "3.temp", sample + '.' + str(freq) +  ".protein_changing.xls", make_list(cols_to_filter[2], non_protein_changing_definitions), make_list('!=', non_protein_changing_definitions), non_protein_changing_definitions)

def make_sample_fq_dict(sample_file):
	sample_fq_dict = {}
	with open(sample_file, "r") as infh:
		line_count = 0
		for line in infh:
			line_count += 1
			if line_count > 1:
				line = line.strip('\n').split(delim)
				sample = line[0]
				fq1 = line[1]
				fq2 = line[2]
				if sample in sample_fq_dict:
					sample_fq_dict[sample][0].append(fq1)
					sample_fq_dict[sample][1].append(fq2)
				else:
					sample_fq_dict[sample] = [[fq1], [fq2]]
	return sample_fq_dict

def filer_by_presence_in_other_sample(samples, infile_suffix, outfile_suffix):
	##make dict of all vars
	var_dict = {}
	for sample in samples:
		infile = sample + infile_suffix
		with open(infile, "r") as in_fh:
			line_count = 0
			for line in in_fh:
				line = line.rstrip().split(delim)
				line_count += 1
				if line_count > 1:
					var = '_'.join(line[:5])
					if var in var_dict:
						if sample not in var_dict[var]:
							var_dict[var].append(sample)
					else:
						var_dict[var] = [sample]
	##check dict
	# for v in var_dict:
	# 	print v, var_dict[v]
	##go through each file and check that var is unique and 
	for sample in samples:
		var_count_dict = {}
		infile = sample + infile_suffix
		outfile = sample + outfile_suffix
		with open(infile, "r") as in_fh, open(outfile, "w") as out_fh:
			line_count, unique_count = 0,0
			for line in in_fh:
				line = line.rstrip().split(delim)
				line_count += 1
				if line_count == 1:
					out_fh.write(delim.join(line + ['in other fetus', '\n']))
				else:
					var2 = '_'.join(line[:5])
					var_sample_no = len(var_dict[var2])
					##so if just in one sample
					if var_sample_no == 1:
						unique_count +=1
						in_other_fetus = 'na'
					else:
						other_samples = var_dict[var2][:]
						if sample in other_samples:
							other_samples.remove(sample)
						# print var_dict[var2], other_samples
						in_other_fetus = ','.join(other_samples)
					##if not a duplicate
					if var2 not in var_count_dict:
						out_fh.write(delim.join(line + [in_other_fetus, '\n']))
					var_count_dict[var2] = 0
			print 'for file %s we have %s variants of which %s were unique to this sample'%(infile, line_count - 1, unique_count,)
##make avinput files


##run methods
##make star index files for hg19 (only need to do once)
# make_star_index_files(star_index_dir, fa_file, gtf_file, threads)



##setup working directory where results will be
working_dir = '/data/atimms/kim_rnaseq_somatic_vars_0418'
os.chdir(working_dir)
project_name = 'somatic_vars_0418'
sample_fq_dict = {}
new_sample_names = ['H24551_EGL', 'H24551_PK', 'H24551_RL', 'H24551_Whole', 'H24616_EGL', 'H24616_PK', 'H24616_RL', 'H24616_Whole']
old_sample_names = ['H26566_PK', 'H26566_RL', 'H26566_Whole', 'H26857_PK', 'H26857_RL', 'H26362_RL', 
			'H26362_Whole', 'H26362_PK', 'H26857_Whole', 'H26857_EGL', 'H26566_EGL', 'H26362_EGL']
sample_names = old_sample_names + new_sample_names
samples = ['H26362', 'H26566', 'H26857', 'H24551', 'H24616']
tissue_names = ['EGL', 'PK', 'RL', 'Whole']
fq_suffixes = ['.r1.fq.gz', '.r2.fq.gz']

##proccess bam files
# star_align_paired_end_2_pass(new_sample_names, star_index_dir, threads, fq_suffixes)
# proccess_bam_files(new_sample_names)
##variant calling
# variant_calling_pisces(project_name, sample_names, proccessed_bam_suffix)

# '''
for sample in samples:
	vcfs, s_names, bam_files = [], [], []
	combined_vcf = sample + '.all_tissues.star.gatk.vcf'

	for tissue in tissue_names:
		vcf =  project_name + '_pisces/' + sample + '_' + tissue + '.star.gatk.vcf'
		full_name = sample + '_' + tissue
		vcfs.append(vcf)
		s_names.append(full_name)
		bam = sample + '_' + tissue + '.star.gatk.bam'
		bam_files.append(bam)
	#first take -- just annotate vars, and filter
	# combine_vcf_files_gatk(vcfs, combined_vcf)
	# annotate_vcf_file(combined_vcf, sample, s_names, '.star.gatk.bam.avinput')
	# filter_ann_txt_files(s_names, filtering_cols, frequencies)
	#take 2 -- compare
	# make_list_of_bams(bam_files, bamlist)
	# genotype_vars_with_gatk(bamlist, sample + '.all_tissues.star.gatk.vcf', sample + '.pisces_gatk_genotypes.vcf')
	annotate_vcf_file_gatk(sample + '.pisces_gatk_genotypes.vcf', sample, s_names)
# '''

##add extra info i.e. in other fetus and remove dups
ann_txt_suffix = '.picses_gatk.annotated.xls'
new_ann_txt_suffix = '.picses_gatk.annotated.reformatted.xls'

filer_by_presence_in_other_sample(samples, ann_txt_suffix, new_ann_txt_suffix)








