#!/tools/BioBuilds-2015.04/bin/python
import sys
import subprocess
import os
import glob
# import dobyns_gemini_pipeline_v2

##parameters
delim = '\t'


##file names etc
ref_dir = '/data/atimms/references/'
fasta = ref_dir + 'human_g1k_v37.fasta'
indels_mills = ref_dir + 'Mills_and_1000G_gold_standard.indels.b37.vcf'
indels_1000g = ref_dir + '1000G_phase1.indels.b37.vcf'
rtc_intervals = ref_dir + 'Mills_1000G.intervals'
dbsnp = ref_dir + 'dbsnp_138.b37.vcf'
final_bam_suffix = '.bwa_gatk.bam'
bamlist = 'bams.list'
exome_bed_for_coverage = '/data/atimms/references/dobyns_exome.in_all_targets.1015.bed'
delly_exclude_regions = '/home/atimms/programs/delly/excludeTemplates/human.hg19.excl.tsv'


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
delly = '/home/atimms/programs/delly_v0.7.6_CentOS5.4_x86_64bit'
##gemini programs and ref files
vt = '/data/atimms/gemini/vt/vt'
gemini = '/data/atimms/gemini/bin/gemini'
fasta = '/data/atimms/references/human_g1k_v37.fasta'
snpeff_jar = '/data/atimms/gemini/snpEff/snpEff.jar'
eadbgzip = '/tools/BioBuilds-2015.04/bin/bgzip'
tabix = '/tools/BioBuilds-2015.04/bin/tabix'



def relign_target_creater_knowns():
	gatk_rtc = subprocess.Popen(['java', '-Xmx3g', '-jar', gatk, '-T', 'RealignerTargetCreator', '-R', fasta, '-nt', '15', '-o', rtc_intervals, '-known', indels_mills, '-known', indels_1000g])
	gatk_rtc.wait()

def convert_bam_fastq_bedtools(bamfile, r1_fastq, r2_fastq):
	read_sorted_bam = bamfile[:-4] + 'n_sorted.bam'
	st_n_sort = subprocess.Popen([samtools, 'sort', '-nO', 'bam', '-o', read_sorted_bam, '-@', '10', '-m', '10G', '-T', 'tempy', bamfile])
	st_n_sort.wait()
	bam_fq = subprocess.Popen([bedtools, 'bamtofastq', '-i', read_sorted_bam, '-fq', r1_fastq, '-fq2', r2_fastq])
	bam_fq.wait()

def convert_bam_fastq_picard(bamfile, r1_fastq, r2_fastq):
	picard_sam_fq = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'SamToFastq', 'INPUT=' + bamfile, 'FASTQ=' + r1_fastq, 'SECOND_END_FASTQ=' + r2_fastq, 'VALIDATION_STRINGENCY=SILENT'])
	picard_sam_fq.wait()

def align_with_bwa_one_at_time(sample, r1_fq, r2_fq):
	rg = '@RG\\tID:' + sample + '_RG\\tSM:' + sample + "\\tPL:ILLUMINA"
	post_bwa_bam = sample + '.bwa.bam'
	sort_bam = sample + '.bwa_sort.bam'
	mkdup_bam = sample + '.bwa_mkdup.bam'
	realigned_bam = sample + '.bwa_religned.bam'
	gatk_bam = sample + final_bam_suffix
	# mkdup_bam = sample + '.bwa_mkdup.bam'
	##bwa and convert to bam
	bwa_pe = subprocess.Popen([bwa, 'mem', '-M', '-t', '30', '-R', rg, fasta, r1_fq, r2_fq], stdout=subprocess.PIPE)
	st_sam_bam_pe = subprocess.Popen([samtools, 'view', '-b', '-@', '5', '-o', post_bwa_bam, '-'], stdin=bwa_pe.stdout)
	st_sam_bam_pe.wait()
	##sort bam
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', sample, '-o', sort_bam, '-@', '10', '-m', '10G', post_bwa_bam])
	st_sort_pe.wait()
	##mark duplicates
	picard_md = subprocess.Popen(['java', '-Xmx100g', '-Djava.io.tmpdir=./tmp', '-jar', picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=LENIENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + mkdup_bam, 'METRICS_FILE=' + sample + '.metrics', 'TMP_DIR=./tmp'])
	picard_md.wait()
	##realign around indels
	gatk_ir = subprocess.Popen(['java', '-Xmx80g', '-jar', gatk, '-T', 'IndelRealigner', '-R', fasta, '-I', mkdup_bam, '-targetIntervals', rtc_intervals, '--consensusDeterminationModel', 'KNOWNS_ONLY', '-LOD', '0.4', '-o', realigned_bam, '-known', indels_mills, '-known', indels_1000g])
	gatk_ir.wait()
	##bqsr
	gatk_br = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'BaseRecalibrator', '-nct', '8', '-R', fasta, '-I', realigned_bam, '-knownSites', indels_mills, '-knownSites', dbsnp, '-knownSites', indels_1000g, '-o', sample + '.recal_data.table'])
	gatk_br.wait()
	gatk_pr = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'PrintReads', '-nct', '8', '-R', fasta, '-I', realigned_bam, '-BQSR', sample + '.recal_data.table', '-o', gatk_bam])
	gatk_pr.wait()


##delete files supplied in list
def delete_unwated_files(file_extensions):
	files_to_delete = []
	for ext in file_extensions:
		files = glob.glob(ext)
		files_to_delete.extend(files)
	for f in files_to_delete:
		os.remove(f)

##make list of all bam files to be analyzed
def make_list_of_bams(sample_dict, bam_suffix, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for sample in sample_dict:
			bam = sample + bam_suffix
			outfh.write(bam + '\n')

##call freebayes in all sample listed in list of bamfiles
def variant_calling_freebayes(bamlist, name_prefix):
	vcf_temp1 = 'temp_fb1.vcf'
	vcf_temp2 = 'temp_fb2.vcf.gz'
	final_vcf = name_prefix + '.freebayes.vcf.gz'
	vcf_fh = open(vcf_temp1, 'w')
	freebayes_run = subprocess.Popen([freebayes,'-f' , fasta, '-L', bamlist], stdout=vcf_fh)
	freebayes_run.wait()
	vcf_fh.close()
	bgzip_run = subprocess.Popen([bgzip, vcf_temp1])
	bgzip_run.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines
	print bcftools, vcf_temp1
	bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
	bcf_index.wait()


def variant_calling_gatk_hc_no_gvcf(bamlist, name_prefix):
	vcf_temp0 = 'temp_gatk0.vcf'
	vcf_raw_snps = 'temp_raw_snps.vcf'
	vcf_filtered_snps = 'temp_filtered_snps.vcf'
	vcf_raw_indels = 'temp_raw_indels.vcf'
	vcf_filtered_indels = 'temp_filtered_indels.vcf'
	vcf_temp1 = 'temp_gatk1.vcf'
	vcf_temp2 = 'temp_gatk2.vcf.gz'
	final_vcf = name_prefix + '.gatkHC.vcf.gz'
	bam_files = []
	##run haplotype caller
	gatk_hc = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fasta, '-I', bamlist, '-stand_call_conf', '30', '-stand_emit_conf', '10', '-o', vcf_temp0])
	gatk_hc.wait()
	##split data in SNPs and indels and apply manual variant filtering (add -L)
	snp_cut = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'SelectVariants', '-R',fasta, '-nt', '5', '-V', vcf_temp0, '-o', vcf_raw_snps, '-selectType', 'SNP'])
	snp_cut.wait()
	snp_vf = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'VariantFiltration', '-R',fasta, '-V', vcf_raw_snps, '-o', vcf_filtered_snps, '--filterExpression ', "QD < 2.0 || FS > 60.0 || SOR > 4.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0", '--filterName', "snp_filter"])
	snp_vf.wait()
	indel_cut = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'SelectVariants', '-R',fasta, '-nt', '5', '-V', vcf_temp0, '-o', vcf_raw_indels, '-selectType', 'INDEL'])
	indel_cut.wait()
	indel_vf = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'VariantFiltration', '-R',fasta, '-V', vcf_raw_indels, '-o', vcf_filtered_indels, '--filterExpression ', "QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0", '--filterName', "indel_filter"])
	indel_vf.wait()
	##combine filtered snps and indels
	combine_var = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineVariants', '-R',fasta, '-nt', '5', '--variant', vcf_filtered_snps,'--variant', vcf_filtered_indels,'-o', vcf_temp1, '-genotypeMergeOptions', 'UNSORTED'])
	combine_var.wait()
	bgzip_cmd = subprocess.Popen([bgzip, vcf_temp1])
	bgzip_cmd.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines, and left normalize indels
	bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-O', 'z', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
	bcf_index.wait()

def intesect_two_vcf_files(vcf1, vcf2, output_dir):
	bcf_isec = subprocess.Popen([bcftools, 'isec', vcf1, vcf2, '-p', output_dir ])
	bcf_isec.wait()

##get exome coverage for bam files
def calculate_exome_coverage(samples, bam_suffix, exome_bed, prefix):
	##parameters and header
	coverage_required = [1,5,10,20,50,100,200]
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
			# bedtools_cov = subprocess.Popen([bedtools, 'coverage', '-a', exome_bed, '-b', bam, '-hist', '-sorted'], stdout=hist_fh)
			# bedtools_cov.wait()
			bedtools_cov = subprocess.Popen([bedtools, 'coverage', '-a', exome_bed, '-b', bam, '-hist', '-sorted'], stdout=subprocess.PIPE)
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

def plink_relatadness_check(vcf, file_prefix):
	##bzgip vcf file
	if os.path.isfile(vcf):
		bgzip_cmd = subprocess.Popen([bgzip, vcf])
		bgzip_cmd.wait()
	##correct filtering??
	bcftools_filter = subprocess.Popen([bcftools, 'view', '-f', 'PASS', '-i', "QUAL>50 & MIN(DP)>50", '-V', 'indels', '-o', 'temp_plink.vcf.gz', '-O', 'z', vcf + '.gz'])
	bcftools_filter.wait()
	##generate plink file from vcf
	make_plink = subprocess.Popen([plink, '--vcf', 'temp_plink.vcf.gz', '--vcf-filter', '--const-fid', '--out', 'temp.pass_q50_dp50'])
	make_plink.wait()
	##check sex -results in .sexcheck file
	plink_sex = subprocess.Popen([plink, '--bfile', 'temp.pass_q50_dp50', '--check-sex', '--out', file_prefix + '.pass_q50_dp50'])
	plink_sex.wait()
	##ibd check
	plink_ibd = subprocess.Popen([plink,  '--bfile', 'temp.pass_q50_dp50', '--genome', '--out', file_prefix + '.pass_q50_dp50'])
	plink_ibd.wait()



##make intervals file for indel realinment step
# relign_target_creater_knowns()

##call all methods
##sample list needed if using scalpel (proband, dad, mom, sib(if available)), but don't use it not using that program
# def call_all_exome_methods(working_dir, prefix, bam_or_fastq, sample_dict, sample_list):
def call_all_genome_methods(working_dir, prefix, bam_or_fastq, sample_dict):
	os.chdir(working_dir)
	##covert bam to fastq and then process
	if bam_or_fastq == 'bam':
		for sample in sample_dict:
			read1_fastq = sample + '.r1.fastq'
			read2_fastq = sample + '.r2.fastq'
			# convert_bam_fastq_bedtools(sample_dict[sample], read1_fastq, read2_fastq)
			convert_bam_fastq_picard(sample_dict[sample], read1_fastq, read2_fastq)
			align_with_bwa_one_at_time(sample, read1_fastq, read2_fastq)
	elif bam_or_fastq == 'messy_bam':
		for sample in sample_dict:
			read1_fastq = sample + '.r1.fastq'
			read2_fastq = sample + '.r2.fastq'
			convert_bam_fastq_bedtools(sample_dict[sample], read1_fastq, read2_fastq)
			align_with_bwa_one_at_time(sample, read1_fastq, read2_fastq)
	elif bam_or_fastq == 'fastq':
		for sample in sample_dict:
			read1_fastq = sample_dict[sample][0]
			read2_fastq = sample_dict[sample][1]
			align_with_bwa_one_at_time(sample, read1_fastq, read2_fastq)

	else:
		print 'must specify if starting with fastq or bam file'
	##variant calling
	make_list_of_bams(sample_dict, final_bam_suffix, bamlist)
	variant_calling_freebayes(bamlist, prefix)
	variant_calling_gatk_hc_no_gvcf(bamlist, prefix)
	intesect_two_vcf_files(prefix + '.gatkHC.vcf.gz', prefix + '.freebayes.vcf.gz', prefix + '.intersected_vcfs')
	##coverage
	calculate_exome_coverage(sample_dict, final_bam_suffix, exome_bed_for_coverage, prefix)
	##plink for sex and relatedness -- work in progress
	plink_relatadness_check(prefix + '.intersected_vcfs/0002.vcf', prefix)
	##scalpel -- not running at the minute as it takes so long
	#run_scalpel(sample_list, final_bam_suffix, prefix)
	##delete unwanted files
	delete_unwated_files(['*.r1.fastq', '*.r2.fastq', '*.bwa.bam', '*.bwa_sort.bam', '*.bwa_mkdup.bam', '*.metrics', '*.bwa_religned.bam', '*.bwa_religned.bai', '*.recal_data.table', '*temp*'])


def combine_fq_file(fq_dir, tissues, final_dir):
	for tissue in tissues:
		print tissue
		r1_fq = final_dir + '/' + tissue + '.r1.fq.gz'
		r2_fq = final_dir + '/' + tissue + '.r2.fq.gz'
		r1_to_combine = glob.glob(tissue + '*/*R1*')
		r2_to_combine = glob.glob(tissue + '*/*R2*')
		with open(r1_fq, 'w') as r1_fh:
			cat_files = subprocess.Popen(['cat'] + r1_to_combine, stdout=r1_fh)
			cat_files.wait()
		with open(r2_fq, 'w') as r2_fh:
			cat_files = subprocess.Popen(['cat'] + r2_to_combine, stdout=r2_fh)
			cat_files.wait()

def make_fq_dict(tissues, r1_suffix, r2_suffix):
	fq_dict = {}
	for tissue in tissues:
		r1_fq = tissue + '.r1.fq.gz'
		r2_fq = tissue + '.r2.fq.gz'
		fq_dict[tissue] = [r1_fq, r2_fq]
	return fq_dict

def load_single_vcf_into_gemini(working_dir, ped, db_name, input_vcf):
	os.chdir(working_dir)
	##dcompress, change a gatk header thing, and decompse and normalize (should be redundant)
	temp_vcf = 'temp1.vcf'
	normalized_vcf = ped + '.int.norm.vcf'
	ped_file = ped + '.ped'
	print ped_file
	with open (temp_vcf, 'w') as tvcf_fh:
		zless_vcf = subprocess.Popen(['zless', input_vcf], stdout=subprocess.PIPE)
		sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
		vt_decompose = subprocess.Popen([vt, 'decompose', '-s', '-'], stdin=sed_vcf.stdout, stdout=subprocess.PIPE)
		vt_normalize = subprocess.Popen([vt, 'normalize', '-r', fasta, '-'], stdin=vt_decompose.stdout, stdout=tvcf_fh)
		vt_normalize.wait()
	##annotate with snpeff and compress and index
	with open (normalized_vcf, 'w') as nvcf_fh:
		snpeff_vcf = subprocess.Popen(['java', '-Xmx20G', '-jar', snpeff_jar, 'GRCh37.75', '-v', '-formatEff', '-classic', temp_vcf], stdout=nvcf_fh)
		snpeff_vcf.wait()
		bgzip_vcf = subprocess.Popen(['bgzip', normalized_vcf])
		bgzip_vcf.wait()
		tabix_vcf = subprocess.Popen(['tabix', '-p', 'vcf', normalized_vcf + '.gz'])
		tabix_vcf.wait()
	##add file to gemini, requires ped file
	gemini_load = subprocess.Popen([gemini, 'load', '--cores', '10', '-t', 'snpEff', '-p', ped_file, '-v', normalized_vcf + '.gz', db_name])
	gemini_load.wait()
	##remove intermediate files
	os.remove(temp_vcf)

def query_gemini_db(query_type, gemini_db_name, outfile):
	with open(outfile, 'w') as out_fh:
		if query_type == 'first':
			query_stuff = ["chrom, start, end, ref, alt, qual, gerp_bp_score, impact, depth, qual_depth, allele_count, max_aaf_all, gt_depths.Brain, gt_ref_depths.Brain, gt_alt_depths.Brain, gt_depths.Heart, gt_ref_depths.Heart, gt_alt_depths.Heart, gt_depths.Liver, gt_ref_depths.Liver, gt_alt_depths.Liver, gt_depths.Lung, gt_ref_depths.Lung, gt_alt_depths.Lung, gt_depths.Muscle, gt_ref_depths.Muscle, gt_alt_depths.Muscle, gt_depths.Placenta, gt_ref_depths.Placenta, gt_alt_depths.Placenta from variants where filter IS NULL and in_dbsnp = 0 and in_segdup = 0 and rmsk IS NULL and type = 'snp'"]
		else:
			print 'query type %s not recognized'%query_type
		query_cmd = [gemini, 'query', '-q'] + query_stuff + ['--header', gemini_db_name]
		gemini_query = subprocess.Popen(query_cmd, stdout=out_fh)

def graphs_aaf_all_tissues(infile, outfile_prefix, cutoff_parameters):
	# 'for_r.txt'
	##make dict wuth var name and lists of, total, ref and alt counts
	var_aaf_dict, var_aaf_dict_cutoff = {}, {}
	with open(infile, 'r') as in_fh:
		line_count = 0
		for line in in_fh:
			line = line.strip('\n').split(delim)
			info = line[12:]
			line_count += 1
			if line_count == 1:
				print info, len(info)
				col_names = info
				samples = []
				for col_name in col_names:
					sample = col_name.split('.')[1]
					if sample not in samples:
						samples.append(sample)
				print samples
			else:
				depths, ref_depths, alt_depths = [], [], []
				var = '_'.join(line[:5])
				chrom = line[0]
				if chrom != 'chrX' and chrom != 'chrY' and chrom != 'chrMT' and chrom[:5] != 'chrGL':
					# print var
					i_count = 0
					for i in info:
						i_description = col_names[i_count].split('.')[0]
						i = int(i)
						##convert missing to 0
						if i == -1:
							i = 0
						# print i_description, i
						if i_description == 'gt_depths':
							depths.append(i)
						elif i_description == 'gt_ref_depths':
							ref_depths.append(i)
						elif i_description == 'gt_alt_depths':
							alt_depths.append(i)
						i_count +=1
					var_aaf_dict[var] = [depths, ref_depths, alt_depths]
					# print var, samples, depths, ref_depths, alt_depths


	##make aaf aaf column for combined sample
	combined_outfile = outfile_prefix + '.' + 'combined' + '.for_r.txt'
	passed_comb_count, passed_cutoff = 0, 0
	with open(combined_outfile, 'w') as cout_fh:
		##header
		cout_fh.write(delim.join(['variant', 'coverage', 'aaf', '\n']))
		##
		for v in var_aaf_dict:
			total = sum(var_aaf_dict[v][0])
			ref = sum(var_aaf_dict[v][1])
			alt = sum(var_aaf_dict[v][2])
			if ref + alt == total:
				aaf = float(alt)/total
				cout_fh.write(delim.join([v, str(total), str(aaf), '\n']))
				passed_comb_count += 1
				##so if coverage >75 and <300 and aaf <0.4 write to dict
				if total >= cutoff_parameters[0] and total <= cutoff_parameters[1] and aaf <= cutoff_parameters[2]:
					var_aaf_dict_cutoff[v] = [['combined'], [str(total)], [str(aaf)]]
					passed_cutoff += 1
	##make aff graphs for all indidivual samples
	s_count = -1
	for sample in samples:
		passed_count = 0
		outfile = outfile_prefix + '.' + sample + '.for_r.txt'
		outfile2 = outfile_prefix + '.' + sample + '.including_aaf0.for_r.txt'
		s_count += 1
		# print sample, s_count
		with open(outfile, 'w') as out_fh, open(outfile2, 'w') as out2_fh:
			##header
			out_fh.write(delim.join(['variant', 'coverage', 'aaf', '\n']))
			out2_fh.write(delim.join(['variant', 'coverage', 'aaf', '\n']))
			##
			for v in var_aaf_dict:
				total = var_aaf_dict[v][0][s_count]
				ref = var_aaf_dict[v][1][s_count]
				alt = var_aaf_dict[v][2][s_count]
				if ref + alt == total:
					if total != 0:
						aaf = float(alt)/total
						out2_fh.write(delim.join([v, str(total), str(aaf), '\n']))
						if v in var_aaf_dict_cutoff:
							var_aaf_dict_cutoff[v][0].append(sample)
							var_aaf_dict_cutoff[v][1].append(total)
							var_aaf_dict_cutoff[v][2].append(aaf)
						if aaf != 0:
							out_fh.write(delim.join([v, str(total), str(aaf), '\n']))
							passed_count += 1
			print sample,s_count, passed_count
	print 'we have %s total vars, of which %s are autosomal, and %s are biallelic, and %s passed the cutoff parameters'%(line_count - 1, len(var_aaf_dict), passed_comb_count, passed_cutoff)

def combine_all_aaf_files(combined_file, tissue_files, outfile, cutoff_parameters):
	passed_cutoff = 0
	header = ['variant', 'combined_coverage', 'combined_aaf']
	line_out_dict = {}
	with open(combined_file, 'r') as comb_fh:
		line_count = 0
		for line in comb_fh:
			line_count += 1
			if line_count >1:
				line = line.strip('\n').split(delim)
				var = line[0]
				cov = int(line[1])
				aaf = float(line[2])
				if cov >= cutoff_parameters[0] and cov <= cutoff_parameters[1] and aaf <= cutoff_parameters[2]:
					passed_cutoff +=1
					line_out = [var, str(cov), str(aaf)]
					# print passed_cutoff, line_out
					for tissue_file in tissue_files:
						tissue = tissue_file.split('.')[2]
						if passed_cutoff == 1:
							header.extend([tissue + '_coverage', tissue + '_aaf'])
						var_found = False
						with open(tissue_file, 'r') as tiss_fh:
							for ind_line in tiss_fh:
								ind_line = ind_line.strip('\n').split(delim)
								ind_var = ind_line[0]
								if ind_var == var:
									ind_cov = int(ind_line[1])
									ind_aaf = float(ind_line[2])
									line_out.extend([str(ind_cov), str(ind_aaf)])
									var_found = True
						if var_found == False:
							line_out.extend(['na', 'na'])
						line_out_dict[var] = line_out
	# 				print line_out
	# print header, passed_cutoff
	with open(outfile, 'w') as out_fh:
		out_fh.write(delim.join(header + ['\n']))
		for v in line_out_dict:
			out_fh.write(delim.join(line_out_dict[v] + ['\n']))

def delly_cnv_analysis(vcf_prefix, samples, bam_suffix):
	bam_files = []
	##make list of bam files
	for sample in samples:
		bam = sample + bam_suffix
		bam_files.append(bam)
	##run delly
	for cnv_type in ['DEL', 'DUP', 'INV', 'TRA']:
		out_bcf = vcf_prefix + '.delly_' + cnv_type + '.bcf'
		out_vcf = vcf_prefix + '.delly_' + cnv_type + '.vcf'
		delly_ali = subprocess.Popen([delly, 'call', '-t', cnv_type, '-x', delly_exclude_regions, '-o', out_bcf, '-g', fasta] + bam_files)
		delly_ali.wait()
		with open(out_vcf, "w") as out_fh:
			bcf_view = subprocess.Popen([bcftools, 'view', out_bcf], stdout=out_fh)
			bcf_view.wait()

##call methods

##parameters for this run
project_name = 'fetal_genomes_0816'
tissue_list = ['Brain', 'Heart', 'Liver', 'Lung', 'Muscle', 'Placenta']
fastq_dir = '/data/atimms/fetal_genomes/Bennett_WGS_20160621-31311285'
# work_dir = '/data/atimms/fetal_genomes'
work_dir = '/data/atimms/fetal_0816' ##now on 26
##make dict from tissue names
tissue_fq_dict = make_fq_dict(tissue_list, '.r1.fq.gz', '.r2.fq.gz')
gemini_db = project_name + '.db'
# in_vcf = project_name + '.intersected_vcfs/0002.vcf.gz'
in_vcf = project_name + '.gatkHC.vcf' ##use gatk
candidate_vars = project_name + '.dbsnp_passed_rpts_snps.vars.xls'
graphing_aafs_prefix = project_name + '.dbsnp_passed_rpts_snps'
##cut off for combined tissues -  min coverage, max coverage, max_aaf
comb_cutoffs = [75,300,0.4]

##combine fqs
# combine_fq_file(fastq_dir, tissue_list, work_dir)

##map, variants calling etc
# call_all_genome_methods(work_dir, project_name, 'fastq', tissue_fq_dict)

##load vcf file into gemini
# load_single_vcf_into_gemini(work_dir, project_name, gemini_db, in_vcf)

##query db
# query_gemini_db('first', gemini_db, candidate_vars)

##calculate and make files for graphing in R
# graphs_aaf_all_tissues(candidate_vars, graphing_aafs_prefix, comb_cutoffs)

##combine the files
# combine_all_aaf_files('fetal_genomes_0816.dbsnp_passed_rpts_snps.combined.for_r.txt', ['fetal_genomes_0816.dbsnp_passed_rpts_snps.Brain.including_aaf0.for_r.txt', 'fetal_genomes_0816.dbsnp_passed_rpts_snps.Heart.including_aaf0.for_r.txt', 'fetal_genomes_0816.dbsnp_passed_rpts_snps.Liver.including_aaf0.for_r.txt', 'fetal_genomes_0816.dbsnp_passed_rpts_snps.Lung.including_aaf0.for_r.txt', 'fetal_genomes_0816.dbsnp_passed_rpts_snps.Muscle.including_aaf0.for_r.txt', 'fetal_genomes_0816.dbsnp_passed_rpts_snps.Placenta.including_aaf0.for_r.txt'], 'fetal_genomes_0816.dbsnp_passed_rpts_snps.aaf_cov.comb_and_all_tissues.xls', comb_cutoffs)


##mv files to /home/atimms/fetal_0816 and graph in R

##run cnv analysis
delly_cnv_analysis(project_name, tissue_list, '.bwa_gatk.bam')




















