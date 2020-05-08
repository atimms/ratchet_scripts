#!/tools/BioBuilds-2015.04/bin/python
import sys
import subprocess
import os
import glob

##parameters
delim = '\t'

##file names etc
ref_dir = '/data/atimms/references/'
fasta = ref_dir + 'human_g1k_v37.fasta'
indels_mills = ref_dir + 'Mills_and_1000G_gold_standard.indels.b37.vcf'
indels_1000g = ref_dir + '1000G_phase1.indels.b37.vcf'
rtc_intervals = ref_dir + 'Mills_1000G.intervals'
dbsnp = ref_dir + 'dbsnp_138.b37.vcf'
exome_capture_bed = ref_dir + 'dobyns_exome.in_any_padded_target.1015.bed'
final_bam_suffix = '.bwa_gatk.bam'
bamlist = 'bams.list'
exome_bed_for_coverage = ref_dir + 'dobyns_exome.in_all_targets.1015.bed'


##programs
samtools = '/home/atimms/programs/samtools-1.3/samtools'
bedtools = '/home/atimms/programs/bedtools2-master/bin/bedtools' ##version 2.25
picard = '/home/atimms/programs/picard-tools-1.139/picard.jar'
bwa = '/home/atimms/programs/bwa-0.7.12/bwa'
gatk = '/home/atimms/programs/GenomeAnalysisTK.jar'
bcftools = '/home/atimms/programs/bcftools-1.2/bcftools'
freebayes = '/home/atimms/programs/freebayes/bin/freebayes'
bgzip = '/home/atimms/programs/htslib-1.3/bgzip'
plink = '/home/atimms/programs/plink'
vt = '/data/atimms/gemini/vt/vt'


def relign_target_creater_knowns():
	gatk_rtc = subprocess.Popen(['java', '-Xmx3g', '-jar', gatk, '-T', 'RealignerTargetCreator', '-R', fasta, '-nt', '15', '-o', rtc_intervals, '-known', indels_mills, '-known', indels_1000g])
	gatk_rtc.wait()

def convert_bam_fastq_bedtools(bamfile, r1_fastq, r2_fastq):
	read_sorted_bam = bamfile[:-4] + '.n_sorted.temp.bam'
	st_n_sort = subprocess.Popen([samtools, 'sort', '-nO', 'bam', '-o', read_sorted_bam, '-@', '10', '-m', '8G', '-T', 'tempy', bamfile])
	st_n_sort.wait()
	bam_fq = subprocess.Popen([bedtools, 'bamtofastq', '-i', read_sorted_bam, '-fq', r1_fastq, '-fq2', r2_fastq])
	bam_fq.wait()

def convert_bam_fastq_picard(bamfile, r1_fastq, r2_fastq):
	picard_sam_fq = subprocess.Popen(['java', '-Xmx50g', '-jar', picard, 'SamToFastq', 'INPUT=' + bamfile, 'FASTQ=' + r1_fastq, 'SECOND_END_FASTQ=' + r2_fastq, 'VALIDATION_STRINGENCY=SILENT'])
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
	bwa_pe = subprocess.Popen([bwa, 'mem', '-M', '-t', '25', '-R', rg, fasta, r1_fq, r2_fq], stdout=subprocess.PIPE)
	st_sam_bam_pe = subprocess.Popen([samtools, 'view', '-b', '-@', '5', '-o', post_bwa_bam, '-'], stdin=bwa_pe.stdout)
	st_sam_bam_pe.wait()
	##sort bam
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', sample, '-o', sort_bam, '-@', '8', '-m', '8G', post_bwa_bam])
	st_sort_pe.wait()
	##mark duplicates
	picard_md = subprocess.Popen(['java', '-Xmx64g', '-jar', picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=LENIENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + mkdup_bam, 'METRICS_FILE=' + sample + '.metrics'])
	picard_md.wait()
	##realign around indels
	gatk_ir = subprocess.Popen(['java', '-Xmx64g', '-jar', gatk, '-T', 'IndelRealigner', '-R', fasta, '-I', mkdup_bam, '-L', exome_capture_bed, '-targetIntervals', rtc_intervals, '--consensusDeterminationModel', 'KNOWNS_ONLY', '-LOD', '0.4', '-o', realigned_bam, '-known', indels_mills, '-known', indels_1000g])
	gatk_ir.wait()
	##bqsr
	gatk_br = subprocess.Popen(['java', '-Xmx64g', '-jar', gatk, '-T', 'BaseRecalibrator', '-nct', '8', '-R', fasta, '-I', realigned_bam, '-L', exome_capture_bed, '-knownSites', indels_mills, '-knownSites', dbsnp, '-knownSites', indels_1000g, '-o', sample + '.recal_data.table'])
	gatk_br.wait()
	gatk_pr = subprocess.Popen(['java', '-Xmx64g', '-jar', gatk, '-T', 'PrintReads', '-nct', '8', '-R', fasta, '-I', realigned_bam, '-L', exome_capture_bed, '-BQSR', sample + '.recal_data.table', '-o', gatk_bam])
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
	freebayes_run = subprocess.Popen([freebayes,'-f' , fasta, '-L', bamlist, '-t', exome_capture_bed], stdout=vcf_fh)
	freebayes_run.wait()
	vcf_fh.close()
	bgzip_run = subprocess.Popen([bgzip, vcf_temp1])
	bgzip_run.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines
	# bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	# bcf_norm1.wait()
	# bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	# bcf_norm2.wait()
	bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-O', 'z', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm2.wait()

	bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-o', final_vcf, vcf_temp2 ])
	bcf_norm1.wait()


	bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
	bcf_index.wait()


##call freebayes in all sample listed in list of bamfiles
def variant_calling_freebayes_with_vt(bamlist, name_prefix):
	vcf_temp1 = 'temp_fb1.vcf'
	# vcf_temp2 = 'temp_fb2.vcf.gz'
	final_vcf = name_prefix + '.freebayes.vcf'
	vcf_fh = open(vcf_temp1, 'w')
	freebayes_run = subprocess.Popen([freebayes,'-f' , fasta, '-L', bamlist, '-t', exome_capture_bed], stdout=vcf_fh)
	freebayes_run.wait()
	vcf_fh.close()
	bgzip_run = subprocess.Popen([bgzip, vcf_temp1])
	bgzip_run.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines
	# bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	# bcf_norm1.wait()
	# bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	# bcf_norm2.wait()

	bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
	bcf_index.wait()
	##dcompress, change a gatk header thing, and decompse and normalize (should be redundant)
	with open (final_vcf, 'w') as tvcf_fh:
		zless_vcf = subprocess.Popen(['zless', vcf_temp1 + '.gz'], stdout=subprocess.PIPE)
		sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
		vt_decompose = subprocess.Popen([vt, 'decompose', '-s', '-'], stdin=sed_vcf.stdout, stdout=subprocess.PIPE)
		vt_normalize = subprocess.Popen([vt, 'normalize', '-r', fasta, '-'], stdin=vt_decompose.stdout, stdout=tvcf_fh)
		vt_normalize.wait()
	bgzip_run = subprocess.Popen([bgzip, final_vcf])
	bgzip_run.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', final_vcf + '.gz'])
	bcf_index.wait()

##not used!!
def variant_calling_gatk_hc(sample_dict, name_prefix):
	vcf_temp0 = 'temp_gatk0.vcf'
	vcf_raw_snps = 'temp_raw_snps.vcf'
	vcf_filtered_snps = 'temp_filtered_snps.vcf'
	vcf_raw_indels = 'temp_raw_indels.vcf'
	vcf_filtered_indels = 'temp_filtered_indels.vcf'
	vcf_temp1 = 'temp_gatk1.vcf'
	vcf_temp2 = 'temp_gatk2.vcf.gz'
	final_vcf = name_prefix + '.gatkHC.vcf.gz'
	gvcf_files = []
	##get g.vcfs for each sample (add -L)
	for sample in sample_dict:
		gatk_hc = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fasta, '-I', sample + final_bam_suffix,'-L', exome_capture_bed, '--emitRefConfidence', 'GVCF', '--variant_index_type', 'LINEAR', '--variant_index_parameter', '128000', '-o', sample + '.g.vcf'])
		gatk_hc.wait()
		gvcf = ['-V', sample + '.g.vcf']
		gvcf_files.extend(gvcf)
	print gvcf_files
	##genotype g.vcfs(add -L)
	command = ['java', '-Xmx20g', '-jar', gatk, '-T', 'GenotypeGVCFs', '-R', fasta, '-nt', '15','-L', exome_capture_bed] + gvcf_files + ['-o', vcf_temp0]
	gatk_gg = subprocess.Popen(command)
	gatk_gg.wait()
	##split data in SNPs and indels and apply manual variant filtering (add -L)
	snp_cut = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'SelectVariants', '-R',fasta, '-nt', '15', '-V', vcf_temp0, '-o', vcf_raw_snps, '-selectType', 'SNP'])
	snp_cut.wait()
	snp_vf = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'VariantFiltration', '-R',fasta, '-V', vcf_raw_snps, '-o', vcf_filtered_snps, '--filterExpression ', "QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0", '--filterName', "manual_filter"])
	snp_vf.wait()
	indel_cut = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'SelectVariants', '-R',fasta, '-nt', '15', '-V', vcf_temp0, '-o', vcf_raw_indels, '-selectType', 'INDEL'])
	indel_cut.wait()
	indel_vf = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'VariantFiltration', '-R',fasta, '-V', vcf_raw_indels, '-o', vcf_filtered_indels, '--filterExpression ', "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0", '--filterName', "manual_filter"])
	indel_vf.wait()
	##combine filtered snps and indels
	combine_var = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineVariants', '-R',fasta, '-nt', '15', '--variant', vcf_filtered_snps,'--variant', vcf_filtered_indels,'-o', vcf_temp1, '-genotypeMergeOptions', 'UNSORTED'])
	combine_var.wait()
	bgzip_cmd = subprocess.Popen([bgzip, vcf_temp1])
	bgzip_cmd.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines, and left normalize indels
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
	gatk_hc = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fasta, '-I', bamlist,'-L', exome_capture_bed, '-stand_call_conf', '30', '-stand_emit_conf', '10', '-o', vcf_temp0])
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


def run_scalpel(samples, bam_suffix, prefix):
	print samples
	if len(samples) == 3:
		scalpel_trio_disco = subprocess.Popen([scalpel_discovery, '--denovo', '--dad', samples[1] + bam_suffix, '--mom', samples[2] + bam_suffix, '--aff', samples[0] + bam_suffix, '--sib', samples[0] + bam_suffix, '--bed', exome_capture_bed, '--ref', fasta, '--numprocs', '10', '--two-pass'])
		scalpel_trio_disco.wait()
		# scalpel_trio_export = subprocess.Popen([scalpel_discovery, '--denovo', '--bed', exome_capture_bed, '--ref', fasta])
		# scalpel_trio_export.wait()

##make intervals file for indel realinment step
# relign_target_creater_knowns()

##call all methods
##sample list needed if using scalpel (proband, dad, mom, sib(if available)), but don't use it not using that program
# def call_all_exome_methods(working_dir, prefix, bam_or_fastq, sample_dict, sample_list):
def call_all_exome_methods(working_dir, prefix, bam_or_fastq, sample_dict):
	os.chdir(working_dir)
	##covert bam to fastq and then process
	# '''
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
			##single fq
			# align_with_bwa_one_at_time_single_end(sample,sample_dict[sample])
	else:
		print 'must specify if starting with fastq or bam file'
	# '''
	##variant calling
	make_list_of_bams(sample_dict, final_bam_suffix, bamlist)
	variant_calling_freebayes_with_vt(bamlist, prefix)
	variant_calling_gatk_hc_no_gvcf(bamlist, prefix)
	intesect_two_vcf_files(prefix + '.gatkHC.vcf.gz', prefix + '.freebayes.vcf.gz', prefix + '.intersected_vcfs')
	##coverage
	calculate_exome_coverage(sample_dict, final_bam_suffix, exome_bed_for_coverage, prefix)
	##plink for sex and relatedness -- work in progress
	plink_relatadness_check(prefix + '.intersected_vcfs/0002.vcf', prefix)
	##scalpel -- not running at the minute as it takes so long
	# #run_scalpel(sample_list, final_bam_suffix, prefix)
	##delete unwanted files
	delete_unwated_files(['*.r1.fastq', '*.r2.fastq', '*.bwa.bam', '*.bwa_sort.bam', '*.metrics', '*.bwa_religned.bam', '*.bwa_religned.bai', '*.recal_data.table', '*temp*'])

def call_exome_methods_freebayes_plus(working_dir, prefix, bam_or_fastq, sample_dict):
	os.chdir(working_dir)
	# '''
	##variant calling
	make_list_of_bams(sample_dict, final_bam_suffix, bamlist)
	variant_calling_freebayes_with_vt(bamlist, prefix)
	# variant_calling_gatk_hc_no_gvcf(bamlist, prefix)
	intesect_two_vcf_files(prefix + '.gatkHC.vcf.gz', prefix + '.freebayes.vcf.gz', prefix + '.intersected_vcfs')
	##coverage
	calculate_exome_coverage(sample_dict, final_bam_suffix, exome_bed_for_coverage, prefix)
	##plink for sex and relatedness -- work in progress
	plink_relatadness_check(prefix + '.intersected_vcfs/0002.vcf', prefix)
	##scalpel -- not running at the minute as it takes so long
	# #run_scalpel(sample_list, final_bam_suffix, prefix)
	##delete unwanted files
	delete_unwated_files(['*.r1.fastq', '*.r2.fastq', '*.bwa.bam', '*.bwa_sort.bam', '*.metrics', '*.bwa_religned.bam', '*.bwa_religned.bai', '*.recal_data.table', '*temp*'])


def call_exome_methods_make_bam(working_dir, prefix, bam_or_fastq, sample_dict):
	os.chdir(working_dir)
	##covert bam to fastq and then process
	# '''
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
			##single fq
			# align_with_bwa_one_at_time_single_end(sample,sample_dict[sample])
	else:
		print 'must specify if starting with fastq or bam file'
	delete_unwated_files(['*.r1.fastq', '*.r2.fastq', '*.bwa.bam', '*.bwa_sort.bam', '*.metrics', '*.bwa_religned.bam', '*.bwa_religned.bai', '*.recal_data.table', '*temp*'])

def variant_calling_plus(working_dir, prefix, bam_or_fastq, sample_dict):
	os.chdir(working_dir)
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
	# #run_scalpel(sample_list, final_bam_suffix, prefix)
	##delete unwanted files
	delete_unwated_files(['*.r1.fastq', '*.r2.fastq', '*.bwa.bam', '*.bwa_sort.bam', '*.metrics', '*.bwa_religned.bam', '*.bwa_religned.bai', '*.recal_data.table', '*temp*'])


##run methods
##bam, trio and single
#call_all_exome_methods(working_dir, '', 'bam', {'':'', '':'', '':''})
#call_all_exome_methods(working_dir, '', 'bam', {'':''})


##all methods, but freebayes didn't work so redo..
##batch1 -- ready
# working_dir = '/data/atimms/microtia_exomes/batch1'
##issue with bam
# call_all_exome_methods(working_dir, '1010001', 'bam', {'101000101':'60017.bam'})
##run
# call_all_exome_methods(working_dir, '1010002', 'bam', {'101000201':'60018.bam', '101000202':'60019.bam', '101000203':'60020.bam'})
# call_all_exome_methods(working_dir, '1010003', 'bam', {'101000301':'60021.bam', '101000302':'60022.bam'})
# call_all_exome_methods(working_dir, '1010004', 'bam', {'101000401':'60023.bam', '101000402':'60024.bam', '101000403':'60025.bam'})
# call_all_exome_methods(working_dir, '1010005', 'bam', {'101000501':'60026.bam', '101000502':'60027.bam', '101000503':'60028.bam'})
# call_all_exome_methods(working_dir, '1010006', 'bam', {'101000601':'60029.bam', '101000602':'60030.bam'})
# call_all_exome_methods(working_dir, '1010007', 'bam', {'101000701':'60031.bam', '101000702':'60032.bam', '101000703':'101000703.94638.bam'})
# call_all_exome_methods(working_dir, '1030001', 'bam', {'103000101':'60033.bam', '103000102':'60034.bam', '103000103':'60035.bam'})

##batch2 -- ready
# working_dir = '/data/atimms/microtia_exomes/batch2'
# call_all_exome_methods(working_dir, '1010008', 'bam', {'101000801':'101000801.94596.bam', '101000802':'101000802.94597.bam', '101000803':'101000803.94598.bam'})
# call_all_exome_methods(working_dir, '1010010', 'bam', {'101001001':'101001001.94599.bam', '101001002':'101001002.94600.bam', '101001003':'101001003.94601.bam'})
# call_all_exome_methods(working_dir, '1010013', 'bam', {'101001301':'101001301.94605.bam', '101001302':'101001302.94606.bam', '101001303':'101001303.94607.bam'})
# call_all_exome_methods(working_dir, '1010019', 'bam', {'101001901':'101001901.94611.bam', '101001902':'101001902.94612.bam', '101001903':'101001903.94613.bam'})
# call_all_exome_methods(working_dir, '1010020', 'bam', {'101002001':'101002001.94614.bam', '101002002':'101002002.94615.bam', '101002003':'101002003.94616.bam'})
# call_all_exome_methods(working_dir, '1010021', 'bam', {'101002101':'101002101.94617.bam', '101002102':'101002102.94618.bam', '101002103':'101002103.94619.bam'})
# call_all_exome_methods(working_dir, '1010022', 'bam', {'101002201':'101002201.94620.bam', '101002202':'101002202.94621.bam', '101002203':'101002203.94622.bam'})
# call_all_exome_methods(working_dir, '1010023', 'bam', {'101002301':'101002301.94623.bam', '101002302':'101002302.94625.bam', '101002303':'101002303.94624.bam'})
# call_all_exome_methods(working_dir, '1010024', 'bam', {'101002401':'101002401.94627.bam', '101002402':'101002402.94626.bam', '101002403':'101002403.94628.bam'})
# call_all_exome_methods(working_dir, '1010025', 'bam', {'101002501':'101002501.94629.bam', '101002502':'101002502.94630.bam', '101002503':'101002503.94631.bam'})
# call_all_exome_methods(working_dir, '1010026', 'bam', {'101002601':'101002601.94632.bam', '101002602':'101002602.94633.bam', '101002603':'101002603.94634.bam'})
# call_all_exome_methods(working_dir, '1010028', 'bam', {'101002801':'101002801.94635.bam', '101002802':'101002802.94636.bam', '101002803':'101002803.94637.bam'})

##batch3 -- ready
# working_dir = '/data/atimms/microtia_exomes/batch3'
# call_all_exome_methods(working_dir, '1030002', 'bam', {'103000201':'103000201.138424.bam', '103000202':'103000202.138425.bam', '103000203':'103000203.138426.bam'})
# call_all_exome_methods(working_dir, '1030006', 'bam', {'103000601':'103000601.138427.bam', '103000602':'103000602.138428.bam', '103000603':'103000603.138429.bam'})
# call_all_exome_methods(working_dir, '1030007', 'bam', {'103000701':'103000701.138430.bam', '103000702':'103000702.138431.bam', '103000703':'103000703.138432.bam'})
# call_all_exome_methods(working_dir, '1060004', 'bam', {'106000401':'106000401.138433.bam', '106000402':'106000402.138434.bam', '106000403':'106000403.138435.bam'})
# call_all_exome_methods(working_dir, '1060015', 'bam', {'106001501':'106001501.138436.bam', '106001502':'106001502.138437.bam', '106001503':'106001503.138438.bam'})
# call_all_exome_methods(working_dir, '1070001', 'bam', {'107000101':'107000101.138439.bam', '107000102':'107000102.138440.bam', '107000103':'107000103.138441.bam'})
# call_all_exome_methods(working_dir, '1070002', 'bam', {'107000201':'107000201.138442.bam', '107000202':'107000202.138443.bam', '107000203':'107000203.138444.bam'})
# call_all_exome_methods(working_dir, '1070003', 'bam', {'107000301':'107000301.138445.bam', '107000302':'107000302.138446.bam', '107000303':'107000303.138447.bam'})
# call_all_exome_methods(working_dir, '1070011', 'bam', {'107001101':'107001101.138448.bam', '107001102':'107001102.138449.bam', '107001103':'107001103.138450.bam'})

##batch1 -- ready
# working_dir = '/data/atimms/microtia_exomes/batch1'
# call_exome_methods_freebayes_plus(working_dir, '1010002', 'bam', {'101000201':'60018.bam', '101000202':'60019.bam', '101000203':'60020.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1010003', 'bam', {'101000301':'60021.bam', '101000302':'60022.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1010004', 'bam', {'101000401':'60023.bam', '101000402':'60024.bam', '101000403':'60025.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1010005', 'bam', {'101000501':'60026.bam', '101000502':'60027.bam', '101000503':'60028.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1010006', 'bam', {'101000601':'60029.bam', '101000602':'60030.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1010007', 'bam', {'101000701':'60031.bam', '101000702':'60032.bam', '101000703':'101000703.94638.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1030001', 'bam', {'103000101':'60033.bam', '103000102':'60034.bam', '103000103':'60035.bam'})

##batch2 -- ready
# working_dir = '/data/atimms/microtia_exomes/batch2'
# call_exome_methods_freebayes_plus(working_dir, '1010008', 'bam', {'101000801':'101000801.94596.bam', '101000802':'101000802.94597.bam', '101000803':'101000803.94598.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1010010', 'bam', {'101001001':'101001001.94599.bam', '101001002':'101001002.94600.bam', '101001003':'101001003.94601.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1010013', 'bam', {'101001301':'101001301.94605.bam', '101001302':'101001302.94606.bam', '101001303':'101001303.94607.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1010019', 'bam', {'101001901':'101001901.94611.bam', '101001902':'101001902.94612.bam', '101001903':'101001903.94613.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1010020', 'bam', {'101002001':'101002001.94614.bam', '101002002':'101002002.94615.bam', '101002003':'101002003.94616.bam'})
##redo all methods, proband missing bam files, so repeat all
# call_exome_methods_freebayes_plus(working_dir, '1010021', 'bam', {'101002101':'101002101.94617.bam', '101002102':'101002102.94618.bam', '101002103':'101002103.94619.bam'})
# call_all_exome_methods(working_dir, '1010021', 'bam', {'101002101':'101002101.94617.bam', '101002102':'101002102.94618.bam', '101002103':'101002103.94619.bam'})


# call_exome_methods_freebayes_plus(working_dir, '1010022', 'bam', {'101002201':'101002201.94620.bam', '101002202':'101002202.94621.bam', '101002203':'101002203.94622.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1010023', 'bam', {'101002301':'101002301.94623.bam', '101002302':'101002302.94625.bam', '101002303':'101002303.94624.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1010024', 'bam', {'101002401':'101002401.94627.bam', '101002402':'101002402.94626.bam', '101002403':'101002403.94628.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1010025', 'bam', {'101002501':'101002501.94629.bam', '101002502':'101002502.94630.bam', '101002503':'101002503.94631.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1010026', 'bam', {'101002601':'101002601.94632.bam', '101002602':'101002602.94633.bam', '101002603':'101002603.94634.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1010028', 'bam', {'101002801':'101002801.94635.bam', '101002802':'101002802.94636.bam', '101002803':'101002803.94637.bam'})

##batch3 -- ready
# working_dir = '/data/atimms/microtia_exomes/batch3'
# call_exome_methods_freebayes_plus(working_dir, '1030002', 'bam', {'103000201':'103000201.138424.bam', '103000202':'103000202.138425.bam', '103000203':'103000203.138426.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1030006', 'bam', {'103000601':'103000601.138427.bam', '103000602':'103000602.138428.bam', '103000603':'103000603.138429.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1030007', 'bam', {'103000701':'103000701.138430.bam', '103000702':'103000702.138431.bam', '103000703':'103000703.138432.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1060004', 'bam', {'106000401':'106000401.138433.bam', '106000402':'106000402.138434.bam', '106000403':'106000403.138435.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1060015', 'bam', {'106001501':'106001501.138436.bam', '106001502':'106001502.138437.bam', '106001503':'106001503.138438.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1070001', 'bam', {'107000101':'107000101.138439.bam', '107000102':'107000102.138440.bam', '107000103':'107000103.138441.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1070002', 'bam', {'107000201':'107000201.138442.bam', '107000202':'107000202.138443.bam', '107000203':'107000203.138444.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1070003', 'bam', {'107000301':'107000301.138445.bam', '107000302':'107000302.138446.bam', '107000303':'107000303.138447.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1070011', 'bam', {'107001101':'107001101.138448.bam', '107001102':'107001102.138449.bam', '107001103':'107001103.138450.bam'})

##batch3 -- ready
##deleted bams for so redo and also rpt 1030007 as father is not related
# working_dir = '/data/atimms/microtia_exomes/batch3'
# call_exome_methods_make_bam(working_dir, '1030002', 'bam', {'103000201':'103000201.138424.bam', '103000202':'103000202.138425.bam', '103000203':'103000203.138426.bam'})
# call_exome_methods_make_bam(working_dir, '1030006', 'bam', {'103000601':'103000601.138427.bam', '103000602':'103000602.138428.bam', '103000603':'103000603.138429.bam'})
# call_all_exome_methods(working_dir, '1030007', 'bam', {'103000701':'103000701.138430.bam', '103000702':'103000702.138431.bam'})
# call_exome_methods_make_bam(working_dir, '1060004', 'bam', {'106000401':'106000401.138433.bam', '106000402':'106000402.138434.bam', '106000403':'106000403.138435.bam'})
# call_exome_methods_make_bam(working_dir, '1060015', 'bam', {'106001501':'106001501.138436.bam', '106001502':'106001502.138437.bam', '106001503':'106001503.138438.bam'})
# call_exome_methods_make_bam(working_dir, '1070001', 'bam', {'107000101':'107000101.138439.bam', '107000102':'107000102.138440.bam', '107000103':'107000103.138441.bam'})
# call_exome_methods_make_bam(working_dir, '1070002', 'bam', {'107000201':'107000201.138442.bam', '107000202':'107000202.138443.bam', '107000203':'107000203.138444.bam'})
# call_exome_methods_make_bam(working_dir, '1070003', 'bam', {'107000301':'107000301.138445.bam', '107000302':'107000302.138446.bam', '107000303':'107000303.138447.bam'})
# call_exome_methods_make_bam(working_dir, '1070011', 'bam', {'107001101':'107001101.138448.bam', '107001102':'107001102.138449.bam', '107001103':'107001103.138450.bam'})
# call_exome_methods_freebayes_plus(working_dir, '1030007', 'bam', {'103000701':'103000701.138430.bam', '103000702':'103000702.138431.bam'})

##batch4 -- ready
# working_dir = '/data/atimms/microtia_exomes/batch4'
# call_all_exome_methods(working_dir, '1010031', 'fastq', {'101003101':['101003101_GTCGTAGA_R1.fq.gz', '101003101_GTCGTAGA_R2.fq.gz'], '101003102':['101003102_GTCTGTCA_R1.fq.gz', '101003102_GTCTGTCA_R2.fq.gz'], '101003103':['101003103_GTGTTCTA_R1.fq.gz', '101003103_GTGTTCTA_R2.fq.gz']})
# call_all_exome_methods(working_dir, '1010032', 'fastq', {'101003201':['101003201_TAGGATGA_R1.fq.gz', '101003201_TAGGATGA_R2.fq.gz'], '101003202':['101003202_TATCAGCA_R1.fq.gz', '101003202_TATCAGCA_R2.fq.gz'], '101003203':['101003203_TCCGTCTA_R1.fq.gz', '101003203_TCCGTCTA_R2.fq.gz']})
# call_all_exome_methods(working_dir, '1030013', 'fastq', {'103001301':['103001301_TCTTCACA_R1.fq.gz', '103001301_TCTTCACA_R2.fq.gz'], '103001302':['103001302_TGAAGAGA_R1.fq.gz', '103001302_TGAAGAGA_R2.fq.gz'], '103001303':['103001303_TGGAACAA_R1.fq.gz', '103001303_TGGAACAA_R2.fq.gz']})
# call_all_exome_methods(working_dir, '1030014', 'fastq', {'103001401':['103001401_TGGCTTCA_R1.fq.gz', '103001401_TGGCTTCA_R2.fq.gz'], '103001402':['103001402_TGGTGGTA_R1.fq.gz', '103001402_TGGTGGTA_R2.fq.gz'], '103001403':['103001403_TTCACGCA_R1.fq.gz', '103001403_TTCACGCA_R2.fq.gz']})
# call_all_exome_methods(working_dir, '1030015', 'fastq', {'103001501':['103001501_AACTCACC_R1.fq.gz', '103001501_AACTCACC_R2.fq.gz'], '103001502':['103001502_AAGAGATC_R1.fq.gz', '103001502_AAGAGATC_R2.fq.gz'], '103001503':['103001503_AAGGACAC_R1.fq.gz', '103001503_AAGGACAC_R2.fq.gz']})
# call_all_exome_methods(working_dir, '1060018', 'fastq', {'106001801':['106001801_AATCCGTC_R1.fq.gz', '106001801_AATCCGTC_R2.fq.gz'], '106001806a':['106001806a_AATGTTGC_R1.fq.gz', '106001806a_AATGTTGC_R2.fq.gz'], '106001806b':['106001806b_ACACGACC_R1.fq.gz', '106001806b_ACACGACC_R2.fq.gz'], '106001809':['106001809_ACAGATTC_R1.fq.gz', '106001809_ACAGATTC_R2.fq.gz']})
# call_all_exome_methods(working_dir, '1060020', 'fastq', {'106002001':['106002001_AGATGTAC_R1.fq.gz', '106002001_AGATGTAC_R2.fq.gz'], '106002002':['106002002_AGCACCTC_R1.fq.gz', '106002002_AGCACCTC_R2.fq.gz'], '106002003':['106002003_AGCCATGC_R1.fq.gz', '106002003_AGCCATGC_R2.fq.gz']})
# call_all_exome_methods(working_dir, '1070012', 'fastq', {'107001201':['107001201_AGGCTAAC_R1.fq.gz', '107001201_AGGCTAAC_R2.fq.gz'], '107001202':['107001202_ATAGCGAC_R1.fq.gz', '107001202_ATAGCGAC_R2.fq.gz'], '107001203':['107001203_ATCATTCC_R1.fq.gz', '107001203_ATCATTCC_R2.fq.gz']})
# call_all_exome_methods(working_dir, '1070013', 'fastq', {'107001301':['107001301_ATTGGCTC_R1.fq.gz', '107001301_ATTGGCTC_R2.fq.gz'], '107001302':['107001302_CAAGGAGC_R1.fq.gz', '107001302_CAAGGAGC_R2.fq.gz'], '107001303':['107001303_CACCTTAC_R1.fq.gz', '107001303_CACCTTAC_R2.fq.gz']})
# call_all_exome_methods(working_dir, '1070015', 'fastq', {'107001501':['107001501_CCATCCTC_R1.fq.gz', '107001501_CCATCCTC_R2.fq.gz'], '107001502':['107001502_CCGACAAC_R1.fq.gz', '107001502_CCGACAAC_R2.fq.gz'], '107001503':['107001503_CCTAATCC_R1.fq.gz', '107001503_CCTAATCC_R2.fq.gz']})

##repeats becaus bams got messed up
# working_dir = '/data/atimms/microtia_exomes/repeats'
# call_exome_methods_make_bam(working_dir, '', 'fastq', {'101003103':['101003103_GTGTTCTA_R1.fq.gz', '101003103_GTGTTCTA_R2.fq.gz'], '101003202':['101003202_TATCAGCA_R1.fq.gz', '101003202_TATCAGCA_R2.fq.gz'], '103001301':['103001301_TCTTCACA_R1.fq.gz', '103001301_TCTTCACA_R2.fq.gz'], '103001302':['103001302_TGAAGAGA_R1.fq.gz', '103001302_TGAAGAGA_R2.fq.gz'], '103001303':['103001303_TGGAACAA_R1.fq.gz', '103001303_TGGAACAA_R2.fq.gz'], '107001203':['107001203_ATCATTCC_R1.fq.gz', '107001203_ATCATTCC_R2.fq.gz'], '107001301':['107001301_ATTGGCTC_R1.fq.gz', '107001301_ATTGGCTC_R2.fq.gz'], '107001501':['107001501_CCATCCTC_R1.fq.gz', '107001501_CCATCCTC_R2.fq.gz'], '107001503':['107001503_CCTAATCC_R1.fq.gz', '107001503_CCTAATCC_R2.fq.gz']})
# call_exome_methods_make_bam(working_dir, '', 'bam', {'101002101':'101002101.94617.bam', '101002103':'101002103.94619.bam', '103000201':'103000201.138424.bam', '103000603':'103000603.138429.bam', '103000701':'103000701.138430.bam', '106000401':'106000401.138433.bam', '106001501':'106001501.138436.bam', '106001502':'106001502.138437.bam', '107000101':'107000101.138439.bam', '107000103':'107000103.138441.bam', '107000303':'107000303.138447.bam', '107001101':'107001101.138448.bam', '107001102':'107001102.138449.bam', '107001103':'107001103.138450.bam'})
# call_exome_methods_make_bam(working_dir, '', 'bam', {'107000302':'107000302.138446.bam'})

##batch5 -- ready
working_dir = '/data/atimms/microtia_exomes/batch5_0117'
# call_all_exome_methods(working_dir, '1030011', 'fastq', {'103001101':['103001101_S26_L005_R1_001.fastq.gz', '103001101_S26_L005_R2_001.fastq.gz'], '103001102':['103001102_S27_L005_R1_001.fastq.gz', '103001102_S27_L005_R2_001.fastq.gz'], '103001103':['103001103_S28_L005_R1_001.fastq.gz', '103001103_S28_L005_R2_001.fastq.gz']})
# call_all_exome_methods(working_dir, '1060023', 'fastq', {'106002301':['106002301_S29_L005_R1_001.fastq.gz', '106002301_S29_L005_R2_001.fastq.gz'], '106002302':['106002302_S30_L005_R1_001.fastq.gz', '106002302_S30_L005_R2_001.fastq.gz'], '106002303':['106002303_S31_L005_R1_001.fastq.gz', '106002303_S31_L005_R2_001.fastq.gz']})
# call_all_exome_methods(working_dir, '1060029', 'fastq', {'106002901':['106002901_S32_L005_R1_001.fastq.gz', '106002901_S32_L005_R2_001.fastq.gz'], '106002902':['106002902_S33_L005_R1_001.fastq.gz', '106002902_S33_L005_R2_001.fastq.gz'], '106002903':['106002903_S34_L005_R1_001.fastq.gz', '106002903_S34_L005_R2_001.fastq.gz'], '106002907':['106002907_S35_L005_R1_001.fastq.gz', '106002907_S35_L005_R2_001.fastq.gz']})
# call_all_exome_methods(working_dir, '4010002', 'fastq', {'401000201':['401000201_S36_L005_R1_001.fastq.gz', '401000201_S36_L005_R2_001.fastq.gz'], '401000202':['401000202_S37_L005_R1_001.fastq.gz', '401000202_S37_L005_R2_001.fastq.gz'], '401000203':['401000203_S38_L005_R1_001.fastq.gz', '401000203_S38_L005_R2_001.fastq.gz']})
# call_all_exome_methods(working_dir, 'CFM-MOS-03', 'fastq', {'CFM-MOS-03-01':['CFM-MOS-03-01_S42_L005_R1_001.fastq.gz', 'CFM-MOS-03-01_S42_L005_R2_001.fastq.gz'], 'CFM-MOS-03-02':['CFM-MOS-03-02_S45_L005_R1_001.fastq.gz', 'CFM-MOS-03-02_S45_L005_R2_001.fastq.gz'], 'CFM-MOS-03-04':['CFM-MOS-03-04_S46_L005_R1_001.fastq.gz', 'CFM-MOS-03-04_S46_L005_R2_001.fastq.gz'], 'CFM-MOS-03-11':['CFM-MOS-03-11_S44_L005_R1_001.fastq.gz', 'CFM-MOS-03-11_S44_L005_R2_001.fastq.gz'], 'CFM-MOS-03-21':['CFM-MOS-03-21_S43_L005_R1_001.fastq.gz', 'CFM-MOS-03-21_S43_L005_R2_001.fastq.gz']})
##make to one large pedigree - convert 1060015 files to fastq and analyze with 1060012 and 1060014 (cousins, so have different ids)
# ##call_all_exome_methods(working_dir, '1060015', 'bam', {'106001501':'106001501.138436.bam', '106001502':'106001502.138437.bam', '106001503':'106001503.138438.bam'}) ##b3 data
# ##call_all_exome_methods(working_dir, '1060012', 'fastq', {'106001201':['106001201_S39_L005_R1_001.fastq.gz', '106001201_S39_L005_R2_001.fastq.gz']}) #new 
# ##call_all_exome_methods(working_dir, '1060014', 'fastq', {'106001401':['106001401_S40_L005_R1_001.fastq.gz', '106001401_S40_L005_R2_001.fastq.gz']}) #new
# call_all_exome_methods(working_dir, '1060012', 'fastq', {'106001201':['106001201_S39_L005_R1_001.fastq.gz', '106001201_S39_L005_R2_001.fastq.gz'], '106001401':['106001401_S40_L005_R1_001.fastq.gz', '106001401_S40_L005_R2_001.fastq.gz'], '106001501':['106001501.r1.fq.gz', '106001501.r2.fq.gz'], '106001502':['106001502.r1.fq.gz', '106001502.r2.fq.gz'], '106001503':['106001503.r1.fq.gz', '106001503.r2.fq.gz']})
##106001401 not related, remove and repeat
call_all_exome_methods(working_dir, '1060012', 'fastq', {'106001201':['106001201_S39_L005_R1_001.fastq.gz', '106001201_S39_L005_R2_001.fastq.gz'], '106001501':['106001501.r1.fq.gz', '106001501.r2.fq.gz'], '106001502':['106001502.r1.fq.gz', '106001502.r2.fq.gz'], '106001503':['106001503.r1.fq.gz', '106001503.r2.fq.gz']})

##106001805 add to previous files from batch4 (had to recombine fastq files)
# ##call_all_exome_methods(working_dir, '1060018', 'fastq', {'106001805':['106001805_S41_L005_R1_001.fastq.gz', '106001805_S41_L005_R2_001.fastq.gz']}) ##just the new one
# call_all_exome_methods(working_dir, '1060018', 'fastq', {'106001805':['106001805_S41_L005_R1_001.fastq.gz', '106001805_S41_L005_R2_001.fastq.gz'], '106001801':['106001801_AATCCGTC_R1.fq.gz', '106001801_AATCCGTC_R2.fq.gz'], '106001806a':['106001806a_AATGTTGC_R1.fq.gz', '106001806a_AATGTTGC_R2.fq.gz'], '106001806b':['106001806b_ACACGACC_R1.fq.gz', '106001806b_ACACGACC_R2.fq.gz'], '106001809':['106001809_ACAGATTC_R1.fq.gz', '106001809_ACAGATTC_R2.fq.gz']})



