#!/tools/BioBuilds-2015.04/bin/python
import subprocess
# import filtering_annotated
# import homozygosity_mapping
import os
import glob
# import calculate_coverage_from_bam

##parameters
delim = '\t'
thread_number = '6'

##working dir
working_dir = '/data/atimms/daredevil_bub'
os.chdir(working_dir)

##programs and files
# gatk = '/Users/atimms/Desktop/ngs/programs/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
# picard_mkdup = '/Users/atimms/Desktop/ngs/programs/picard-tools-1.113/picard-tools-1.113/MarkDuplicates.jar'
fasta = '/data/atimms/daredevil_bub/ref/mm10.fa'
bwa = '/tools/BioBuilds-2015.04/bin/bwa'
samtools = '/home/atimms/programs/samtools-bcftools-htslib-1.0_x64-linux/bin/samtools'
bcftools = '/home/atimms/programs/samtools-bcftools-htslib-1.0_x64-linux/bin/bcftools'
picard = '/home/atimms/programs/picard-tools-1.139/picard.jar'


##dictionaries to combine fq files into individual mice, into phenotyoe, and then dict of final files
dd_fqs_to_combine = {'daredevil.r1.fastq.gz': ['HL5KTCCXX_s1_1_GSLv3-7_06_SL143671.fastq.gz', 'HL5KTCCXX_s1_1_GSLv3-7_07_SL143672.fastq.gz', 'HL5KTCCXX_s1_1_GSLv3-7_08_SL143673.fastq.gz'], 
		'daredevil.r2.fastq.gz': ['HL5KTCCXX_s1_2_GSLv3-7_06_SL143671.fastq.gz', 'HL5KTCCXX_s1_2_GSLv3-7_07_SL143672.fastq.gz', 'HL5KTCCXX_s1_2_GSLv3-7_08_SL143673.fastq.gz']}
dd_fq_dict = {'K416_0006_5': ['HL5KTCCXX_s1_1_GSLv3-7_06_SL143671.fastq.gz', 'HL5KTCCXX_s1_2_GSLv3-7_06_SL143671.fastq.gz'], 'K416_0007_9': ['HL5KTCCXX_s1_1_GSLv3-7_07_SL143672.fastq.gz', 'HL5KTCCXX_s1_2_GSLv3-7_07_SL143672.fastq.gz'], 
		'K417_0024_8': ['HL5KTCCXX_s1_1_GSLv3-7_08_SL143673.fastq.gz', 'HL5KTCCXX_s1_2_GSLv3-7_08_SL143673.fastq.gz'], 'daredevil': ['daredevil.r1.fastq.gz', 'daredevil.r2.fastq.gz']}

bub_fqs_to_combine = {'bub.r1.fastq.gz': ['HL5KTCCXX_s1_1_GSLv3-7_09_SL143674.fastq.gz', 'HL5KTCCXX_s1_1_GSLv3-7_10_SL143675.fastq.gz', 'HL5KTCCXX_s1_1_GSLv3-7_11_SL143676.fastq.gz'], 
		'bub.r2.fastq.gz': ['HL5KTCCXX_s1_2_GSLv3-7_09_SL143674.fastq.gz', 'HL5KTCCXX_s1_2_GSLv3-7_10_SL143675.fastq.gz', 'HL5KTCCXX_s1_2_GSLv3-7_11_SL143676.fastq.gz']}
bub_fq_dict = {'K402_0014_6': ['HL5KTCCXX_s1_1_GSLv3-7_09_SL143674.fastq.gz', 'HL5KTCCXX_s1_2_GSLv3-7_09_SL143674.fastq.gz'], 'K402_0035_1': ['HL5KTCCXX_s1_1_GSLv3-7_10_SL143675.fastq.gz', 'HL5KTCCXX_s1_2_GSLv3-7_10_SL143675.fastq.gz'], 
		'K416_0035_2': ['HL5KTCCXX_s1_1_GSLv3-7_11_SL143676.fastq.gz', 'HL5KTCCXX_s1_2_GSLv3-7_11_SL143676.fastq.gz'], 'bub': ['bub.r1.fastq.gz', 'bub.r2.fastq.gz']}


##annovar parameters
av_genome = 'mm10'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/Users/atimms/Desktop/ngs/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,rmsk,snp137,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic', '-genericdbfile', 'J308.hc.avinput,J308.st.avinput,J318.hc.avinput,J318.st.avinput,J320.hc.avinput,J320.st.avinput,J327.hc.avinput,J327.st.avinput,J328.hc.avinput,J328.st.avinput,J329.hc.avinput,J329.st.avinput,fosse_k421.avinput,fosse_k407.avinput,mgp.v4.snps.dbSNP.avinput']
av_operation = ['-operation', 'g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f']
av_options = ['-otherinfo', '-remove','-arg', '-splicing 10 ,,,,,,,,,,,,,,,,,']
#av_options = ['-otherinfo', '-remove', '-vcfinput']

##created file names
post_bwa_bam = '.bwa.bam'
sorted_bam = '.bwa_sorted.bam'
mkdup_bam = '.bwa_mkdup.bam'
bamslist_file = 'bams.list'
dd_st_vcf = 'daredevil.0316.vcf.gz'
bub_st_vcf = 'bub.0316.vcf.gz'
dd_st_avinputs = []
bub_st_avinputs = []


##filtering and hom mapping parameters
##filtering
col_exon = 6
exon_definition = ['exonic', 'splicing']
col_function = 9
syn_definition = 'synonymous SNV'
zygosity_col = 27
cov_col = 29
cov_definition = 5
qual_col = 28
qual_definition = 30
##hom mapping 
genome_fai = '/Users/atimms/Desktop/ngs/references/mm10/mm10.fa.fai'
window_size = [10000000,5000000,2000000]
step_size = 1000000
info_col = 39
naf_values = [0.8,0.9,0.95]


##combine fq files
def combine_fq_files(file_list, out_fq):
	fq_files = ' '.join(file_list)
	with open(out_fq, 'w') as fq_fh:
		cat_files = subprocess.Popen(['cat'] + file_list, stdout=fq_fh)
		cat_files.wait()

##align with bwa
def align_with_bwa(sample_dict):
	for sample in sample_dict:
		r1_fq = sample_dict[sample][0]
		r2_fq = sample_dict[sample][1]
		print sample, r1_fq, r2_fq
		rg = '@RG\\tID:' + sample + '_RG\\tSM:' + sample + "\\tPL:ILLUMINA"
		pe_bam = sample + post_bwa_bam
		sort_bam = sample + sorted_bam
		pic_dup_bam = sample + mkdup_bam
		bwa_pe = subprocess.Popen([bwa, 'mem', '-M', '-t', '30', '-R', rg, fasta, r1_fq, r2_fq], stdout=subprocess.PIPE)
		st_sam_bam_pe = subprocess.Popen([samtools, 'view', '-q', '20', '-b', '-@', '5', '-o', pe_bam, '-'], stdin=bwa_pe.stdout)
		st_sam_bam_pe.wait()
		# st_sort_pe = subprocess.Popen(['samtools', 'sort','-O', 'bam', '-o', sort_bam, '-@', thread_number, '-m', '2G', '-T', sample, pe_bam])
		st_sort_pe = subprocess.Popen([samtools, 'sort','-O', 'bam', '-o', sort_bam, '-T', sample,'-@', '10', '-m', '10G', pe_bam])
		st_sort_pe.wait()
		# picard_md = subprocess.Popen(['java', '-Xmx24g', '-jar', picard_mkdup, 'VALIDATION_STRINGENCY=SILENT', 'INPUT=' + sort_bam, 'OUTPUT=' + pic_dup_bam, 'CREATE_INDEX=true', 'METRICS_FILE=' + sample + '.metrics'])
		# picard_md.wait()
		##mark duplicates
		picard_md = subprocess.Popen(['java', '-Xmx80g', '-jar', picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=SILENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + pic_dup_bam, 'METRICS_FILE=' + sample + '.metrics', 'TMP_DIR=' + working_dir])
		picard_md.wait()

##make list of all bam files to be analyzed
def make_list_of_bams(sample_dict, bam_suffix, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for sample in sample_dict:
			bam = sample + bam_suffix
			outfh.write(bam + '\n')


##call samtools on bamfiles
def variant_calling_samtools(bamlist, final_vcf):
	vcf_temp1 = 'temp_st.vcf.gz'
	stmp = subprocess.Popen([samtools,'mpileup', '-ug','-t', 'DP,DV,DPR', '-q', '20', '-C', '50', '-f', fasta, '-b', bamlist], stdout=subprocess.PIPE)
	bcft = subprocess.Popen([bcftools,'call', '-vmO', 'z', '-V', 'indels', '-o', vcf_temp1], stdin=stmp.stdout)
	bcft.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines
	bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-O', 'z', '-o', final_vcf, vcf_temp1])
	bcf_norm1.wait()

##make avinput files
def convert_to_annovar(vcf):
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', 'temp'])
	con_ann.wait()
	temp_files = glob.glob('temp*.avinput')
	for temp_file in temp_files:
		os.rename(temp_file, temp_file[5:])


def run_table_annovar(avinputs):
	for avinput in avinputs:
		av_prefix = avinput.rsplit('.',1)[0]
		command = [table_annovar] + av_buildver + [avinput] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', av_prefix]
		annovar = subprocess.Popen(command)
		annovar.wait()

def multianno_to_annotated(avinputs): 
	head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 'snp137', 'J308.hc', 'J308.st', 'J318.hc', 'J318.st', 'J320.hc', 'J320.st', 'J327.hc', 'J327.st', 'J328.hc', 'J328.st', 'J329.hc', 'J329.st', 'fosse_k421', 'fosse_k407', 'mgp.v4.snps.dbSNP', 'Zygosity', 'Qual', 'Coverage', 'Chr', 'Pos', 'Filter', 'Ref2', 'Alt2','Qual', 'Filter', 'GT_info', 'Format', 'Info']
	head_out = delim.join(head + ['\n'])
	for avinput in avinputs:
		av_prefix = avinput.rsplit('.',1)[0]
		multianno = av_prefix + '.mm10_multianno.txt'
		annotated = av_prefix + '.annotated.txt'
		with open(multianno, "r") as multi, open(annotated, "w") as final:
			final.write(head_out)
			line_count = 0
			for line in multi:
				line_count += 1
				if line_count > 1:
					final.write(line)





##call methods, samtools, annovar and format multianno

##daredevil
##combine fqs for each phenotype
# for fq in dd_fqs_to_combine:
# 	print 'combining files: %s to make file %s'% (dd_fqs_to_combine[fq], fq)
# 	combine_fq_files(dd_fqs_to_combine[fq], fq)

##map with bwa and process with samtools etc
# align_with_bwa(dd_fq_dict)
# make_list_of_bams(dd_fq_dict, mkdup_bam, bamslist_file)
# variant_calling_samtools(bamslist_file, dd_st_vcf)
# convert_to_annovar(dd_st_vcf)
# run_table_annovar(dd_st_avinputs)
# multianno_to_annotated(dd_st_avinputs)

##bub
##combine fqs for each phenotype
# for fq in bub_fqs_to_combine:
# 	print 'combining files: %s to make file %s'% (bub_fqs_to_combine[fq], fq)
# 	combine_fq_files(bub_fqs_to_combine[fq], fq)


##map with bwa and process with samtools etc
# align_with_bwa(bub_fq_dict)
# make_list_of_bams(bub_fq_dict, mkdup_bam, bamslist_file)
variant_calling_samtools(bamslist_file, bub_st_vcf)
convert_to_annovar(bub_st_vcf)
# run_table_annovar(bub_st_avinputs)
# multianno_to_annotated(bub_st_avinputs)



'''
##filter variants for candidates snps
for sample in individual_dict:
	##exonic_variants
	filtering_annotated.filter(working_dir, "or", sample + '.annotated.txt' , sample + "_1.temp", [col_exon, col_exon], ['==','=='], [exon_definition[0],exon_definition[1]])
	##remove synonymous
	filtering_annotated.filter(working_dir, "and", sample + "_1.temp", sample + "_2.temp", [col_function], ['!='], [syn_definition])
	##remove if in dbsnp, sanger, or other mouse line
	filtering_annotated.filter(working_dir, "and", sample + "_2.temp", sample + "_3.temp", [12,13,14,15,16,17,18,19,20,21,22,23,24,25,26], ['==','==','==','==','==','==','==','==','==','==','==','==','==','==','=='], ['','','','','','','','','','','','','','',''])
	##keep if hom
	filtering_annotated.filter(working_dir, "and", sample + "_3.temp", sample + '.hom_exonic_rare.xls', [zygosity_col], ['=='], ['hom'])
	##filter variants by coverage and quality 
	filtering_annotated.filter(working_dir, "and", sample + '.hom_exonic_rare.xls', sample + '.hom_exonic_rare_qual_filtered.xls', [cov_col,qual_col], ['>=','>='], [cov_definition,qual_definition])
'''

	
'''
##homozygosity mapping
for ws in window_size:
	for sample in individual_dict:
		##remove if in dbsnp, sanger, other ped or rmsk
		filtering_annotated.filter(working_dir, "and", sample + '.annotated.txt', sample + "11.temp", [11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26], ['==','==','==','==','==','==','==','==','==','==','==','==','==','==','==','=='], ['','','','','','','','','','','','','','','',''])

		##filter variants by coverage and quality 
		filtering_annotated.filter(working_dir, "and", sample + "11.temp", sample + '.hom_temp.txt', [cov_col,qual_col], ['>=','>='], [cov_definition,qual_definition])

		#make bed file with windows and returns genome name and window size variable
		genome_and_window = homozygosity_mapping.make_windows(working_dir, genome_fai, ws, step_size).split('.')[0]

		##make bed file from variants
		homozygosity_mapping.make_bed_from_ann(working_dir, 'samtools', sample +  '.hom_temp.txt', zygosity_col, info_col)
		##hom and het count and hom percentage
		homozygosity_mapping.count_and_percentage(working_dir, genome_and_window, sample + '.bed')
		##naf
		homozygosity_mapping.naf_in_window(working_dir, genome_and_window, sample + '.bed')
		##total snp number
		homozygosity_mapping.total_snp_in_window(working_dir, genome_and_window, sample + '.bed')

		##combine bedgraphs for graphing in r
		homozygosity_mapping.combine_bedgraphs_for_r(working_dir, sample, genome_and_window)
'''
