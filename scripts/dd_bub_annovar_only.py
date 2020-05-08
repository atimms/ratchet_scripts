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
av_ref_dir = ['/data/atimms/daredevil_bub/ref']
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
dd_st_avinputs = ['K416_0006_5.avinput', 'K416_0007_9.avinput', 'K417_0024_8.avinput', 'daredevil.avinput']
bub_st_avinputs = ['K402_0014_6.avinput', 'K402_0035_1.avinput', 'K416_0035_2.avinput', 'bub.avinput']


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

#convert_to_annovar(dd_st_vcf)
#run_table_annovar(dd_st_avinputs)
#multianno_to_annotated(dd_st_avinputs)

##bub
convert_to_annovar(bub_st_vcf)
run_table_annovar(bub_st_avinputs)
multianno_to_annotated(bub_st_avinputs)



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
