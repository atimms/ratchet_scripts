#!/usr/bin/env python
import subprocess
import os
import glob
import shutil
import filtering_annotated
import homozygosity_mapping_ub26

##parameters
delim = '\t'
thread_number = '6'

##working dir
# working_dir = '/data/atimms/timon_0317'
working_dir = '/data/atimms/krista_ul_0118'
os.chdir(working_dir)

##programs and files
gatk = '/home/atimms/programs/GenomeAnalysisTK.jar'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
fasta = '/data/atimms/references/mm10/mm10.fa'
bwa = '/home/atimms/programs/bwa-0.7.12/bwa'
samtools = '/home/atimms/programs/samtools-1.3/samtools'
bcftools = '/home/atimms/programs/bcftools-1.3/bcftools'
picard = '/home/atimms/programs/picard-tools-1.139/picard.jar'
delly = '/home/atimms/programs/delly_v0.7.6_CentOS5.4_x86_64bit'
bedtools = '/home/atimms/programs/bedtools2-master/bin/bedtools'
vt = '/data/atimms/gemini/vt/vt'
bgzip = '/home/atimms/programs/htslib-1.3/bgzip'
freebayes = '/home/atimms/programs/freebayes/bin/freebayes'


##files

# st_avinputs = ['K50000022.avinput', 'K50000024.avinput', 'K50000095.avinput', 'combined.avinput']


##annovar parameters
av_genome = 'mm10'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/data/atimms/references/annovar/' + av_genome]
##filter on daredevil, bub, pleather and some of jabiers mice
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,snp142,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic', '-genericdbfile', 'pleather_sc_st.K416.avinput,daredevil.avinput,bub.avinput,J308.st.avinput,J318.st.avinput,J320.st.avinput,J327.st.avinput,J328.st.avinput,J329.st.avinput,mgp.v5.merged.snps_all.dbSNP142.avinput']
av_operation = ['-operation', 'g,r,r,f,f,f,f,f,f,f,f,f,f,f']
av_options = ['-otherinfo', '-remove','-arg', '-splicing 10 ,,,,,,,,,,,,,']
#av_options = ['-otherinfo', '-remove', '-vcfinput']


##filtering and hom mapping parameters
##filtering
col_exon = 6
exon_definition = ['exonic', 'splicing']
col_function = 9
syn_definition = 'synonymous SNV'
zygosity_col = 24
cov_col = 26
cov_definition = 5
qual_col = 25
qual_definition = 30
##het snp mapping 
genome_fai = '/data/atimms/references/mm10/mm10.fa.fai'
window_size = [1000000,5000000,2000000]
# window_size = [10000000]
step_size = 1000000
info_col = 36
naf_values = [0.8,0.9,0.95]



##methods
##call samtools on bamfiles
def variant_calling_samtools(samples, bam_suffix, final_vcf_suffix):
	for sample in samples:
		bamfile = sample + bam_suffix
		vcf_temp1 = sample + '.temp_st.vcf.gz'
		vcf_temp2 = sample + '.temp_st2.vcf.gz'
		final_vcf = sample + final_vcf_suffix
		stmp = subprocess.Popen([samtools,'mpileup', '-ug','-t', 'DP,DV,DPR', '-q', '20', '-C', '50', '-f', fasta, bamfile], stdout=subprocess.PIPE)
		bcft = subprocess.Popen([bcftools,'call', '-vmO', 'z', '-V', 'indels', '-o', vcf_temp1], stdin=stmp.stdout)
		bcft.wait()
		bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1])
		bcf_index.wait()
		#split multi-allelic variants calls in separate lines
		bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-o', vcf_temp2, vcf_temp1])
		bcf_norm1.wait()
		bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf +'.gz', vcf_temp2])
		bcf_norm2.wait()
		bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
		bcf_index.wait()

def variant_calling_gatk(samples, bam_suffix, final_vcf_suffix):
	for sample in samples:
		bamfile = sample + bam_suffix
		vcf_temp0 = sample + '.temp_gatk.vcf'
		vcf_raw_snps = sample + '.temp_gatk_snps.vcf'
		vcf_filtered_snps = sample + '.temp_gatk_snps_filtered.vcf'
		vcf_temp1 = sample + '.temp_gatk_snps2.vcf.gz'
		final_vcf = sample + final_vcf_suffix
		##run haplotype caller
		gatk_hc = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fasta, '-I', bamfile, '-stand_call_conf', '30', '-stand_emit_conf', '10', '-nct', '20', '-o', vcf_temp0])
		gatk_hc.wait()
		##split data in SNPs and indels and apply manual variant filtering (add -L)
		snp_cut = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'SelectVariants', '-R',fasta, '-nt', '5', '-V', vcf_temp0, '-o', vcf_raw_snps, '-selectType', 'SNP'])
		snp_cut.wait()
		snp_vf = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'VariantFiltration', '-R',fasta, '-V', vcf_raw_snps, '-o', vcf_filtered_snps, '--filterExpression ', "QD < 2.0 || FS > 60.0 || SOR > 4.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0", '--filterName', "snp_filter"])
		snp_vf.wait()
		bgzip_cmd = subprocess.Popen([bgzip, vcf_filtered_snps])
		bgzip_cmd.wait()
		bcf_index = subprocess.Popen([bcftools, 'index', vcf_filtered_snps + '.gz'])
		bcf_index.wait()
		#split multi-allelic variants calls in separate lines, and left normalize indels
		bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-O', 'z', '-o', vcf_temp1, vcf_filtered_snps + '.gz'])
		bcf_norm1.wait()
		bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf + '.gz', vcf_temp1])
		bcf_norm2.wait()
		bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
		bcf_index.wait()

##call freebayes in all sample listed in list of bamfiles
def variant_calling_freebayes(samples, bam_suffix, final_vcf_suffix):
	for sample in samples:
		bamfile = sample + bam_suffix
		vcf_temp1 = sample + '.temp_fb.vcf'
		final_vcf = sample + final_vcf_suffix
		with open(vcf_temp1, 'w') as vcf_fh:
			##iXu means no indles, multiallelic or complex
			freebayes_run = subprocess.Popen([freebayes,'-f' , fasta, '-iXu', bamfile], stdout=vcf_fh)
			freebayes_run.wait()
		bgzip_run = subprocess.Popen([bgzip, vcf_temp1])
		bgzip_run.wait()
		bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1 + '.gz'])
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


def intesect_two_vcf_files(vcf1, vcf2, output_dir):
	bcf_isec = subprocess.Popen([bcftools, 'isec', vcf1, vcf2, '-p', output_dir ])
	bcf_isec.wait()

def intersect_variants(samples, freebayes_suffix, gatk_suffix, intersected_suffix):
	for sample in samples:
		vcf1 = sample + freebayes_suffix
		vcf2 = sample + gatk_suffix
		out_dir = sample + intersected_suffix
		bcf_index = subprocess.Popen([bcftools, 'index', vcf2])
		bcf_index.wait()
		intesect_two_vcf_files(vcf1, vcf2, out_dir)


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

##make avinput files
def convert_to_annovar(samples, vcf_prefix):
	for sample in samples:
		vcf = sample + vcf_prefix
		con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', 'temp'])
		con_ann.wait()
	temp_files = glob.glob('temp*.avinput')
	for temp_file in temp_files:
		real_file = temp_file[5:]
		os.rename(temp_file, real_file)
		# shutil.copy(real_file, str(av_ref_dir[0]))

def run_table_annovar(samples):
	for sample in samples:
		avinput = sample + '.avinput'
		av_prefix = avinput.rsplit('.',1)[0]
		command = [table_annovar] + av_buildver + [avinput] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', av_prefix]
		annovar = subprocess.Popen(command)
		annovar.wait()

def multianno_to_annotated(avinputs): 
	head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 'genomicSuperDups', 'snp142','pleather','daredevil','bub','J308','J318','J320','J327','J328','J329', 'mgp.v5.snps', 'Zygosity', 'Qual', 'Coverage', 'Chr', 'Pos', 'Filter', 'Ref2', 'Alt2','Qual', 'Filter', 'GT_info', 'Format', 'Info']
	head_out = delim.join(head + ['\n'])
	for sample in samples:
		avinput = sample + '.avinput'
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

def get_shared_snps(samples, file_suffix, out_file):
	snp_dict = {}
	##get all vars and put in a dict to get counts
	for sample in samples:
		in_file = sample + file_suffix
		with open(in_file, "r") as in_fh:
			line_count = 0
			for line in in_fh:
				line = line.split(delim)
				line_count += 1
				if line_count > 1:
					snp = '_'.join(line[:5])
					if snp in snp_dict:
						snp_dict[snp].append(sample)
					else:
						snp_dict[snp] = [sample]
	##get snps in all three sample
	snp_list = []
	for s in snp_dict:
		if len(snp_dict[s]) == 3:
			snp_list.append(s)
			# print s, snp_dict[s]
	print str(len(snp_list)) + ' snps shared in all sample'
	##make final file
	infile = samples[0] + file_suffix
	with open(infile, "r") as infh, open(out_file, "w") as outfh:
		lc = 0
		for line in infh:
			line = line.split(delim)
			lc += 1
			if lc == 1:
				outfh.write(delim.join(line))
			else:
				snp = '_'.join(line[:5])
				if snp in snp_list:
					outfh.write(delim.join(line))

def make_bed_from_ann_txt(infile, outfile):
	with open(outfile, "w") as out_fh, open(infile, "r") as in_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			if line_count >1:
				line = line.rstrip().split(delim)
				out_fh.write(delim.join(line[:3] + ['\n']))





##call methods
##parameters
project_name = 'krista_ol_comparison'
samples = ['K50000022', 'K50000024', 'K50000095', 'combined_ol', 'fosse_cp', 'fosse_ex_comb', 'timon_comb']
# samples = ['fosse_cp', 'fosse_ex_comb', 'timon_comb']
mkdup_bam = '.bwa_mkdup.bam'
st_vcf_suffix = '.st.vcf'
freebayes_vcf_suffix = '.freebayes.vcf'
gatk_vcf_suffix = '.gatk.vcf'
samtools_vcf_suffix = '.samtools.vcf'
intersected_vcf_suffix = '.intersected_vcfs'
exome_bed_file = '/data/atimms/references/mm10/mm10_refGene_coding_exons.bed'

##coverage
# calculate_exome_coverage(samples, mkdup_bam, exome_bed_file, project_name)

##snp calling etc fb and gatk
# variant_calling_freebayes(samples, mkdup_bam, freebayes_vcf_suffix)
# variant_calling_gatk(samples, mkdup_bam, gatk_vcf_suffix)
##interect variants
# intersect_variants(samples, freebayes_vcf_suffix + '.gz', gatk_vcf_suffix + '.gz', intersected_vcf_suffix)
##annotate -- need to update and do all samples!!!!
# convert_to_annovar(samples, intersected_vcf_suffix + '/0003.vcf')##needed to manually change name of combined avinput file
# run_table_annovar(samples)
# multianno_to_annotated(samples)

##samtools seems better so use old data to get variants
# run_table_annovar(samples)
# multianno_to_annotated(samples)



##filter variants for variants,  homozygsity mapping then counts
'''
##filter variants for candidates snps
for sample in samples:
	##exonic_variants
	filtering_annotated.filter(working_dir, "or", sample + '.annotated.txt' , sample + "_1.temp", [col_exon, col_exon], ['==','=='], [exon_definition[0],exon_definition[1]])
	##remove synonymous
	filtering_annotated.filter(working_dir, "and", sample + "_1.temp", sample + "_2.temp", [col_function], ['!='], [syn_definition])
	##remove if in dbsnp, sanger, or other mouse line
	filtering_annotated.filter(working_dir, "and", sample + "_2.temp", sample + "_3.temp", [13,14,15,16,17,18,19,20,21,22,23], ['==','==','==','==','==','==','==','==','==','==','=='], ['','','','','','','','','','',''])
	##keep if hom
	filtering_annotated.filter(working_dir, "and", sample + "_3.temp", sample + '.hom_exonic_rare.xls', [zygosity_col], ['=='], ['hom'])
	##filter variants by coverage and quality 
	filtering_annotated.filter(working_dir, "and", sample + '.hom_exonic_rare.xls', sample + '.hom_exonic_rare_qual_filtered.xls', [cov_col,qual_col], ['>=','>='], [cov_definition,qual_definition])
'''

	
'''
##homozygosity mapping
# window_size = [100000,500000,1000000,2000000]
window_size = [10000000]
step_size = 100000
for ws in window_size:
	for sample in samples:
		##remove if in dbsnp, sanger, other ped or rmsk
		filtering_annotated.filter(working_dir, "and", sample + '.annotated.txt', sample + "11.temp", [11,12,13,14,15,16,17,18,19,20,21,22,23], ['==','==','==','==','==','==','==','==','==','==','==','==','=='], ['','','','','','','','','','','','',''])
		# filtering_annotated.filter(working_dir, "and", sample + '.annotated.txt', sample + "11.temp", [11,12,23], ['==','==','=='], ['','',''])

		##filter variants by coverage and quality 
		filtering_annotated.filter(working_dir, "and", sample + "11.temp", sample + '.hom_temp.txt', [cov_col,qual_col], ['>=','>='], [cov_definition,qual_definition])

		#make bed file with windows and returns genome name and window size variable
		genome_and_window = homozygosity_mapping_ub26.make_windows(working_dir, genome_fai, ws, step_size).split('.')[0]

		##make bed file from variants
		homozygosity_mapping_ub26.make_bed_from_ann(working_dir, 'gatk_hc', sample + '.hom_temp.txt', zygosity_col, info_col)
		##hom and het count and hom percentage
		homozygosity_mapping_ub26.count_and_percentage(working_dir, genome_and_window, sample + '.bed')
		##naf
		homozygosity_mapping_ub26.naf_in_window(working_dir, genome_and_window, sample + '.bed')
		##total snp number
		homozygosity_mapping_ub26.total_snp_in_window(working_dir, genome_and_window, sample + '.bed')

		##combine bedgraphs for graphing in r
		homozygosity_mapping_ub26.combine_bedgraphs_for_r(working_dir, sample, genome_and_window)
'''
'''
##filter variants for counts
for sample in samples:
	##remove if in dbsnp, sanger, or other mouse lineor rpt region
	filtering_annotated.filter(working_dir, "and", sample + '.annotated.txt', sample + "21.temp", [11,12,13,14,15,16,17,18,19,20,21,22,23], ['==','==','==','==','==','==','==','==','==','==','==','==','=='], ['','','','','','','','','','','','',''])
	##filter variants by coverage and quality 
	filtering_annotated.filter(working_dir, "and", sample + "21.temp", sample + '.all_enu_vars.xls', [cov_col,qual_col], ['>=','>='], [cov_definition,qual_definition])
	##exonic_variants
	filtering_annotated.filter(working_dir, "or", sample + '.all_enu_vars.xls' , sample + '.exonic_enu_vars.xls', [col_exon, col_exon], ['==','=='], [exon_definition[0],exon_definition[1]])
'''



##het modifiers

##parameters
sample_to_get_shared_snps = ['K50000022', 'K50000024', 'K50000095']
enu_file_suffix = '.all_enu_vars.xls'
enu_bed_suffix = '.all_enu_vars.bed'
share_file = 'ol_shared' + enu_file_suffix
window_size = [10000000]
enu_exonic_file_suffix = '.exonic_enu_vars.xls'
share_exonic_file = 'ol_shared' + enu_exonic_file_suffix
samples_to_graph = ['K50000022', 'K50000024', 'K50000095', 'combined_ol', 'fosse_cp', 'fosse_ex_comb', 'timon_comb', 'ol_shared']

##get snps that are in all three samples
# get_shared_snps(sample_to_get_shared_snps, enu_file_suffix, share_file)

##get exonic snps that are in all three samples
get_shared_snps(sample_to_get_shared_snps, enu_exonic_file_suffix, share_exonic_file)

##make bed from enu vars
# for sample in samples_to_graph:
# 	make_bed_from_ann_txt(sample + enu_file_suffix, sample + enu_bed_suffix)

##graph all samples, including in all three
'''
for sample in samples_to_graph:
	for ws in window_size:
		#make bed file with windows and returns genome name and window size variable
		genome_and_window = homozygosity_mapping_ub26.make_windows(working_dir, genome_fai, ws, step_size).split('.')[0]
		print genome_and_window
		window_bed = genome_and_window + '.bed'
		##all shared snps
		out_bed = sample + '.' + genome_and_window + '.bed'
		##bedtools intersect 
		with open('temp.bed', "w") as naf_fh: 
			hom_bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', window_bed, '-b', sample + enu_bed_suffix, '-c'], stdout=naf_fh)
			hom_bt_intersect.wait()
		##add header
		with open(out_bed, "w") as out_fh, open('temp.bed', "r") as in_fh:
			out_fh.write(delim.join(['chr', 'start', 'end', 'snp_number', '\n']))
			for line in in_fh:
				out_fh.write(line)
'''










