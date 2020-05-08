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
working_dir = '/data/atimms/greb1l_0317'
fastq_dir = '/data/atimms/greb1l_0317/fq_files/'
os.chdir(working_dir)

##programs and files
# gatk = '/Users/atimms/Desktop/ngs/programs/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
# picard_mkdup = '/Users/atimms/Desktop/ngs/programs/picard-tools-1.113/picard-tools-1.113/MarkDuplicates.jar'
fasta = '/data/atimms/references/mm10/mm10.fa'
bwa = '/home/atimms/programs/bwa-0.7.12/bwa'
samtools = '/home/atimms/programs/samtools-1.3/samtools'
bcftools = '/home/atimms/programs/bcftools-1.3/bcftools'
picard = '/home/atimms/programs/picard-tools-1.139/picard.jar'
delly = '/home/atimms/programs/delly_v0.7.6_CentOS5.4_x86_64bit'
##files etc
fq_dict = {'J321_756_5': ['J321_756_5.r1.fq.gz', 'J321_756_5.r2.fq.gz']}
post_bwa_bam = '.bwa.bam'
sorted_bam = '.bwa_sorted.bam'
mkdup_bam = '.bwa_mkdup.bam'
bamslist_file = 'bams.list'
st_vcf = 'greb1l_0317.vcf.gz'
st_avinputs = ['J321_756_5.avinput']
delly_exclude_regions = '/data/atimms/references/annovar/mm10/mouse.mm10.excl.tsv'
bedtools = '/home/atimms/programs/bedtools2-master/bin/bedtools'


##annovar parameters
av_genome = 'mm10'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/data/atimms/references/annovar/' + av_genome]
##use krista's mice to filter
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,snp142,generic,generic,generic,generic,generic,generic', '-genericdbfile', 'pleather_sc_st.K416.avinput,daredevil.avinput,bub.avinput,fosse_k421.avinput,fosse_k407.avinput,mgp.v5.merged.snps_all.dbSNP142.avinput']
av_operation = ['-operation', 'g,r,r,f,f,f,f,f,f,f']
av_options = ['-otherinfo', '-remove','-arg', '-splicing 10 ,,,,,,,,,']
#av_options = ['-otherinfo', '-remove', '-vcfinput']


##filtering and hom mapping parameters
##filtering
col_exon = 6
exon_definition = ['exonic', 'splicing']
col_function = 9
syn_definition = 'synonymous SNV'
zygosity_col = 20
cov_col = 22
cov_definition = 3
qual_col = 21
qual_definition = 30
##het snp mapping 
genome_fai = '/data/atimms/references/mm10/mm10.fa.fai'
# window_size = [10000000,5000000,2000000]
window_size = [10000000]
step_size = 1000000
info_col = 32
naf_values = [0.8,0.9,0.95]


def combine_fq_file(fq_dir, samples, final_dir):
	for sample in samples:
		print sample
		r1_fq = final_dir + '/' + sample + '.r1.fq.gz'
		r2_fq = final_dir + '/' + sample + '.r2.fq.gz'
		r1_to_combine = glob.glob(fq_dir + '*' + sample + '*_1.fq.gz')
		r2_to_combine = glob.glob(fq_dir + '*' + sample + '*_2.fq.gz')
		print r1_fq, r1_to_combine
		print r2_fq, r2_to_combine
		with open(r1_fq, 'w') as r1_fh:
			cat_files = subprocess.Popen(['cat'] + r1_to_combine, stdout=r1_fh)
			cat_files.wait()
		with open(r2_fq, 'w') as r2_fh:
			cat_files = subprocess.Popen(['cat'] + r2_to_combine, stdout=r2_fh)
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
		bwa_pe = subprocess.Popen([bwa, 'mem', '-M', '-t', '20', '-R', rg, fasta, r1_fq, r2_fq], stdout=subprocess.PIPE)
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
	vcf_temp2 = 'temp_fb2.vcf.gz'
	stmp = subprocess.Popen([samtools,'mpileup', '-ug','-t', 'DP,DV,DPR', '-q', '20', '-C', '50', '-f', fasta, '-b', bamlist], stdout=subprocess.PIPE)
	bcft = subprocess.Popen([bcftools,'call', '-vmO', 'z', '-V', 'indels', '-o', vcf_temp1], stdin=stmp.stdout)
	bcft.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines
	bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-o', vcf_temp2, vcf_temp1])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
	bcf_index.wait()

##make avinput files
def convert_to_annovar(vcf):
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', 'temp'])
	con_ann.wait()
	temp_files = glob.glob('temp*.avinput')
	for temp_file in temp_files:
		real_file = temp_file[5:]
		os.rename(temp_file, real_file)
		shutil.copy(real_file, str(av_ref_dir[0]))

def run_table_annovar(avinputs):
	for avinput in avinputs:
		av_prefix = avinput.rsplit('.',1)[0]
		command = [table_annovar] + av_buildver + [avinput] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', av_prefix]
		annovar = subprocess.Popen(command)
		annovar.wait()

def multianno_to_annotated(avinputs): 
	head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 'genomicSuperDups', 'snp142', 'pleather','daredevil','bub','fosse_k421','fosse_k407', 'mgp.v5.snps', 'Zygosity', 'Qual', 'Coverage', 'Chr', 'Pos', 'Filter', 'Ref2', 'Alt2','Qual', 'Filter', 'GT_info', 'Format', 'Info']
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
		delly_ali = subprocess.Popen([delly, 'call', '-t', cnv_type, '-x', delly_exclude_regions, '-o', out_vcf, '-g', fasta] + bam_files)
		delly_ali.wait()
		with open(out_vcf, "w") as out_fh:
			bcf_view = subprocess.Popen([bcftools, 'view', out_bcf], stdout=out_fh)
			bcf_view.wait()
 # bcftools view tippy_mouse_1016.delly_DEL.bcf > tippy_mouse_1016.delly_DEL.vcf

##get candidate cnvs
def parse_delly_cnv_vcfs(vcfs, chr_req, pos_req, outfile):
	with open(outfile, "w") as out_fh:
		vcf_count = 0
		for vcf in vcfs:
			vcf_count += 1
			with open(vcf, "r") as in_fh:
				line_count = 0
				for line in in_fh:
					if line[:2] != '##':
						line_count += 1
						if line_count == 1 and vcf_count ==1:
							out_fh.write(line)
						elif line_count > 1:
							line = line.rstrip().split(delim)
							chrom = line[0]
							pos = int(line[1])
							m2_genotype = line[9].split(':')[0]
							m1_genotype = line[10].split(':')[0]
							fvb_genotype = line[11].split(':')[0]
							# print fvb_genotype
							if chrom == chr_req and pos >= pos_req[0] and fvb_genotype == '0/0' and m2_genotype == '1/1' and m1_genotype == '1/1':
								out_fh.write(delim.join(line + ['\n']))

def make_bed_from_ann_txt(infile, outfile):
	with open(outfile, "w") as out_fh, open(infile, "r") as in_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			if line_count >1:
				line = line.rstrip().split(delim)
				out_fh.write(delim.join(line[:3] + ['\n']))


##call methods, samtools, annovar and format multianno


##combine fqs
# combine_fq_file(fastq_dir, fq_dict, working_dir)

##map with bwa and process with samtools etc
# align_with_bwa(fq_dict)
# make_list_of_bams(fq_dict, mkdup_bam, bamslist_file)
# variant_calling_samtools(bamslist_file, st_vcf)
# convert_to_annovar(st_vcf)
# run_table_annovar(st_avinputs)
# multianno_to_annotated(st_avinputs)

##filter vars and hom mapping

##variants seen in all samples
'''
##filter variants for candidates snps
for sample in fq_dict:
	##exonic_variants
	filtering_annotated.filter(working_dir, "or", sample + '.annotated.txt' , sample + "_1.temp", [col_exon, col_exon], ['==','=='], [exon_definition[0],exon_definition[1]])
	##remove synonymous
	filtering_annotated.filter(working_dir, "and", sample + "_1.temp", sample + "_2.temp", [col_function], ['!='], [syn_definition])
	##remove if in dbsnp, sanger, or other mouse line
	filtering_annotated.filter(working_dir, "and", sample + "_2.temp", sample + "_3.temp", [13,14,15,16,17,18,19], ['==','==','==','==','==','==','=='], ['','','','','','',''])
	##keep if hom
	filtering_annotated.filter(working_dir, "and", sample + "_3.temp", sample + '.hom_exonic_rare.xls', [zygosity_col], ['=='], ['hom'])
	##filter variants by coverage and quality 
	filtering_annotated.filter(working_dir, "and", sample + '.hom_exonic_rare.xls', sample + '.hom_exonic_rare_qual_filtered.xls', [cov_col,qual_col], ['>=','>='], [cov_definition,qual_definition])
'''

	
# '''
##homozygosity mapping
for ws in window_size:
	for sample in fq_dict:
		##remove if in dbsnp, sanger, other ped or rmsk
		filtering_annotated.filter(working_dir, "and", sample + '.annotated.txt', sample + "11.temp", [11,12,13,14,15,16,17,18,19], ['==','==','==','==','==','==','==','==','=='], ['','','','','','','','',''])

		##filter variants by coverage and quality 
		filtering_annotated.filter(working_dir, "and", sample + "11.temp", sample + '.hom_temp.txt', [cov_col,qual_col], ['>=','>='], [cov_definition,qual_definition])

		#make bed file with windows and returns genome name and window size variable
		genome_and_window = homozygosity_mapping_ub26.make_windows(working_dir, genome_fai, ws, step_size).split('.')[0]

		##make bed file from variants
		homozygosity_mapping_ub26.make_bed_from_ann(working_dir, 'samtools', sample +  '.hom_temp.txt', zygosity_col, info_col)
		##hom and het count and hom percentage
		homozygosity_mapping_ub26.count_and_percentage(working_dir, genome_and_window, sample + '.bed')
		##naf
		homozygosity_mapping_ub26.naf_in_window(working_dir, genome_and_window, sample + '.bed')
		##total snp number
		homozygosity_mapping_ub26.total_snp_in_window(working_dir, genome_and_window, sample + '.bed')

		##combine bedgraphs for graphing in r
		homozygosity_mapping_ub26.combine_bedgraphs_for_r(working_dir, sample, genome_and_window)
# '''
