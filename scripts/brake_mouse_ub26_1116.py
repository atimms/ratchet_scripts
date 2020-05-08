#!/tools/BioBuilds-2015.04/bin/python
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
working_dir = '/data/atimms/brake_mouse_1116'
fastq_dir = '/data/atimms/brake_mouse_1116/fq_files/'
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
delly = '/home/atimms/programs/delly_v0.7.7_CentOS5.4_x86_64bit'
##files etc
# fq_dict = {'39044': ['39044.r1.fq.gz', '39044.r2.fq.gz'], '58546': ['58546.r1.fq.gz', '58546.r2.fq.gz'], 
# 		'81598': ['81598.r1.fq.gz', '81598.r2.fq.gz'], '16889': ['16889.r1.fq.gz', '16889.r2.fq.gz']}
##new version with bl6
fq_dict = {'39044': ['39044.r1.fq.gz', '39044.r2.fq.gz'], '58546': ['58546.r1.fq.gz', '58546.r2.fq.gz'], 
		'81598': ['81598.r1.fq.gz', '81598.r2.fq.gz'], '16889': ['16889.r1.fq.gz', '16889.r2.fq.gz'], 'C57BL_6NJ': ['C57BL_6NJ.r1.fq', 'C57BL_6NJ.r2.fq']}
bl6_fq_dict = {'C57BL_6NJ': ['C57BL_6NJ.r1.fq', 'C57BL_6NJ.r2.fq']}
post_bwa_bam = '.bwa.bam'


sorted_bam = '.bwa_sorted.bam'
mkdup_bam = '.bwa_mkdup.bam'
bamslist_file = 'bams.list'
st_vcf = 'brake.1116.vcf.gz'
st_avinputs = ['39044.avinput', '58546.avinput', '81598.avinput', '16889.avinput']
delly_exclude_regions = '/data/atimms/references/annovar/mm10/mouse.mm10.excl.tsv'
bedtools = '/home/atimms/programs/bedtools2-master/bin/bedtools'

##annovar parameters
av_genome = 'mm10'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/data/atimms/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,snp142,generic,generic,generic,generic,generic', '-genericdbfile', '39044.avinput,58546.avinput,81598.avinput,16889.avinput,mgp.v5.merged.snps_all.dbSNP142.avinput']
av_operation = ['-operation', 'g,r,r,f,f,f,f,f,f']
av_options = ['-otherinfo', '-remove','-arg', '-splicing 10 ,,,,,,,,']
#av_options = ['-otherinfo', '-remove', '-vcfinput']




##filtering and hom mapping parameters
##filtering
col_exon = 6
exon_definition = ['exonic', 'splicing']
col_function = 9
syn_definition = 'synonymous SNV'
zygosity_col = 19
cov_col = 21
cov_definition = 5
qual_col = 20
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
		r1_to_combine = glob.glob(fq_dir + '*' + sample + '*R1*')
		r2_to_combine = glob.glob(fq_dir + '*' + sample + '*R2*')
		# print r1_fq, r1_to_combine
		# print r2_fq, r2_to_combine
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
		# bwa_pe = subprocess.Popen([bwa, 'mem', '-M', '-t', '20', '-R', rg, fasta, r1_fq, r2_fq], stdout=subprocess.PIPE)
		# st_sam_bam_pe = subprocess.Popen([samtools, 'view', '-q', '20', '-b', '-@', '5', '-o', pe_bam, '-'], stdin=bwa_pe.stdout)
		# st_sam_bam_pe.wait()
		# st_sort_pe = subprocess.Popen([samtools, 'sort','-O', 'bam', '-o', sort_bam, '-T', sample,'-@', '10', '-m', '10G', pe_bam])
		# st_sort_pe.wait()
		# picard_md = subprocess.Popen(['java', '-Xmx24g', '-jar', picard_mkdup, 'VALIDATION_STRINGENCY=SILENT', 'INPUT=' + sort_bam, 'OUTPUT=' + pic_dup_bam, 'CREATE_INDEX=true', 'METRICS_FILE=' + sample + '.metrics'])
		# picard_md.wait()
		##mark duplicates
		picard_md = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=SILENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + pic_dup_bam, 'METRICS_FILE=' + sample + '.metrics', 'TMP_DIR=' + working_dir])
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
	bcft = subprocess.Popen([bcftools,'call', '-vmO', 'z', '-o', vcf_temp1], stdin=stmp.stdout)
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
	head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 'genomicSuperDups', 'snp142', '39044','58546','81598','16889', 'mgp.v5.snps', 'Zygosity', 'Qual', 'Coverage', 'Chr', 'Pos', 'Filter', 'Ref2', 'Alt2','Qual', 'Filter', 'GT_info', 'Format', 'Info']
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

def filter_delly_vcfs(vcfs, outfile):
	with open(outfile, "w") as out_fh:
		vcf_count = 0
		for vcf in vcfs:
			passed_count = 0
			vcf_count += 1
			with open(vcf, "r") as in_fh:
				line_count = 0
				for line in in_fh:
					if line[0] == '#':
						if vcf_count ==1:
							out_fh.write(line)
					else:
						line_count += 1
						line = line.rstrip().split(delim)
						p_filter = line[6]
						brake_genotype = line[9:13]
						brake_genotype = [g.split(':')[0] for g in brake_genotype]
						bl6_genotype = line[13].split(':')[0]
						print brake_genotype, bl6_genotype
						if p_filter == 'PASS' and brake_genotype == ['0/1', '0/1', '0/1', '0/1'] and bl6_genotype != '0/1' and bl6_genotype != '1/1':
						# if p_filter == 'PASS' and brake_genotype == ['0/1', '0/1', '0/1', '0/1']:
							passed_count +=1
							out_fh.write(delim.join(line + ['\n']))
			print 'for vcf %s we have %s varaints of which %s passed the filter'%(vcf,line_count, passed_count)

def make_refgene_dict(refgene_file):
	rg_dict, rg_dict2 = {}, {}
	with open(refgene_file, "r") as in_fh:
		for line in in_fh:
			line = line.split(delim)
			nm = line[1]
			genename = line[12]
			# print nm, genename
			if nm in rg_dict:
				rg_dict[nm].append(genename)
			else:
				rg_dict[nm] = [genename]
	for g in rg_dict:
		print g, rg_dict[g], list(set(rg_dict[g]))
		rg_dict2[g] = list(set(rg_dict[g]))
	return rg_dict2


def get_varaiants_that_intersect_with_exons(in_vcf, outfile):
	av_prefix = in_vcf.rsplit('.',1)[0]
	avinput = av_prefix + '.avinput'
	##make bed file, convert to avinpu then make temp.bed
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', in_vcf, '-includeinfo', '-withfreq', '-allsample', '-outfile', avinput])
	con_ann.wait()
	with open('temp.bed', "w") as out_fh, open(avinput, "r") as in_fh:
		for line in in_fh:
			line = line.split(delim)
			out_fh.write(delim.join(line[:3]) + '\n')
	##complete intersection with bedtools
	with open('temp_int.bed', 'w') as out2_fh:
		bed_int = subprocess.Popen(['bedtools', 'intersect', '-a', 'temp.bed', '-b', 'mm10_refGene_coding_exons.std_chr.bed', '-wa', '-wb'], stdout=out2_fh)
		bed_int.wait()
	##make refgene dict to convery nms to genenames
	refg_dict = make_refgene_dict('mm10_refgene.txt')
	##make dict of intersections
	int_dict, int_dict2 = {},{}
	nms_not_found = 0
	with open('temp_int.bed', 'r') as tint_fh:
		for line in tint_fh:
			line = line.split(delim)
			var = '_'.join(line[:3])
			# print line[6]
			nm = line[6].split('_')[0] + '_' + line[6].split('_')[1]
			if nm in refg_dict:
				gene = refg_dict[nm]
				if var in int_dict:
					if nm not in int_dict[var][0]:
						int_dict[var][0].append(nm)
						int_dict[var][1].extend(gene)
				else:
					int_dict[var] = [[nm],gene]
			else:
				nms_not_found += 1
	print len(int_dict)
	print nms_not_found
	for g in int_dict:
		print g, int_dict[g]
		int_dict2[g] = [list(set(int_dict[g][0])),list(set(int_dict[g][1]))]
	##combine avinput with dict and write final file
	with open(outfile, "w") as out_fh, open(avinput, "r") as in_fh:
		header = ['chr', 'start', 'end', 'id', 'ref', 'alt', 'delly_info', 'format', '16889','81598','39044','58546', 'NMs', 'genes_affected' '\n']
		out_fh.write(delim.join(header))
		for line in in_fh:
			line = line.rstrip().split(delim)
			var = '_'.join(line[:3])
			if var in int_dict2:
				line_out = line[:3] + line[10:13] + line[15:] + [','.join(int_dict2[var][0]), ','.join(int_dict2[var][1]), '\n']
				out_fh.write(delim.join(line_out))

def get_varaiants_that_intersect_with_het_regions(in_vcf, outfile):
	av_prefix = in_vcf.rsplit('.',1)[0]
	avinput = av_prefix + '.avinput'
	##make bed file, convert to avinpu then make temp.bed
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', in_vcf, '-includeinfo', '-withfreq', '-allsample', '-outfile', avinput])
	con_ann.wait()
	with open('temp.bed', "w") as out_fh, open(avinput, "r") as in_fh:
		for line in in_fh:
			line = line.split(delim)
			out_fh.write(delim.join(line[:3]) + '\n')
	##complete intersection with bedtools
	with open('temp_int.bed', 'w') as out2_fh:
		bed_int = subprocess.Popen(['bedtools', 'intersect', '-a', 'temp.bed', '-b', 'enu.het_in_all_samples.qual_30.cov_5_100.ge_5.mm10_10000kb_1000kb.bed', '-wa', '-wb'], stdout=out2_fh)
		bed_int.wait()
	##make dict of intersections
	int_dict = {}
	with open('temp_int.bed', 'r') as tint_fh:
		for line in tint_fh:
			line = line.split(delim)
			var = '_'.join(line[:3])
			# print line[6]
			window = '_'.join(line[3:])
			if var in int_dict:
				int_dict[var].append(window)
			else:
				int_dict[var] = [window]
	print len(int_dict)
	# for g in int_dict:
	# 	print g, int_dict[g]
	##combine avinput with dict and write final file
	with open(outfile, "w") as out_fh, open(avinput, "r") as in_fh:
		header = ['chr', 'start', 'end', 'id', 'ref', 'alt', 'delly_info', 'format', '16889','81598','39044','58546', '\n']
		out_fh.write(delim.join(header))
		for line in in_fh:
			line = line.rstrip().split(delim)
			var = '_'.join(line[:3])
			if var in int_dict:
				line_out = line[:3] + line[10:13] + line[15:] + ['\n']
				out_fh.write(delim.join(line_out))

def delly_cnv_analysis(vcf_prefix, samples, bam_suffix):
	bam_files = []
	filtered_vcf = vcf_prefix + 'filtered.vcf'
	exonic_filtered_vars = vcf_prefix + 'filtered.exonic.xls'
	het_regions_filtered_vars = vcf_prefix + 'filtered.het_regions.xls'
	##make list of bam files
	for sample in samples:
		bam = sample + bam_suffix
		bam_files.append(bam)
	##run delly
	delly_vcfs = []
	for cnv_type in ['DEL', 'DUP', 'INV', 'BND', 'INS']:
		out_bcf = vcf_prefix + 'delly_' + cnv_type + '.bcf'
		out_vcf = vcf_prefix + 'delly_' + cnv_type + '.vcf'
		##temp name
		delly_vcfs.append(out_vcf)
		'''
		delly_ali = subprocess.Popen([delly, 'call', '-t', cnv_type, '-x', delly_exclude_regions, '-o', out_bcf, '-g', fasta] + bam_files)
		delly_ali.wait()
		bcf_view = subprocess.Popen([bcftools, 'view', '-o',out_vcf , out_bcf])
		bcf_view.wait()
		'''
	filter_delly_vcfs(delly_vcfs, filtered_vcf)
	get_varaiants_that_intersect_with_exons(filtered_vcf, exonic_filtered_vars)
	get_varaiants_that_intersect_with_het_regions(filtered_vcf, het_regions_filtered_vars)

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
##just running bl6
# align_with_bwa(bl6_fq_dict)
# make_list_of_bams(fq_dict, mkdup_bam, bamslist_file)
# variant_calling_samtools(bamslist_file, st_vcf)
# convert_to_annovar(st_vcf)
# run_table_annovar(st_avinputs)
# multianno_to_annotated(st_avinputs)


##run delly
samples_for_delly = ['39044', '58546', '81598', '16889', 'C57BL_6NJ']
delly_cnv_analysis('brake_mouse.delly_cnvs.', samples_for_delly, '.bwa_mkdup.bam')

##filter vars and hom mapping
# samples = ['39044', '58546', '81598', '16889']

##variants seen in all samples
'''
samples = ['16889']
for sample in samples:
	##remove indels
	filtering_annotated.filter(working_dir, "and", sample + ".annotated.txt", sample + "_1.temp", [4,5], ['!=','!='], ['-','-'])
	##remove if in rmsk, segdup, dbsnp, sanger, or other mouse line
	filtering_annotated.filter(working_dir, "and", sample + "_1.temp", sample + "_2.temp", [11,12,13,18], ['==','==','==','=='], ['','','',''])
	# filtering_annotated.filter(working_dir, "and", sample + "_1.temp", sample + "_2.temp", [13,18,19], ['==','==','=='], ['','',''])
	##keep if het in all samples
	filtering_annotated.filter(working_dir, "and", sample + "_2.temp", "enu.het_in_all_samples.xls", [14,15,16,17], ['==','==','==','=='], ['het','het','het','het'])
	##filter variants by coverage and quality  -- not used
	# filtering_annotated.filter(working_dir, "and", sample + "_2.temp", sample + ".enu.het_in_all_samples.xls", [cov_col,qual_col], ['>=','>='], [cov_definition,qual_definition])
	##exonic_variants
	filtering_annotated.filter(working_dir, "or", "enu.het_in_all_samples.xls" , sample + "_3.temp", [col_exon, col_exon], ['==','=='], [exon_definition[0],exon_definition[1]])
	# filtering_annotated.filter(working_dir, "or", sample + ".enu.het_in_all_samples.xls" , sample + "_3.temp", [col_exon, col_exon], ['==','=='], [exon_definition[0],exon_definition[1]])
	##remove synonymous
	filtering_annotated.filter(working_dir, "and", sample + "_3.temp", "enu.het_in_all_samples.exonic.xls", [col_function], ['!='], [syn_definition])
	##filter variants by coverage and quality -- not used
	# filtering_annotated.filter(working_dir, "and", sample + "_2.temp", sample + ".enu.het_in_all_samples.xls", [cov_col,qual_col], ['>=','>='], [cov_definition,qual_definition])
'''

'''
##take all het snps in all mice and filter by coverage and quality
##filter variants by coverage >5 and <100 and quality >30
filtering_annotated.filter(working_dir, "and", "enu.het_in_all_samples.xls", "enu.het_in_all_samples.qual_30.cov_5_100.xls", [cov_col,cov_col,qual_col], ['>=','<=','>='], [5,100,30])
filtering_annotated.filter(working_dir, "and", "enu.het_in_all_samples.xls", "enu.het_in_all_samples.qual_30.cov_10_100.xls", [cov_col,cov_col,qual_col], ['>=','<=','>='], [10,100,30])
'''

'''
##don't filter using rmsk/segdups
samples = ['16889']
for sample in samples:
	##remove indels
	filtering_annotated.filter(working_dir, "and", sample + ".annotated.txt", sample + "_1.temp", [4,5], ['!=','!='], ['-','-'])
	##remove if in rmsk, segdup, dbsnp, sanger, or other mouse line
	filtering_annotated.filter(working_dir, "and", sample + "_1.temp", sample + "_2.temp", [13,18,19], ['==','==','=='], ['','',''])
	##keep if het in all samples
	filtering_annotated.filter(working_dir, "and", sample + "_2.temp", "extra.temp", [14,15,16,17], ['==','==','==','=='], ['het','het','het','het'])
	##filter variants by coverage and quality 
	# filtering_annotated.filter(working_dir, "and", sample + "_2.temp", sample + ".enu.het_in_all_samples.xls", [cov_col,qual_col], ['>=','>='], [cov_definition,qual_definition])
	##exonic_variants
	filtering_annotated.filter(working_dir, "or", "extra.temp" , sample + "_3.temp", [col_exon, col_exon], ['==','=='], [exon_definition[0],exon_definition[1]])
	# filtering_annotated.filter(working_dir, "or", sample + ".enu.het_in_all_samples.xls" , sample + "_3.temp", [col_exon, col_exon], ['==','=='], [exon_definition[0],exon_definition[1]])
	##remove synonymous
	filtering_annotated.filter(working_dir, "and", sample + "_3.temp", "enu.het_in_all_samples.exonic.not_rpt_filtered.xls", [col_function], ['!='], [syn_definition])

'''


##graph het snps
##make bedfile from all het snps
# make_bed_from_ann_txt('enu.het_in_all_samples.xls', 'enu.het_in_all_samples.bed')
# make_bed_from_ann_txt('enu.het_in_all_samples.qual_30.cov_10_100.xls', 'enu.het_in_all_samples.qual_30.cov_10_100.bed')
# make_bed_from_ann_txt('enu.het_in_all_samples.qual_30.cov_5_100.xls', 'enu.het_in_all_samples.qual_30.cov_5_100.bed')


#bedtools intersect comparing snps to genome
'''
for ws in window_size:
	#make bed file with windows and returns genome name and window size variable
	genome_and_window = homozygosity_mapping_ub26.make_windows(working_dir, genome_fai, ws, step_size).split('.')[0]
	print genome_and_window
	window_bed = genome_and_window + '.bed'
	##all shared snps
	out_bed = 'enu.het_in_all_samples.' + genome_and_window + '.bed'
	##bedtools intersect 
	with open('temp.bed', "w") as naf_fh: 
		hom_bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', window_bed, '-b', 'enu.het_in_all_samples.bed', '-c'], stdout=naf_fh)
		hom_bt_intersect.wait()
	##add header
	with open(out_bed, "w") as out_fh, open('temp.bed', "r") as in_fh:
		out_fh.write(delim.join(['chr', 'start', 'end', 'snp_number', '\n']))
		for line in in_fh:
			out_fh.write(line)
	##cov 5/100
	out_bed = 'enu.het_in_all_samples.qual_30.cov_10_100.' + genome_and_window + '.bed'
	##bedtools intersect 
	with open('temp.bed', "w") as naf_fh: 
		hom_bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', window_bed, '-b', 'enu.het_in_all_samples.qual_30.cov_10_100.bed', '-c'], stdout=naf_fh)
		hom_bt_intersect.wait()
	##add header
	with open(out_bed, "w") as out_fh, open('temp.bed', "r") as in_fh:
		out_fh.write(delim.join(['chr', 'start', 'end', 'snp_number', '\n']))
		for line in in_fh:
			out_fh.write(line)
	##cov 10/100
	out_bed = 'enu.het_in_all_samples.qual_30.cov_5_100.' + genome_and_window + '.bed'
	##bedtools intersect 
	with open('temp.bed', "w") as naf_fh: 
		hom_bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', window_bed, '-b', 'enu.het_in_all_samples.qual_30.cov_5_100.bed', '-c'], stdout=naf_fh)
		hom_bt_intersect.wait()
	##add header
	with open(out_bed, "w") as out_fh, open('temp.bed', "r") as in_fh:
		out_fh.write(delim.join(['chr', 'start', 'end', 'snp_number', '\n']))
		for line in in_fh:
			out_fh.write(line)
'''
