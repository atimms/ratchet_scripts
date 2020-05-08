#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import filtering_annotated

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
amgen_sample_info = 'amgen_sample_fq_info.txt'
# amgen_sample_info = 'amgen_sample_fq_info_3.txt'

##samples
bri_samples = ['N353', 'T112', 'T147', 'T148', 'T149', 'T150', 'T151', 'T152', 'T154', 'T155', 'T241', 'T352', 'T354', 'T355', 
		'T356', 'T359', 'T361', 'T362', 'T363', 'T365', 'T366', 'T368', 'T370', 'T372', 'T373']

#annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/data/atimms/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,dbnsfp30a,popfreq_all_20150413,avsnp147']
av_operation = ['-operation', 'g,r,r,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput']

##definitions for determining var function
exonic_definitions = ['exonic', 'splicing']
non_protein_changing_definitions =['synonymous SNV', 'synonymous_SNV']
no_record_definition = ['', '.']
##list of columns for both methods (all start at 1):
##for filtering variants: refgene function, frequency dbs, refgene exonic function, refgene geneame, quality and coverage
filtering_cols = [15, 56, 18, 16, [7,8]]

frequencies = [0.01]



def combine_fq_file(samples):
	for sample in samples:
		print sample
		r1_fq = sample + '_1.fastq.gz'
		r2_fq = sample + '_2.fastq.gz'
		r1_to_combine = glob.glob(sample + '*R1*fastq.gz')
		r2_to_combine = glob.glob(sample + '*R2*fastq.gz')
		print r1_fq, r1_to_combine
		print r2_fq, r2_to_combine
		with open(r1_fq, 'w') as r1_fh:
			cat_files = subprocess.Popen(['cat'] + r1_to_combine, stdout=r1_fh)
			cat_files.wait()
		with open(r2_fq, 'w') as r2_fh:
			cat_files = subprocess.Popen(['cat'] + r2_to_combine, stdout=r2_fh)
			cat_files.wait()

def combine_amgen_fq_files(sample_dict):
	for sample in sample_dict:
		print sample
		r1_fq = sample + '_1.fq.gz'
		r2_fq = sample + '_2.fq.gz'
		r1_to_combine = sample_dict[sample][0]
		r2_to_combine = sample_dict[sample][1]
		print r1_fq, r1_to_combine
		print r2_fq, r2_to_combine
		with open(r1_fq, 'w') as r1_fh:
			cat_files = subprocess.Popen(['cat'] + r1_to_combine, stdout=r1_fh)
			cat_files.wait()
			for f in r1_to_combine:
				os.remove(f)
		with open(r2_fq, 'w') as r2_fh:
			cat_files = subprocess.Popen(['cat'] + r2_to_combine, stdout=r2_fh)
			cat_files.wait()
			for f in r2_to_combine:
				os.remove(f)

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


def variant_calling_gatk_hc(samples, bam_suffix, vcf_suffix):
	#java -Xmx32g -jar gatk -T UnifiedGenotyper -R hg19_fasta -rf BadCigar -allowPotentiallyMisencodedQuals -L MIPtargets.intervals -I bam.realigned.list -o vcf -dcov 5000 -dt NONE -glm both -A AlleleBalance -stand_emit_conf 10 -minIndelFrac 0.005 -G Standard
	#java -Xmx8g -jar gatk -T VariantFiltration -R hg19_fasta -rf BadCigar -allowPotentiallyMisencodedQuals -V in_vcf -o out.vcf -window 20 -cluster 5 -filterName ABFilter -filter "ABHet > 0.75" -filterName QDFilter -filter "QD < 5.0" -filterName QUALFilter -filter "QUAL < 30.0" -filterName LowCoverage -filter "DP < 5"
	vcf_files = []
	for sample in samples:
		bam = sample + bam_suffix
		vcf_temp1 = 'temp_ug1' + sample + vcf_suffix
		vcf_temp2 = 'temp_ug2' + sample + vcf_suffix
		final_vcf = sample + vcf_suffix
		vcf_files.append(final_vcf)
		gatk_hc = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fa_file, '-I', bam, '-o', vcf_temp1, '-dontUseSoftClippedBases', '-stand_call_conf', '20.0', '-stand_emit_conf', '20.0'])
		gatk_hc.wait()
		##manual filtering
		gatk_filter = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'VariantFiltration', '-R', fa_file, '-V', vcf_temp1, '-o', vcf_temp2, '-window', '35', '-cluster', '3', '-filterName', 'FS', '-filter' ,'FS>30.0', '-filterName', 'QD', '-filter', 'QD<2.0'])
		gatk_filter.wait()
		bgzip_run = subprocess.Popen(['bgzip', vcf_temp2])
		bgzip_run.wait()
		bcf_index = subprocess.Popen([bcftools_12, 'index', vcf_temp2 + '.gz'])
		bcf_index.wait()
		with open (final_vcf, 'w') as tvcf_fh:
			zless_vcf = subprocess.Popen(['zless', vcf_temp2 + '.gz'], stdout=subprocess.PIPE)
			sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
			vt_decompose = subprocess.Popen([vt, 'decompose', '-s', '-'], stdin=sed_vcf.stdout, stdout=subprocess.PIPE)
			vt_normalize = subprocess.Popen([vt, 'normalize', '-r', fa_file, '-'], stdin=vt_decompose.stdout, stdout=tvcf_fh)
			vt_normalize.wait()
	return vcf_files


def variant_calling_freebayes(samples, bam_suffix, vcf_suffix):    #if >~100 bams used increase ulimit i.e. ulimit -n 10000
	#freebayes -f ref.fa aln.bam >var.vcf
	vcf_files = []
	for sample in samples:
		bam = sample + bam_suffix
		vcf_temp1 = sample +'temp_fb1' + vcf_suffix
		vcf_temp2 = sample +'temp_fb2' + vcf_suffix
		final_vcf = sample + vcf_suffix
		vcf_files.append(final_vcf)
		with open(vcf_temp1, 'w') as vcf_fh:
			run_freebayes = subprocess.Popen([freebayes, '-f', fa_file, bam], stdout=vcf_fh)
			run_freebayes.wait()
		bgzip_run = subprocess.Popen(['bgzip', vcf_temp1])
		bgzip_run.wait()
		bcf_index = subprocess.Popen([bcftools_12, 'index', vcf_temp1 + '.gz'])
		bcf_index.wait()
		with open (final_vcf, 'w') as tvcf_fh:
			zless_vcf = subprocess.Popen(['zless', vcf_temp1 + '.gz'], stdout=subprocess.PIPE)
			sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
			vt_decompose = subprocess.Popen([vt, 'decompose', '-s', '-'], stdin=sed_vcf.stdout, stdout=subprocess.PIPE)
			vt_normalize = subprocess.Popen([vt, 'normalize', '-r', fa_file, '-'], stdin=vt_decompose.stdout, stdout=tvcf_fh)
			vt_normalize.wait()
	return vcf_files

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
def table_annovar_vcf(vcf, project_prefix):
	out_prefix = project_prefix
	command = [table_annovar] + av_buildver + [vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', out_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()

def split_info_field(info_list):
	indices = [i for i, s in enumerate(info_list) if 'ANNOVAR_DATE' in s]
	# print indices
	i_count = 0
	final_list = []
	for info in info_list:
		# print info
		if i_count > indices[0] and info != 'ALLELE_END':
			info2 = info.split('=')[1]
			#print info2
			final_list.append(info2)
		i_count += 1
	return final_list
	
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

def format_avinput(description, sample_list):
	for sample in sample_list:
		head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Zygosity', 'Qual', 'Coverage', 'Filter', 'Pos', 'Ref2', 'Alt2', 'Format', 'Info', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 'genomicSuperDups', 'SIFT_score', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score', 'Polyphen2_HVAR_pred', 'LRT_score', 'LRT_pred', 'MutationTaster_score', 'MutationTaster_pred', 'MutationAssessor_score', 'MutationAssessor_pred', 'FATHMM_score', 'FATHMM_pred', 'PROVEAN_score', 'PROVEAN_pred', 'VEST3_score', 'CADD_raw', 'CADD_phred', 'DANN_score', 'fathmm-MKL_coding_score', 'fathmm-MKL_coding_pred', 'MetaSVM_score', 'MetaSVM_pred', 'MetaLR_score', 'MetaLR_pred', 'integrated_fitCons_score', 'integrated_confidence_value', 'GERP++_RS', 'phyloP7way_vertebrate', 'phyloP20way_mammalian', 'phastCons7way_vertebrate', 'phastCons20way_mammalian', 'SiPhy_29way_logOdds', 'PopFreqMax', '1000G_ALL', '1000G_AFR', '1000G_AMR', '1000G_EAS', '1000G_EUR', '1000G_SAS', 'ExAC_ALL', 'ExAC_AFR', 'ExAC_AMR', 'ExAC_EAS', 'ExAC_FIN', 'ExAC_NFE', 'ExAC_OTH', 'ExAC_SAS', 'ESP6500siv2_ALL', 'ESP6500siv2_AA', 'ESP6500siv2_EA', 'CG46', 'avsnp147']
		head_out = delim.join(head + ['\n'])
		avinput = description + '.' + sample + '.avinput'
		outfile = description + '.' + sample + '.annotated.txt'
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


##calls all annovar methods
def annotate_vcf_file(vcf, project_prefix, samples):
	table_annovar_vcf(vcf, project_prefix)
	post_annovar_vcf = project_prefix + '.' + av_genome + '_multianno.vcf'
	convert_to_annovar(post_annovar_vcf, project_prefix)
	format_avinput(project_prefix, samples)

##simple list maker -- give string and dictionary to get length from
def make_list(string, dict):
	l = []
	for i in range(len(dict)):
		l.append(string)
	return l

def filter_ann_txt_files(samples, cols_to_filter, freq_req, outfile_prefix):
	##filter rnaseq data
	for sample in samples:
		for freq in freq_req:
			##only 'rare' (<=1%)
			filtering_annotated.filter(working_dir, "and", outfile_prefix + '.' + sample + '.annotated.txt', "1.temp", [cols_to_filter[1]], ['<='], [freq])
			##q>=30 and coverage >=5
			filtering_annotated.filter(working_dir, "and", "1.temp", "2.temp", cols_to_filter[4], ['>=', '>='], [30,5])
			##exonic_variants in refGene
			filtering_annotated.filter(working_dir, "or", "2.temp", outfile_prefix + '.' + str(freq) +  "3.temp", [cols_to_filter[0], cols_to_filter[0]], ['==','=='], [exonic_definitions[0],exonic_definitions[1]])
			##get all protein changing
			filtering_annotated.filter(working_dir, "and", outfile_prefix + '.' + str(freq) +  "3.temp", outfile_prefix+ '.' + sample + '.' + str(freq) +  ".protein_changing.xls", make_list(cols_to_filter[2], non_protein_changing_definitions), make_list('!=', non_protein_changing_definitions), non_protein_changing_definitions)

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

##make avinput files


##run methods
##make star index files for hg19 (only need to do once)
# make_star_index_files(star_index_dir, fa_file, gtf_file, threads)


##for bri batch
##setup working directory where results will be
# working_dir = '/data/atimms/rich_rnaseq/asthma_rnaseq_0217'
# os.chdir(working_dir)
# project_name = 'asthma_bri_0217'
##combine fq files, amp and proccess bam files
# combine_fq_file(bri_samples)
# star_align_paired_end_2_pass(bri_samples, star_index_dir, threads, ['_1.fastq.gz', '_2.fastq.gz'])
# proccess_bam_files(bri_samples)
##variant calling
# vcfs_to_combine = variant_calling_freebayes(bri_samples, proccessed_bam_suffix, '.freebayes.vcf')
# combine_vcf_files_gatk(vcfs_to_combine, project_name + '.freebayes.vcf')
# vcfs_to_combine = variant_calling_gatk_hc(bri_samples, proccessed_bam_suffix, '.gatkHC.vcf')	
# combine_vcf_files_gatk(vcfs_to_combine, project_name + '.gatkHC.vcf')
##intersect vcfs
# intesect_two_vcf_files(project_name + '.gatkHC.vcf', project_name + '.freebayes.vcf', project_name + '.intersected_vcfs')
##annotate vcfs and then filter
# annotate_vcf_file(project_name + '.intersected_vcfs/0002.vcf', project_name, bri_samples)
# filter_ann_txt_files(bri_samples, filtering_cols, frequencies, project_name)
##remove intermediate files


##for bri batch
##setup working directory where results will be
# working_dir = '/data/atimms/rich_rnaseq/asthma_rnaseq_0217'
working_dir = '/data/atimms/rich_rnaseq/amgen_asthma_rnaseq_0317'
os.chdir(working_dir)
project_name = 'asthma_amgen_0417'
##make dict from sample and fq file

amgen_sample_fq_dict = make_sample_fq_dict(amgen_sample_info)
amgen_samples = amgen_sample_fq_dict.keys()
print amgen_sample_fq_dict
print len(amgen_samples), '\n', amgen_samples

# amgen_samples = ['T319', 'T117', 'T114', 'T115', 'T111', 'T119']
# amgen_samples = ['T140']

##combine fq files, amp and proccess bam files
# combine_amgen_fq_files(amgen_sample_fq_dict)
# star_align_paired_end_2_pass(amgen_samples, star_index_dir, threads, ['_1.fq.gz', '_2.fq.gz'])
# proccess_bam_files(amgen_samples)
##variant calling
vcfs_to_combine = variant_calling_freebayes(amgen_samples, proccessed_bam_suffix, '.freebayes.vcf')
print len(vcfs_to_combine)
combine_vcf_files_gatk(vcfs_to_combine, project_name + '.freebayes.vcf')
# vcfs_to_combine = variant_calling_gatk_hc(amgen_samples, proccessed_bam_suffix, '.gatkHC.vcf')	
# vcfs_to_combine = ['N135.gatkHC.vcf', 'T104.gatkHC.vcf', 'T108.gatkHC.vcf', 'T111.gatkHC.vcf', 'T114.gatkHC.vcf', 'T115.gatkHC.vcf', 'T117.gatkHC.vcf', 'T119.gatkHC.vcf', 'T122.gatkHC.vcf', 'T124.gatkHC.vcf', 'T125.gatkHC.vcf', 'T126.gatkHC.vcf', 'T127.gatkHC.vcf', 'T128.gatkHC.vcf', 'T130.gatkHC.vcf', 'T131.gatkHC.vcf', 'T132.gatkHC.vcf', 'T133.gatkHC.vcf', 'T134.gatkHC.vcf', 'T136.gatkHC.vcf', 'T137.gatkHC.vcf', 'T140.gatkHC.vcf', 'T143.gatkHC.vcf', 'T144.gatkHC.vcf', 'T145.gatkHC.vcf', 'T216.gatkHC.vcf', 'T219.gatkHC.vcf', 'T222.gatkHC.vcf', 'T223.gatkHC.vcf', 'T303.gatkHC.vcf', 'T312.gatkHC.vcf', 'T313.gatkHC.vcf', 'T314.gatkHC.vcf', 'T319.gatkHC.vcf', 'T320.gatkHC.vcf', 'T321.gatkHC.vcf', 'T323.gatkHC.vcf', 'T324.gatkHC.vcf', 'T325.gatkHC.vcf', 'T329.gatkHC.vcf', 'T332.gatkHC.vcf']
# print len(vcfs_to_combine)
# combine_vcf_files_gatk(vcfs_to_combine, project_name + '.gatkHC.vcf')
##intersect vcfs
intesect_two_vcf_files(project_name + '.gatkHC.vcf', project_name + '.freebayes.vcf', project_name + '.intersected_vcfs')
##annotate vcfs and then filter
annotate_vcf_file(project_name + '.intersected_vcfs/0002.vcf', project_name, amgen_samples)
filter_ann_txt_files(amgen_samples, filtering_cols, frequencies, project_name)
##remove intermediate files












