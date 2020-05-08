#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import filtering_annotated
import shutil

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
genome_name = 'mm10'
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
imprinted_bed = 'imprinted_genes_merged.bed'
mono_bed = 'mono_alleleic_genes_merged.bed'
sig_imprinted_bed = 'sig_imprinted_merged.bed'

def combine_fqs(fqs_to_combine, fq_name):
	with open(fq_name, 'w') as cfq_fh:
		print "combining %s files named %s to make file %s"%(len(fqs_to_combine), fqs_to_combine, fq_name)
		cat_files = subprocess.Popen(['cat'] + fqs_to_combine, stdout=cfq_fh)
		cat_files.wait()


def star_align_paired_end_2_pass(sample_name, star_genome_dir, threads_to_use, fq_prefixes):
	r1_fq = sample_name + fq_prefixes[0]
	r2_fq = sample_name + fq_prefixes[1]
	print sample_name, r1_fq, r2_fq
	star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', r1_fq, r2_fq, '--outSAMtype', 'BAM', 'Unsorted', '--outSAMmapqUnique', '60', '--readFilesCommand', 'zcat', '--twopassMode', 'Basic', '--outFileNamePrefix', sample_name + '.', '--runThreadN', threads_to_use])
	star_align.wait()



def proccess_bam_files(sample_name):
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

def add_rg(sample_name):
	in_bam = sample_name + star_bam_suffix
	arg_bam = sample_name + '.with_rg.bam'
	picard_arg = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'AddOrReplaceReadGroups', 'INPUT=' + in_bam, 'OUTPUT=' + arg_bam, 'SO=coordinate', 'RGID=' + sample_name, 'RGLB=' + sample_name, 'RGPL=Illumina', 'RGPU=machine1', 'RGSM=' + sample_name, 'CREATE_INDEX=true' ])
	picard_arg.wait()	



##make list of all bam files to be analyzed
def make_list_of_bams(sample_dict, bam_suffix, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for sample in sample_dict:
			bam = sample + bam_suffix
			outfh.write(bam + '\n')



def variant_calling_gatk_hc(bamlist, name_prefix, vcf_suffix, gene_bed):
	#java -Xmx32g -jar gatk -T UnifiedGenotyper -R hg19_fasta -rf BadCigar -allowPotentiallyMisencodedQuals -L MIPtargets.intervals -I bam.realigned.list -o vcf -dcov 5000 -dt NONE -glm both -A AlleleBalance -stand_emit_conf 10 -minIndelFrac 0.005 -G Standard
	#java -Xmx8g -jar gatk -T VariantFiltration -R hg19_fasta -rf BadCigar -allowPotentiallyMisencodedQuals -V in_vcf -o out.vcf -window 20 -cluster 5 -filterName ABFilter -filter "ABHet > 0.75" -filterName QDFilter -filter "QD < 5.0" -filterName QUALFilter -filter "QUAL < 30.0" -filterName LowCoverage -filter "DP < 5"
	vcf_files = []
	vcf_temp1 = name_prefix + 'temp_ug1' + vcf_suffix
	vcf_temp2 = name_prefix + 'temp_ug2' + vcf_suffix
	final_vcf = name_prefix + vcf_suffix
	vcf_files.append(final_vcf)
	# gatk_hc = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fa_file, '-I', bamlist, '-L', gene_bed, '-o', vcf_temp1, '-dontUseSoftClippedBases', '-stand_call_conf', '20.0', '-stand_emit_conf', '20.0'])
	gatk_hc = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fa_file, '-I', bamlist, '-L', gene_bed, '-o', vcf_temp1, '-stand_call_conf', '20.0', '-stand_emit_conf', '20.0'])

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

def variant_calling_gatk_hc_star_bams(bamlist, name_prefix, vcf_suffix, gene_bed):
	#java -Xmx32g -jar gatk -T UnifiedGenotyper -R hg19_fasta -rf BadCigar -allowPotentiallyMisencodedQuals -L MIPtargets.intervals -I bam.realigned.list -o vcf -dcov 5000 -dt NONE -glm both -A AlleleBalance -stand_emit_conf 10 -minIndelFrac 0.005 -G Standard
	#java -Xmx8g -jar gatk -T VariantFiltration -R hg19_fasta -rf BadCigar -allowPotentiallyMisencodedQuals -V in_vcf -o out.vcf -window 20 -cluster 5 -filterName ABFilter -filter "ABHet > 0.75" -filterName QDFilter -filter "QD < 5.0" -filterName QUALFilter -filter "QUAL < 30.0" -filterName LowCoverage -filter "DP < 5"
	vcf_files = []
	vcf_temp1 = name_prefix + 'temp_ug1' + vcf_suffix
	vcf_temp2 = name_prefix + 'temp_ug2' + vcf_suffix
	final_vcf = name_prefix + vcf_suffix
	vcf_files.append(final_vcf)
	# gatk_hc = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fa_file, '-I', bamlist, '-L', gene_bed, '-o', vcf_temp1, '-dontUseSoftClippedBases', '-stand_call_conf', '20.0', '-stand_emit_conf', '20.0'])
	gatk_hc = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fa_file, '-I', bamlist, '-L', gene_bed, '-o', vcf_temp1, '-stand_call_conf', '20.0', '-stand_emit_conf', '20.0', '--filter_reads_with_N_cigar', '-dontUseSoftClippedBases'])

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



def filter_vars(name_prefix, vcf_suffix, out_suffix):
	vcf = name_prefix + vcf_suffix
	outfile = name_prefix + out_suffix
	vars_picked = []
	with open(vcf, 'r') as vcf_fh, open(outfile, 'w') as out_fh:
		for line in vcf_fh:
			if line[0:2] != '##':
				if line[0] == '#':
					out_fh.write(line)
				else:
					line = line.rstrip().split(delim)
					genotypes = line[9:]
					var = line[:5]
					filter_info = line[6]
					ref_len = len(line[3])
					alt_len = len(line[4])
					qual = float(line[5])
					if filter_info == 'PASS' and ref_len == 1 and alt_len == 1 and qual > 50:
						for genotype in genotypes:
							# print genotype
							gt = genotype.split(':')[0]
							if gt == '0/1':
								cov = genotype.split(':')[2]
								if int(cov) >= 20:
									if var not in vars_picked:
										out_fh.write(delim.join(line) + '\n')
										vars_picked.append(var)


##run methods
sample_fq_dict = {'C7_D0': [['C7_D0_1_1.fq.gz', 'C7_D0_2_1.fq.gz', 'C7_D0_3_1.fq.gz'], ['C7_D0_1_2.fq.gz', 'C7_D0_2_2.fq.gz', 'C7_D0_3_2.fq.gz']],
				'C7_D4': [['C7_D4_1_1.fq.gz', 'C7_D4_2_1.fq.gz', 'C7_D4_3_1.fq.gz'], ['C7_D4_1_2.fq.gz', 'C7_D4_2_2.fq.gz', 'C7_D4_3_2.fq.gz']],
				'C7_D8': [['C7_1_1.fq.gz', 'C7_2_1.fq.gz', 'C7_3_1.fq.gz'], ['C7_1_2.fq.gz', 'C7_2_2.fq.gz', 'C7_3_2.fq.gz']],
				'F1_D0': [['F1_D0_1_1.fq.gz', 'F1_D0_2_1.fq.gz', 'F1_D0_3_1.fq.gz'], ['F1_D0_1_2.fq.gz', 'F1_D0_2_2.fq.gz', 'F1_D0_3_2.fq.gz']],
				'F1_D4': [['F1_D4_1_1.fq.gz', 'F1_D4_2_1.fq.gz', 'F1_D4_3_1.fq.gz'], ['F1_D4_1_2.fq.gz', 'F1_D4_2_2.fq.gz', 'F1_D4_3_2.fq.gz']],
				'F1_D8': [['F1_1_1.fq.gz', 'F1_2_1.fq.gz', 'F1_3_1.fq.gz'], ['F1_1_2.fq.gz', 'F1_2_2.fq.gz', 'F1_3_2.fq.gz']]}

##run methods


##combine fq files, map and process bams
for sample in sample_fq_dict:
	r1_fq = sample + '.r1.fastq.gz'
	r2_fq = sample + '.r2.fastq.gz'
	##combine the 3 sample for each timepoint
	# combine_fqs(sample_fq_dict[sample][0], r1_fq)
	# combine_fqs(sample_fq_dict[sample][1], r2_fq)
	# star_align_paired_end_2_pass(sample, star_index_dir, threads, ['.r1.fastq.gz', '.r2.fastq.gz'])
	# proccess_bam_files(sample)
	##just add the read group
	# add_rg(sample)

##make bams.list
# make_list_of_bams(sample_fq_dict.keys(), proccessed_bam_suffix, bamlist)
make_list_of_bams(sample_fq_dict.keys(), '.with_rg.bam', bamlist)


##run gatk on imprinted gene
# project_name = 'greb1l_eb_imprinted_vars_1017'

# variant_calling_gatk_hc(bamlist,project_name, '.gatkHC.vcf', imprinted_bed)
# filter_vars(project_name, '.gatkHC.vcf', '.het_vars.xls')
##repeat without the soft clipped bases setting
# variant_calling_gatk_hc(bamlist,project_name, '.gatkHC_2.vcf', imprinted_bed)
# filter_vars(project_name, '.gatkHC_2.vcf', '.het_vars_2.xls')
##from star bams
# variant_calling_gatk_hc_star_bams(bamlist,project_name, '.gatkHC_3.vcf', imprinted_bed)
# filter_vars(project_name, '.gatkHC_3.vcf', '.het_vars_3.xls')


##run gatk on mono allele gene
##get bed file.. in xl use genes list to get cds start and end, then sort/merge with bedtools
# project_name = 'greb1l_eb_mono_alleleic_vars_1017'
# variant_calling_gatk_hc(bamlist,project_name, '.gatkHC.vcf', mono_bed)
# filter_vars(project_name, '.gatkHC.vcf', '.het_vars.xls')
##repeat without the soft clipped bases setting
# variant_calling_gatk_hc(bamlist,project_name, '.gatkHC_2.vcf', mono_bed)
# filter_vars(project_name, '.gatkHC_2.vcf', '.het_vars_2.xls')
##from star bams
# variant_calling_gatk_hc_star_bams(bamlist,project_name, '.gatkHC_3.vcf', mono_bed)
# filter_vars(project_name, '.gatkHC_3.vcf', '.het_vars_3.xls')

##reduced list of imprinted genes
project_name = 'greb1l_eb_sig_imprinted_vars_reduced_1017'
variant_calling_gatk_hc_star_bams(bamlist,project_name, '.gatkHC_3.vcf', sig_imprinted_bed)
filter_vars(project_name, '.gatkHC_3.vcf', '.het_vars_3.xls')



