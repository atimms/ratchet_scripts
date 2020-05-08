#!/usr/bin/env python
import subprocess
import os
import glob
import shutil
import filtering_annotated
import homozygosity_mapping_ub26
import math

##parameters
delim = '\t'
thread_number = '6'

##working dir
# working_dir = '/data/atimms/timon_0317'
working_dir = '/data/atimms/laura_cell_line_0707'
fastq_dir = '/data/atimms/laura_cell_line_0707'
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
gatk = '/home/atimms/programs/gatk_3.6/GenomeAnalysisTK.jar'
##files etc
fq_dict = {'tumor': ['Sample_1_S1_R1_001.fastq.gz', 'Sample_1_S1_R2_001.fastq.gz'], 
		'spleen': ['Sample_2_S2_R1_001.fastq.gz', 'Sample_2_S2_R2_001.fastq.gz']}
post_bwa_bam = '.bwa.bam'
sorted_bam = '.bwa_sorted.bam'
mkdup_bam = '.bwa_mkdup.bam'
bamslist_file = 'bams.list'
# st_vcf = 'timon_0317.vcf.gz'
# st_avinputs = []
delly_exclude_regions = '/data/atimms/references/annovar/mm10/mouse.mm10.excl.tsv'
bedtools = '/home/atimms/programs/bedtools2-master/bin/bedtools'
gsea = '/home/atimms/programs/GSEA/gsea-3.0.jar'
human_mouse_file = '/data/atimms/references/human_mouse_hcop_fifteen_column_0517.txt'


##annovar parameters
av_genome = 'mm10'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/data/atimms/references/annovar/' + av_genome]
##filter on public dbs
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,snp142,generic,generic,generic', '-genericdbfile', 'mgp.v4.snps.dbSNP.avinput,mgp.v4.indels.dbSNP.avinput,mgp.v5.merged.snps_all.dbSNP142.avinput']
av_operation = ['-operation', 'g,r,r,f,f,f,f']
# av_options = ['-otherinfo', '-remove','-arg', '-splicing 10 ,,,,,,']
av_options = ['-otherinfo', '-remove', '-vcfinput']


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


def combine_fq_file(fq_dir, sample, final_dir):
	print sample
	r1_fq = final_dir + '/' + sample + '.r1.fq.gz'
	r2_fq = final_dir + '/' + sample + '.r2.fq.gz'
	r1_to_combine = glob.glob(fq_dir + '/K*_1.fq.gz')
	r2_to_combine = glob.glob(fq_dir + '/K*_2.fq.gz')
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


def variant_calling_paired_mutect2(bams, name_prefix):
	tumor_bam = bams[0]
	normal_bam = bams[1]
	vcf_temp1 = name_prefix + '.temp_1.vcf'
	vcf_temp2 = name_prefix + '.temp_2.vcf.gz'
	final_vcf = name_prefix + '.mutect.vcf.gz'
	##run mutect caller
	call_mutect = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'MuTect2', '-R', fasta, '-I:tumor', tumor_bam, '-I:normal', normal_bam, '-o', vcf_temp1])
	call_mutect.wait()
	bgzip_run = subprocess.Popen(['bgzip', vcf_temp1])
	bgzip_run.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines, and left normalize indels
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-O', 'z', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', final_vcf])
	bcf_index.wait()

def filter_by_pass(infile, outfile):
	if infile.endswith('.gz'):
		bgzip_run = subprocess.Popen(['bgzip', '-d', infile])
		bgzip_run.wait()
		infile = infile.rsplit('.',1)[0]
	snp_cut = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'SelectVariants', '-R',fasta, '-nt', '15', '-V', infile, '-o', outfile, '--excludeFiltered'])
	snp_cut.wait()

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

def run_table_annovar_vcf(vcf):
	av_prefix = vcf.rsplit('.',2)[0]
	command = [table_annovar] + av_buildver + [vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', av_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()

def multianno_to_annotated(avinputs): 
	head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 'genomicSuperDups', 'snp142','mgp.v4.snps','mgp.v4.indels','mgp.v5.snps', 'Zygosity', 'Qual', 'Coverage', 'Chr', 'Pos', 'Filter', 'Ref2', 'Alt2','Qual', 'Filter', 'GT_info', 'Format', 'spleen', 'tumor']
	head_out = delim.join(head + ['\n'])
	for avinput in avinputs:
		av_prefix = avinput.rsplit('.',2)[0]
		multianno = av_prefix + '.mm10_multianno.txt'
		annotated = av_prefix + '.annotated.xls'
		with open(multianno, "r") as multi, open(annotated, "w") as final:
			final.write(head_out)
			line_count = 0
			for line in multi:
				line_count += 1
				if line_count > 1:
					final.write(line)


def make_bed_from_ann_txt(infile, outfile):
	with open(outfile, "w") as out_fh, open(infile, "r") as in_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			if line_count >1:
				line = line.rstrip().split(delim)
				out_fh.write(delim.join(line[:3] + ['\n']))

def filter_exonic_variants(infile, outfile):
	col_exon = 6
	exon_definition = ['exonic', 'splicing']
	col_function = 9
	syn_definition = 'synonymous SNV'
	filtering_annotated.filter(working_dir, "or", infile , 'temp1.txt', [col_exon, col_exon], ['==','=='], [exon_definition[0],exon_definition[1]])
	##remove synonymous
	filtering_annotated.filter(working_dir, "and", 'temp1.txt', outfile, [col_function], ['!='], [syn_definition])

def add_rnaseq_data(infile, outfile, rnaseq_file):
	rnaseq_dict = {}
	with open(rnaseq_file, "r") as rna_fh:
		lc = 0
		for line in rna_fh:
			lc += 1
			if lc >= 1:
				line = line.rstrip().split(',')
				gene = line[0].replace('"', '')
				logfc = line[2]
				padj = line[6]
				if gene in rnaseq_dict:
					print 'gene %s seen twice'%gene
				else:
					rnaseq_dict[gene] = [logfc, padj]
	with open(outfile, "w") as out_fh, open(infile, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.rstrip().split(delim)
			gene = line[6]
			if lc == 1:
				header = line[:9] + ['log2FoldChange', 'padj'] + line[9:] + ['\n']
				out_fh.write(delim.join(header))
			else:
				if gene in rnaseq_dict:
					line_out = line[:9] + rnaseq_dict[gene] + line[9:] + ['\n']
				else:
					line_out = line[:9] + ['na', 'na'] + line[9:] + ['\n']
				out_fh.write(delim.join(line_out))



def rank_deseq_results_for_gsea_from_geo_microarray(deseq_de_results, ranked_de_results):
	gene_list = []
	with open(deseq_de_results, "r") as inph, open(ranked_de_results, "w") as outph:
		outph.write(delim.join(['gene', 'signed_log10_pvalue', '\n']))
		line_count = 0
		for line in inph:
			line_count += 1
			line = line.strip().split(delim)
			if line_count > 1:
				gene = line[6].translate(None, '"').split("///")[0]
				# print line
				if gene in gene_list:
					print 'gene %s seen multiple times'%gene
				else:
					if gene != '':
						gene_list.append(gene)
						log2fc = line[5].translate(None, '"')
						pvalue = line[1].translate(None, '"')
						if pvalue != 'NA':
							if pvalue == '0':
								pvalue = '1e-300'
							# print gene, log2fc, pvalue
							log10_p = abs(math.log10(float(pvalue)))
							if float(log2fc) >= 0:
								final = log10_p/1
							else:
								final = log10_p/-1
							# print gene, log2fc, pvalue, log10_p, final
							outph.write(delim.join([gene, str(final), '\n']))


def rank_deseq_results_for_gsea(deseq_de_results, ranked_de_results):
	with open(deseq_de_results, "r") as inph, open(ranked_de_results, "w") as outph:
		outph.write(delim.join(['gene', 'signed_log10_pvalue', '\n']))
		line_count = 0
		for line in inph:
			line_count += 1
			line = line.strip().split(',')
			if line_count > 1:
				gene = line[0].translate(None, '"')
				log2fc = line[2]
				pvalue = line[5]
				if pvalue != 'NA':
					if pvalue == '0':
						pvalue = '1e-300'
					# print gene, log2fc, pvalue
					log10_p = abs(math.log10(float(pvalue)))
					if float(log2fc) >= 0:
						final = log10_p/1
					else:
						final = log10_p/-1
					print gene, log2fc, pvalue, log10_p, final
					outph.write(delim.join([gene, str(final), '\n']))

def convert_mouse_to_human_ranked_file(infile, outfile, hcop_file):
	##make dict from mouse and human gene names
	gename_dict = {}
	with open(hcop_file, 'r') as hcop_fh:
		line_count = 0
		for line in hcop_fh:
			line_count += 1
			line = line.rstrip().split(delim)
			if line_count >1:
				human_symbol = line[4]
				mouse_symbol = line[11]
				# print human_symbol, mouse_symbol
				if mouse_symbol in gename_dict:
					gename_dict[mouse_symbol].append(human_symbol)
				else:
					gename_dict[mouse_symbol] = [human_symbol]
	##check dict
	# for g in gename_dict:
	# 	print g, gename_dict[g]
	##update file
	# '''
	with open(infile, 'r') as in_fh, open(outfile, 'w') as out_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			if line_count == 1:
				out_fh.write(line)
			else:
				line = line.rstrip().split(delim)
				m_gene = line[0]
				if m_gene in gename_dict:
					h_gene = gename_dict[m_gene][0]
					signed_value = line[1]
					# print h_gene, m_gene
					out_fh.write(delim.join([h_gene,signed_value]) + '\n')
				else:
					print 'gene %s not in file'%m_gene
	# '''

def run_gsea(ranked_file, gene_sets, test_id, sets_to_plot, geneset_min_max):
	##java -Xmx512m xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c2.cp.biocarta.v5.2.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk /Users/atimms/gsea_home/ranked_files/greb1l_rnaseq_0817.wt_ko.de.cap.rnk -scoring_scheme weighted -rpt_label test -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /Users/atimms/gsea_home/output/sep21 -gui false
	##java -cp ~/programs/gsea-3.0.jar -Xmx512m xtools.gsea.GseaPreranked
	run_gsea = subprocess.Popen(['java', '-cp', gsea, '-Xmx5000m', 'xtools.gsea.GseaPreranked', '-gmx', gene_sets,
			'-norm', 'meandiv', '-nperm', '10000', '-rnk', ranked_file, '-scoring_scheme', 'weighted', '-rpt_label', test_id,
			'-make_sets', 'true', '-plot_top_x', sets_to_plot, '-rnd_seed', 'timestamp' '-set_max', geneset_min_max[1],
			'-set_min', geneset_min_max[0], '-zip_report', 'false', '-out', working_dir, '-gui', 'false'])
	run_gsea.wait()

def get_sig_genes_write_log(pos_res, neg_res, res_dir, logfile, fdr_req):
	sig_genes = []
	with open(logfile, 'a') as log_fh:
		log_fh.write('\n' + 'for results' + res_dir + ':\n')
		print 'for results in dir: ' + res_dir
		lc = 0
		##get sig genes sets
		sig_gene_sets = []
		with open(pos_res, 'r') as pos_fh:
			for line in pos_fh:
				lc += 1
				if lc >1:
					line = line.rstrip().split(delim)
					fdr = float(line[7])
					gset = line[0]
					if fdr <= fdr_req:
						log_fh.write(gset + ' is positively enriched with an fdr of ' + str(fdr) + '\n')
						print gset + ' is positively enriched with an fdr of ' + str(fdr)
						sig_gene_sets.append(gset)
		lc = 0
		with open(neg_res, 'r') as neg_fh:
			for line in neg_fh:
				lc += 1
				if lc >1:
					line = line.rstrip().split(delim)
					fdr = float(line[7])
					gset = line[0]
					if fdr <= fdr_req:
						log_fh.write(gset + ' is negatively enriched with an fdr of ' + str(fdr) + '\n')
						print gset + ' is negatively enriched with an fdr of ' + str(fdr)
						sig_gene_sets.append(gset)
		##get genes for all sig gene sets
		for sig_gene_set in sig_gene_sets:
			log_fh.write('\n' + 'in gene set ' + sig_gene_set + ':\n')
			log_fh.write('we have the following genes:\n')
			infile = res_dir + '/' + sig_gene_set + '.xls'
			with open(infile, 'r') as in_fh:
				for line in in_fh:
					lc += 1
					if lc >1:
						line = line.rstrip().split(delim)
						gene = line[1]
						sig_genes.append(gene)
						log_fh.write(gene + '\n')
		return sig_genes

def convert_human_mouse_genelist(human_genelist, hcop_file):
	##make dict from mouse and human gene names
	gename_dict = {}
	with open(hcop_file, 'r') as hcop_fh:
		line_count = 0
		for line in hcop_fh:
			line_count += 1
			line = line.rstrip().split(delim)
			if line_count >1:
				human_symbol = line[4]
				mouse_symbol = line[11]
				# print human_symbol, mouse_symbol
				if human_symbol in gename_dict:
					gename_dict[human_symbol].append(mouse_symbol)
				else:
					gename_dict[human_symbol] = [mouse_symbol]
	mouse_genelist = []
	for hgene in human_genelist:
		if hgene in gename_dict:
			mgene = gename_dict[hgene][0]
			# print h_gene, m_gene
			mouse_genelist.append(mgene)
		else:
			print 'gene %s not in file'%hgene
	return mouse_genelist

def add_gsea_data_to_var_file(infile, outfile, sig_genes):
	with open(outfile, "w") as out_fh, open(infile, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.rstrip().split(delim)
			gene = line[6]
			if lc == 1:
				header = line[:9] + ['in_gsea_pathway'] + line[9:] + ['\n']
				out_fh.write(delim.join(header))
			else:
				if gene in sig_genes:
					line_out = line[:9] + ['yes'] + line[9:] + ['\n']
				else:
					line_out = line[:9] + ['no'] + line[9:] + ['\n']
				out_fh.write(delim.join(line_out))

def add_sig_gene_sets_to_variant_results(infile, outfile_prefix, results_dirs, log_prefix, fdrs_req):
	for fdr_req in fdrs_req:
		logfile = log_prefix + '.fdr_' + str(fdr_req) + '.txt'
		outfile = outfile_prefix + '.' + str(fdr_req) + '.xls'
		all_genes = []
		##get all genes in significantly altered pathways
		for results_dir in results_dirs:
			pos_results = results_dir + '/gsea_report_for_na_pos_' + results_dir.rsplit('.', 1)[1] + '.xls'
			neg_results = results_dir + '/gsea_report_for_na_neg_' + results_dir.rsplit('.', 1)[1] + '.xls'
			genes = get_sig_genes_write_log(pos_results, neg_results, results_dir, logfile, fdr_req)
			all_genes.extend(genes)
		print fdr_req, len(list(set(all_genes)))
		##convert human genes to mouse
		mouse_genes = convert_human_mouse_genelist(list(set(all_genes)), human_mouse_file)
		print len(list(set(mouse_genes)))
		##add info to annotated file
		add_gsea_data_to_var_file(infile, outfile, mouse_genes)



##call methods, samtools, annovar and format multianno



##map with bwa and process with samtools etc
project = 'laura_0717'
bam_files = ['tumor.bwa_mkdup.bam', 'spleen.bwa_mkdup.bam']
mutect_vcf = 'laura_0717.mutect.vcf.gz'
mutect_filtered_vcf = 'laura_0717.mutect.filtered.vcf'
st_avinputs = ['TUMOR.avinput']
annoated_file = 'laura_0717.mutect.annotated.xls'
exonic_annoated_file = 'laura_0717.mutect.annotated.exonic.xls'
exonic_rna_annoated_file = 'laura_0717.mutect.annotated.exonic.rnaseq.xls'
rnaseq_data = 'laura_0417.results.tumor_vs_normal.csv'
# align_with_bwa(fq_dict)
# variant_calling_paired_mutect2(bam_files, project)
# filter_by_pass(mutect_vcf, mutect_filtered_vcf)
# run_table_annovar(st_avinputs)
# run_table_annovar_vcf(mutect_filtered_vcf)
# multianno_to_annotated([mutect_filtered_vcf])
# filter_exonic_variants(annoated_file, exonic_annoated_file)
# add_rnaseq_data(exonic_annoated_file, exonic_rna_annoated_file, rnaseq_data)

##to do - compare genes against rnaseq and microarray
##get GSEA running on 26
##run gsea on microarray data x2 -- make rnk file
##run gsea on rnaseq data x2 -- make rnk file and use mouse
##annoate vars with if gene in a enriched pathway, list direction and pathway name

##parameters
kegg_gene_sets = 'gseaftp.broadinstitute.org://pub/gsea/gene_sets/c2.cp.kegg.v5.2.symbols.gmt'
biocarta_gene_sets = 'gseaftp.broadinstitute.org://pub/gsea/gene_sets/c2.cp.biocarta.v5.2.symbols.gmt'

##human micrarray data
name = 'gse1987_scc_all_vs_normal'
deseq_de_results = name + '.txt'
ranked_de_results = name + '.rnk'
# rank_deseq_results_for_gsea_from_geo_microarray(deseq_de_results, ranked_de_results)
# run_gsea(ranked_de_results, biocarta_gene_sets, name + '.biocarta', '200', ['15','500'])
# run_gsea(ranked_de_results, kegg_gene_sets, name + '.kegg', '200', ['15','500'])



##mouse rnaseq
name = 'laura_0417.results.tumor_vs_normal'
deseq_de_results = name + '.csv'
ranked_de_results = name + '.rnk'
ranked_de_results_human = name + '.human.rnk'
# rank_deseq_results_for_gsea(deseq_de_results, ranked_de_results)
# convert_mouse_to_human_ranked_file(ranked_de_results, ranked_de_results_human, human_mouse_file)
# run_gsea(ranked_de_results_human, biocarta_gene_sets, name + '.biocarta', '200', ['15','500'])
# run_gsea(ranked_de_results_human, kegg_gene_sets, name + '.kegg', '200', ['15','500'])


##
results_folders = ['gse1987_scc_all_vs_normal.biocarta.GseaPreranked.1506372707240', 'gse1987_scc_all_vs_normal.kegg.GseaPreranked.1506372813720',
			'laura_0417.results.tumor_vs_normal.biocarta.GseaPreranked.1506372968835', 'laura_0417.results.tumor_vs_normal.kegg.GseaPreranked.1506373070918']

in_variant_file = 'laura_0717.mutect.annotated.exonic.rnaseq.xls'
out_variant_file_prefix = 'laura_0717.mutect.annotated.exonic.rnaseq_and_pathways'
log_file_prefix = 'laura_0717.added_gsea_genes_sets.0917'
add_sig_gene_sets_to_variant_results(in_variant_file, out_variant_file_prefix, results_folders, log_file_prefix, [0.01, 0.05, 0.1])




