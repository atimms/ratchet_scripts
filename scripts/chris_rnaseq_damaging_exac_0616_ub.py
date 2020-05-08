#!/tools/BioBuilds-2015.04/bin/python
import filtering_annotated
import os
import subprocess
import glob
import shutil
# from scipy import stats

##set input variables and parameters
delim = '\t'
##working directory
working_dir = '/data/atimms/chris_rnaseq'
os.chdir(working_dir)
# css_genes_prefix = 'esp6500_css_genes_0815'
# all_prefix = 'esp6500_all_0915'
exac_prefix = 'exac3_1115'
exac_unfiltered_prefix  = exac_prefix + '_unfiltered'
##cases and controls lists
cases_list = ['95445', '95446', '95447', '95448', '95450', '95452', '95453', '95454', '95456', '95457', '95458', '95459', '95460', '95462', '95463', '95464', '95466', '95467', '95468', '95469', '95470', '95471', '95472', '95474', '95475', '95476', '95477', '95479', '95481', '95482', '95484', '95485', '95486', '95487', '95488', '95490', '95491', '95492', '95493', '95494', '95495', '95496', '95498', '95500', '95501', '95502', '95505', '95506', '95507', '95509', '95510', '95511', '95512', '95514', '95515', '95516', '95518', '95519', '95520', '95521', '95522', '95524', '95525', '95526', '95527', '95528', '95529', '95530', '95531', '95534', '95535', '95536', '95537', '95540', '95542', '95543', '95544', '95545', '95546', '95547', '95549', '95550', '95551', '95552', '95555', '95556', '95557', '95558', '95559', '95560', '95562', '95563', '95564', '95565', '95566', '95567', '95569', '95570', '95571', '95572', '95573', '95575', '95576', '95578', '95579', '95580', '95581', '95583', '95584', '95585', '95586', '95587', '95588', '95589', '95591', '95592', '95593', '95594', '95595', '95596', '95597', '95598', '95599', '95600', '95601', '95602', '95603', '95604', '95605', '95606', '95609', '95610', '95611', '95612', '95614', '95615', '95616', '95617', '95619', '95620', '95622', '95623', '95624', '95626', '95627', '95628', '95629', '95630', '95631', '95633', '95634', '95635', '95636', '95637', '95638', '95639', '95640', '95642', '95644', '95645', '95646', '95649', '95650', '95651', '95652', '95653', '95654', '95655', '95656', '95657', '95658', '95661', '95662', '95663', '95665', '95668', '95669', '95670', '95671', '95673', '95675', '95676', '95677', '95678', '95679', '95680', '95681', '95682', '95683', '95684', '95685', '95686', '95688', '95689', '95690', '95691', '95693', '95694', '95696', '95697', '95698', '95699', '95700', '95702', '95703']
controls_list = ['95449', '95451', '95455', '95461', '95465', '95473', '95478', '95480', '95483', '95489', '95497', '95499', '95504', '95508', '95513', '95517', '95523', '95532', '95533', '95538', '95539', '95541', '95548', '95553', '95554', '95561', '95568', '95574', '95577', '95582', '95590', '95607', '95608', '95613', '95618', '95621', '95625', '95647', '95648', '95659', '95660', '95666', '95667', '95672', '95674', '95692', '95695', '95704']

##file names etc
# refgene_file = '/Users/atimms/Desktop/ngs/references/annovar/hg19/hg19_refGene.txt'
# gene_regions_file = 'css_gene_regions.txt'
# gene_regions_vcf = 'esp6500.css_genes.0815.vcf.gz'
# project_avinput = css_genes_prefix + '.avinput'
cases_annotated = 'css_cases_0815.annotated.txt'
cases_passed_annotated = 'css_cases_0815.gatk_passed.annotated.txt'
controls_passed_annotated = 'css_controls_0815.gatk_passed.annotated.txt'
exac_vcf = 'ExAC.r0.3.sites.vep.tidy.vcf.gz'
exac_passed_vcf = 'ExAC.r0.3.sites.vep.tidy.passed.vcf.gz'
exac_avinput =  exac_prefix + '.avinput'
exac_annotated = exac_prefix + '.annotated.txt'
exac_unfiltered_annotated = exac_unfiltered_prefix + '.annotated.txt'
fasta = '/data/atimms/references/human_g1k_v37.fasta'
##oufiles
rvis_lite_exac_out_file_suffix = '.exac.rvis_lite.xls'

##annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/data/atimms/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,knownGene,ensGene,rmsk,ljb26_all,esp6500si_all,esp6500si_aa,esp6500si_ea,exac02,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_amr,1000g2014oct_eas,1000g2014oct_eur,1000g2014oct_sas,snp138']
av_operation = ['-operation', 'g,g,g,r,f,f,f,f,f,f,f,f,f,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.']

##programs
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
gatk = '/tools/GATK/GenomeAnalysisTK.jar'


##variables for damaging counts

##definitions for determining var function
exonic_definitions = ['exonic', 'splicing']
disruptive_definitions = ['','.', 'stopgain', 'stoploss','frameshift insertion', 'frameshift deletion', 'frameshift_insertion', 'frameshift_deletion']
protein_changing_definitions =['','.', 'stopgain', 'stoploss', 'nonsynonymous SNV','frameshift insertion', 'frameshift deletion', 'nonsynonymous_SNV','frameshift_insertion', 'frameshift_deletion']
nosyn_definitions =['nonsynonymous SNV','nonsynonymous_SNV']
pp2_score = [0.9, 0.9]
cadd_gerp_score = [15,3]
damaging_definitions = [0.9,0.9,15,3]
no_record_definition = ['', '.']
##list of columns for both methods (all start at 1):
##for filtering variants: refgene function, frequency dbs, refgene exonic function, pp2 x2/cadd/gerp, refgene gename
case_ctl_rvis_col = [16, [57,58,59,60,61,62,63,64,65,66], 19, [34,36,52,53], 17]
esp6500_rvis_col = [6, [47,48,49,50,58,59,60,61,62,63], 9, [24,26,42,43], 7]
case_ctl_qual_cov = [8,9]
esp6500_info_col = 72
# esp6500_col_list = [6, 5, 8, 0, 5, 23, 25, 42, 41, 70]
# case_ctl_col_list = [16, 15, 18, 1, 6, 33, 35, 52, 51, 9]
frequencies = [0.01, 0.005, 0.001]



##values we what >= to call variant damaging: pp2, gerp, cadd
# damaging_numbers_1 = [0.9,3,12]
# all_var_counts = all_prefix + '.pp2_' + str(damaging_numbers_1[0]) + '.gerp_' + str(damaging_numbers_1[1]) + '.cadd_' + str(damaging_numbers_1[2]) + '.var_counts.xls'


##remove not passed snps fom vcf
def remove_not_passed_from_vcf(in_vcf, out_vcf):
	# tabix_vcf = subprocess.Popen(['tabix', in_vcf])
	# tabix_vcf.wait()
	gatk_cut = subprocess.Popen(['java', '-Xmx50g', '-jar', gatk, '-T', 'SelectVariants', '-R', fasta, '-V', in_vcf, '-o', out_vcf, '--excludeFiltered'])
	gatk_cut.wait()


##convert vcf file to individual avinput file
def convert_to_annovar(vcf, av_name):
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', vcf, '-includeinfo', '-outfile', av_name])
	con_ann.wait()

##annoate avinput
def run_table_annovar(avinput, av_prefix):
	command = [table_annovar] + av_buildver + [avinput] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', av_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()
	

def multianno_to_annotated(av_prefix):
	head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'Func.knownGene', 'Gene.knownGene', 'GeneDetail.knownGene', 'ExonicFunc.knownGene', 'AAChange.knownGene', 'Func.ensGene', 'Gene.ensGene', 'GeneDetail.ensGene', 'ExonicFunc.ensGene', 'AAChange.ensGene', 'rmsk', 'SIFT_score', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score', 'Polyphen2_HVAR_pred', 'LRT_score', 'LRT_pred', 'MutationTaster_score', 'MutationTaster_pred', 'MutationAssessor_score', 'MutationAssessor_pred', 'FATHMM_score', 'FATHMM_pred', 'RadialSVM_score', 'RadialSVM_pred', 'LR_score', 'LR_pred', 'VEST3_score', 'CADD_raw', 'CADD_phred', 'GERP++_RS', 'phyloP46way_placental', 'phyloP100way_vertebrate', 'SiPhy_29way_logOdds', 'esp6500si_all', 'esp6500si_aa', 'esp6500si_ea', 'ExAC_Freq', 'ExAC_AFR', 'ExAC_AMR', 'ExAC_EAS', 'ExAC_FIN', 'ExAC_NFE', 'ExAC_OTH', 'ExAC_SAS', '1000g2014oct_all', '1000g2014oct_afr', '1000g2014oct_amr', '1000g2014oct_eas', '1000g2014oct_eur', '1000g2014oct_sas', 'snp138', 'chr', 'pos', 'id', 'ref', 'alt', 'id2', 'filter', 'info']
	head_out = delim.join(head + ['\n'])
	multianno = av_prefix + '.hg19_multianno.txt'
	annotated = av_prefix + '.annotated.txt'
	with open(multianno, "r") as multi, open(annotated, "w") as final:
		final.write(head_out)
		line_count = 0
		for line in multi:
			line_count += 1
			if line_count > 1:
				final.write(line)

##calls all annovar methods
def annotate_vcf_file(vcf_file, av_input, avprefix):
	convert_to_annovar(vcf_file, av_input)
	run_table_annovar(av_input, avprefix)
	multianno_to_annotated(avprefix)

##simple list maker -- give string and dictionary to get length from
def make_list(string, dict):
	l = []
	for i in range(len(dict)):
		l.append(string)
	return l

##method to combine files
def combine_ann_txt(sample_list, prefix, suffix, outfile):
	with open(outfile, "w") as final_file:
		sample_count = 0
		for sample in sample_list:
			ifile = prefix + sample + suffix
			sample_count += 1
			with open(ifile, "r") as open_file:
				line_count = 0
				if sample_count == 1:
					for line in open_file:
						line_count += 1
						if line_count == 1:
							final_file.write('Proband' + delim + line)
						else:
							final_file.write(sample + delim + line)
				else:
					for line in open_file:
						line_count += 1
						if line_count > 1:
							final_file.write(sample + delim + line)

def remove_dup_vars_and_add_count(infile, outfile):
	var_dict = {}
	with open(infile, "r") as in_fh:
		line_count = 0
		for line in in_fh:
			line = line.split(delim)
			line_count += 1
			if line_count == 1:
				header = delim.join(['Count'] + line[1:])
			else:
				var = '.'.join(line[1:6])
				if var in var_dict:
					var_dict[var][0] += 1
				else:
					var_dict[var] = [1] + line[1:]
	with open(outfile, "w") as out_fh:
		out_fh.write(header)
		for variant in var_dict:
			var_dict[variant][0] = str(var_dict[variant][0])
			out_fh.write(delim.join(var_dict[variant]))


def remove_rows_with_no_data(infile, outfile, cols):
	with open(outfile, "w") as out_fh:
		with open(infile, "r") as in_fh:
			line_count = 0
			vars_passed_count = 0
			for line in in_fh:
				line_count += 1
				empty_count = 0
				if line_count == 1:
					out_fh.write(line)
				else:
					line = line.split(delim)
					for col in cols:
						if line[col -1] == '' or line[col -1] == '.':
							empty_count += 1
					if empty_count == 0:
						out_fh.write(delim.join(line))
						vars_passed_count +=1
	print 'filtering variants from the file', infile
	print 'removing variants with empty fields'
	print '%s kept out of %s checked'%(vars_passed_count, line_count -1)
	print 'writen to file:', outfile, '\n'



def filter_ann_txt_files(exac_txt, cases_txt, control_txt, esp6500_cols, case_controls_col, freq_req):
	##remove duplicate vars from cases and controls and add count in first col
	remove_dup_vars_and_add_count(cases_txt,'cases.nodup_temp.txt')
	remove_dup_vars_and_add_count(control_txt,'controls.nodup_temp.txt')
	##filter rnaseq data
	for freq in freq_req:
		for group in ['cases', 'controls']:
			##only 'rare' (<=1%)
			filtering_annotated.filter(working_dir, "and", group + '.nodup_temp.txt', group + '.' + str(freq) +  ".1.temp", case_controls_col[1], make_list('<=', case_controls_col[1]), make_list(freq, case_controls_col[1]))
			##q>=30 and coverage >=5
			filtering_annotated.filter(working_dir, "and", group + '.' + str(freq) + ".1.temp", group + '.' + str(freq) +  ".2.temp", case_ctl_qual_cov, ['>=', '>='], [30,5])
			##exonic_variants in refGene
			filtering_annotated.filter(working_dir, "or", group + '.' + str(freq) +  ".2.temp", group + '.' + str(freq) +  ".exonic.xls", [case_controls_col[0], case_controls_col[0]], ['==','=='], [exonic_definitions[0],exonic_definitions[1]])
			##get dispuptive
			filtering_annotated.filter(working_dir, "or", group + '.' + str(freq) +  ".exonic.xls", group + '.' + str(freq) +  ".disruptive.xls", make_list(case_controls_col[2], disruptive_definitions), make_list('==', disruptive_definitions), disruptive_definitions)
			##get all protein changing
			filtering_annotated.filter(working_dir, "or", group + '.' + str(freq) +  ".exonic.xls", group + '.' + str(freq) +  ".protein_changing.xls", make_list(case_controls_col[2], protein_changing_definitions), make_list('==', protein_changing_definitions), protein_changing_definitions)
			##get damaging - pp2_hdiv, pp2_hvar, cadd_phred, gerp - in all or in any
			#get non synonymous snps
			filtering_annotated.filter(working_dir, "or", group + '.' + str(freq) +  ".exonic.xls", group + '.' + str(freq) +  ".3.temp", [case_controls_col[2], case_controls_col[2]], ['==', '=='], nosyn_definitions)
			#get if all are positive
			# filtering_annotated.filter(working_dir, "and", group + '.' + str(freq) +  ".3.temp", group + '.' + str(freq) +  ".4.temp", case_controls_col[3], make_list('>=', case_controls_col[3]), damaging_definitions)
			filtering_annotated.filter(working_dir, "or", group + '.' + str(freq) +  ".3.temp", group + '.' + str(freq) +  ".4a.temp", case_controls_col[3][:2], make_list('>=', case_controls_col[3][:2]), pp2_score)
			filtering_annotated.filter(working_dir, "and", group + '.' + str(freq) +  ".4a.temp", group + '.' + str(freq) +  ".4b.temp", case_controls_col[3][2:], make_list('>=', case_controls_col[3][2:]), cadd_gerp_score)
			remove_rows_with_no_data(group + '.' + str(freq) +  ".4b.temp", group + '.' + str(freq) +  ".damaging_all.xls", case_controls_col[3])
			#get if any are positive
			filtering_annotated.filter(working_dir, "or", group + '.' + str(freq) +  ".3.temp", group + '.' + str(freq) +  ".5.temp", case_controls_col[3], make_list('>=', case_controls_col[3]), damaging_definitions)
			remove_rows_with_no_data(group + '.' + str(freq) +  ".5.temp", group + '.' + str(freq) +  ".damaging_any.xls", case_controls_col[3])
		##and exac
		for group in ['exac3_1115']:
			##only 'rare' (<=1%)
			filtering_annotated.filter(working_dir, "and", exac_txt, group + '.' + str(freq) +  ".1.temp", esp6500_rvis_col[1], make_list('<=', esp6500_rvis_col[1]), make_list(freq, esp6500_rvis_col[1]))
			##exonic_variants in refGene
			filtering_annotated.filter(working_dir, "or", group + '.' + str(freq) +  ".1.temp", group + '.' + str(freq) +  ".exonic.xls", [esp6500_rvis_col[0], esp6500_rvis_col[0]], ['==','=='], [exonic_definitions[0],exonic_definitions[1]])
			##get dispuptive
			filtering_annotated.filter(working_dir, "or", group + '.' + str(freq) +  ".exonic.xls", group + '.' + str(freq) +  ".disruptive.xls", make_list(esp6500_rvis_col[2], disruptive_definitions), make_list('==', disruptive_definitions), disruptive_definitions)
			##get all protein changing
			filtering_annotated.filter(working_dir, "or", group + '.' + str(freq) +  ".exonic.xls", group + '.' + str(freq) +  ".protein_changing.xls", make_list(esp6500_rvis_col[2], protein_changing_definitions), make_list('==', protein_changing_definitions), protein_changing_definitions)
			##get damaging - pp2_hdiv, pp2_hvar, cadd_phred, gerp - in all or in any
			#get non synonymous snps
			filtering_annotated.filter(working_dir, "or", group + '.' + str(freq) +  ".exonic.xls", group + '.' + str(freq) +  ".3.temp", [esp6500_rvis_col[2], esp6500_rvis_col[2]], ['==', '=='], nosyn_definitions)
			#get if all are positive
			# filtering_annotated.filter(working_dir, "and", group + '.' + str(freq) +  ".3.temp", group + '.' + str(freq) +  ".4.temp", esp6500_rvis_col[3], make_list('>=', esp6500_rvis_col[3]), damaging_definitions)
			filtering_annotated.filter(working_dir, "or", group + '.' + str(freq) +  ".3.temp", group + '.' + str(freq) +  ".4a.temp", esp6500_rvis_col[3][:2], make_list('>=', esp6500_rvis_col[3][:2]), pp2_score)
			filtering_annotated.filter(working_dir, "and", group + '.' + str(freq) +  ".4a.temp", group + '.' + str(freq) +  ".4b.temp", esp6500_rvis_col[3][2:], make_list('>=', esp6500_rvis_col[3][2:]), cadd_gerp_score)
			remove_rows_with_no_data(group + '.' + str(freq) +  ".4b.temp", group + '.' + str(freq) +  ".damaging_all.xls", esp6500_rvis_col[3])
			#get if any are positive
			filtering_annotated.filter(working_dir, "or", group + '.' + str(freq) +  ".3.temp", group + '.' + str(freq) +  ".5.temp", esp6500_rvis_col[3], make_list('>=', esp6500_rvis_col[3]), damaging_definitions)
			remove_rows_with_no_data(group + '.' + str(freq) +  ".5.temp", group + '.' + str(freq) +  ".damaging_any.xls", esp6500_rvis_col[3])

def filter_ann_txt_files_just_exac(exac_txt, esp6500_cols, freq_req):
	for freq in freq_req:
		for group in [exac_unfiltered_prefix]:
			# ##only 'rare' (<=1%)
			# filtering_annotated.filter(working_dir, "and", exac_txt, group + '.' + str(freq) +  ".1.temp", esp6500_rvis_col[1], make_list('<=', esp6500_rvis_col[1]), make_list(freq, esp6500_rvis_col[1]))
			# ##exonic_variants in refGene
			# filtering_annotated.filter(working_dir, "or", group + '.' + str(freq) +  ".1.temp", group + '.' + str(freq) +  ".exonic.xls", [esp6500_rvis_col[0], esp6500_rvis_col[0]], ['==','=='], [exonic_definitions[0],exonic_definitions[1]])

			##only 'rare' (<=1%)
			filtering_annotated.filter(working_dir, "and", exac_txt, group + '.' + str(freq) +  ".1.temp", esp6500_rvis_col[1], make_list('<=', esp6500_rvis_col[1]), make_list(freq, esp6500_rvis_col[1]))
			##exonic_variants in refGene
			filtering_annotated.filter(working_dir, "or", group + '.' + str(freq) +  ".1.temp", group + '.' + str(freq) +  ".exonic.xls", [esp6500_rvis_col[0], esp6500_rvis_col[0]], ['==','=='], [exonic_definitions[0],exonic_definitions[1]])
			##get dispuptive
			filtering_annotated.filter(working_dir, "or", group + '.' + str(freq) +  ".exonic.xls", group + '.' + str(freq) +  ".disruptive.xls", make_list(esp6500_rvis_col[2], disruptive_definitions), make_list('==', disruptive_definitions), disruptive_definitions)
			##get all protein changing
			filtering_annotated.filter(working_dir, "or", group + '.' + str(freq) +  ".exonic.xls", group + '.' + str(freq) +  ".protein_changing.xls", make_list(esp6500_rvis_col[2], protein_changing_definitions), make_list('==', protein_changing_definitions), protein_changing_definitions)
			##get damaging - pp2_hdiv, pp2_hvar, cadd_phred, gerp - in all or in any
			#get non synonymous snps
			filtering_annotated.filter(working_dir, "or", group + '.' + str(freq) +  ".exonic.xls", group + '.' + str(freq) +  ".3.temp", [esp6500_rvis_col[2], esp6500_rvis_col[2]], ['==', '=='], nosyn_definitions)
			#get if all are positive
			# filtering_annotated.filter(working_dir, "and", group + '.' + str(freq) +  ".3.temp", group + '.' + str(freq) +  ".4.temp", esp6500_rvis_col[3], make_list('>=', esp6500_rvis_col[3]), damaging_definitions)
			filtering_annotated.filter(working_dir, "or", group + '.' + str(freq) +  ".3.temp", group + '.' + str(freq) +  ".4a.temp", esp6500_rvis_col[3][:2], make_list('>=', esp6500_rvis_col[3][:2]), pp2_score)
			filtering_annotated.filter(working_dir, "and", group + '.' + str(freq) +  ".4a.temp", group + '.' + str(freq) +  ".4b.temp", esp6500_rvis_col[3][2:], make_list('>=', esp6500_rvis_col[3][2:]), cadd_gerp_score)
			remove_rows_with_no_data(group + '.' + str(freq) +  ".4b.temp", group + '.' + str(freq) +  ".damaging_all.xls", esp6500_rvis_col[3])
			#get if any are positive
			filtering_annotated.filter(working_dir, "or", group + '.' + str(freq) +  ".3.temp", group + '.' + str(freq) +  ".5.temp", esp6500_rvis_col[3], make_list('>=', esp6500_rvis_col[3]), damaging_definitions)
			remove_rows_with_no_data(group + '.' + str(freq) +  ".5.temp", group + '.' + str(freq) +  ".damaging_any.xls", esp6500_rvis_col[3])


def remove_dup_data_esp6500(infile, outfile):
	with open(outfile, "w") as out_fh:
		info_list = []
		with open(infile, "r") as in_fh:
			for line in in_fh:
				line = line.split(delim)
				info = line[esp6500_info_col -1]
				if info not in info_list:
					out_fh.write(delim.join(line))
					info_list.append(info)
				# else:
				# 	print info

def get_var_counts_per_gene(infile, genename_col, data_type):
	var_count_dict = {}
	line_count = 0
	if data_type == 'case_control' or data_type == 'exac3_1115':
		varfile = infile
	##esp6500 contain duplicate data when i split multialleleic split into different lines
	elif data_type == 'esp6500':
		varfile = 'esp6500.info_dups.temp'
		remove_dup_data_esp6500(infile, varfile)
	with open(varfile, "r") as in_fh:
		for line in in_fh:
			chrom = line[0]
			line_count += 1
			if line_count != 1:
				line = line.split(delim)
				genes = line[genename_col -1]
				for gene in genes.split(','):
					if data_type == 'case_control':
						count = int(line[0])
					elif data_type == 'esp6500':
						allelle_counts = line[esp6500_info_col -1].split(';')[8].split('=')[1].split(',')
						if chrom == 'X':
							count = sum(map(int, allelle_counts[:-2]))
							# print chrom, gene, allelle_counts, count
						else:
							count = sum(map(int, allelle_counts[:-1]))
					elif data_type == 'exac3_1115':
						allelle_counts = line[esp6500_info_col -1].split(';')[6:8]
						count = 0 
						for ac in allelle_counts:
							ac = ac.split('=')[1].split(',')
							# print ac
							for a in ac:
								a = int(a)
								count += a
						# print allelle_counts, count
					if gene in var_count_dict:
						var_count_dict[gene] += count
					else:
						var_count_dict[gene] = count
	# for g in var_count_dict:
	# 	print g, var_count_dict[g]
	##put count into a list
	for g in var_count_dict:
		var_count_dict[g] = [var_count_dict[g]]
	return var_count_dict

def fischers_test(cases, controls, esp6500):
		results = []
		print cases, controls, esp6500
		cases_oddsratio, cases_pvalue = stats.fisher_exact([esp6500,cases])
		controls_oddsratio, controls_pvalue = stats.fisher_exact([esp6500,controls])
		case_control_oddsratio, case_control_pvalue = stats.fisher_exact([controls,cases])
		results.extend([cases_pvalue, controls_pvalue, case_control_pvalue])
		if sum(cases) == 0 or sum(esp6500) == 0:
			cases_a = 'na'
		elif (float(cases[1]) / sum(cases)) > (float(esp6500[1]) / sum(esp6500)):
			cases_a = 'yes'
		else:
			cases_a = 'no'
		if sum(controls) == 0 or sum(esp6500) == 0:
			controls_a = 'na'
		elif (float(controls[1]) / sum(controls)) > (float(esp6500[1]) / sum(esp6500)):
			controls_a = 'yes'
		else:
			controls_a = 'no'
		if sum(cases) == 0 or sum(controls) == 0:
			case_control_a = 'na'
		elif (float(cases[1]) / sum(cases)) > (float(controls[1]) / sum(controls)):
			case_control_a = 'yes'
		else:
			case_control_a = 'no'
		results.extend([cases_a, controls_a, case_control_a])
		return results
# get_var_counts_per_gene('cases.exonic.xls', case_ctl_rvis_col[4], 'case_control')
# get_var_counts_per_gene('esp6500.disruptive.xls', esp6500_rvis_col[4], 'esp6500')


def rvis_lite(outfile):
	big_dict = get_var_counts_per_gene('cases.exonic.xls', case_ctl_rvis_col[4], 'case_control')
	cases_dis = get_var_counts_per_gene('cases.disruptive.xls', case_ctl_rvis_col[4], 'case_control')
	cases_dam_all = get_var_counts_per_gene('cases.damaging_all.xls', case_ctl_rvis_col[4], 'case_control')
	cases_dam_any= get_var_counts_per_gene('cases.damaging_any.xls', case_ctl_rvis_col[4], 'case_control')
	cases_pc = get_var_counts_per_gene('cases.protein_changing.xls', case_ctl_rvis_col[4], 'case_control')
	control_ex = get_var_counts_per_gene('controls.exonic.xls', case_ctl_rvis_col[4], 'case_control')
	control_dis = get_var_counts_per_gene('controls.disruptive.xls', case_ctl_rvis_col[4], 'case_control')
	control_dam_all = get_var_counts_per_gene('controls.damaging_all.xls', case_ctl_rvis_col[4], 'case_control')
	control_dam_any= get_var_counts_per_gene('controls.damaging_any.xls', case_ctl_rvis_col[4], 'case_control')
	control_pc = get_var_counts_per_gene('controls.protein_changing.xls', case_ctl_rvis_col[4], 'case_control')
	esp6500_ex = get_var_counts_per_gene('esp6500.exonic.xls', esp6500_rvis_col[4], 'esp6500')
	esp6500_dis = get_var_counts_per_gene('esp6500.disruptive.xls', esp6500_rvis_col[4], 'esp6500')
	esp6500_dam_all = get_var_counts_per_gene('esp6500.damaging_all.xls', esp6500_rvis_col[4], 'esp6500')
	esp6500_dam_any = get_var_counts_per_gene('esp6500.damaging_any.xls', esp6500_rvis_col[4], 'esp6500')
	esp6500_pc = get_var_counts_per_gene('esp6500.protein_changing.xls', esp6500_rvis_col[4], 'esp6500')
	##combine results from all dicts
	for gene in big_dict:
		for genedict in [cases_dis, cases_dam_all, cases_dam_any, cases_pc, control_ex, control_dis, control_dam_all, control_dam_any, control_pc, esp6500_ex, esp6500_dis, esp6500_dam_all, esp6500_dam_any, esp6500_pc]:
			if gene in genedict:
				big_dict[gene].extend(genedict[gene])
			else:
				big_dict[gene].extend([0])
	##do tests and add
	for gene in big_dict:
		cases_disruptive_vs = [big_dict[gene][0] - big_dict[gene][1], big_dict[gene][1]]
		cases_disruptive_plus_dall_vs = [(big_dict[gene][0] - (big_dict[gene][1] + big_dict[gene][2])), big_dict[gene][1] + big_dict[gene][2]]
		cases_disruptive_plus_dany_vs = [(big_dict[gene][0] - (big_dict[gene][1] + big_dict[gene][3])), big_dict[gene][1] + big_dict[gene][3]]
		cases_pc_vs = [(big_dict[gene][0] - big_dict[gene][4]), big_dict[gene][4]]

		controls_disruptive_vs = [(big_dict[gene][5] - big_dict[gene][6]), big_dict[gene][6]]
		controls_disruptive_plus_dall_vs = [(big_dict[gene][5] - (big_dict[gene][6] + big_dict[gene][7])), big_dict[gene][6] + big_dict[gene][7]]
		controls_disruptive_plus_dany_vs = [(big_dict[gene][5] - (big_dict[gene][6] + big_dict[gene][8])), big_dict[gene][6] + big_dict[gene][8]]
		controls_pc_vs = [(big_dict[gene][5] - big_dict[gene][9]), big_dict[gene][9]]

		esp6500_disruptive_vs = [(big_dict[gene][10] - big_dict[gene][11]), big_dict[gene][11]]
		esp6500_disruptive_plus_dall_vs = [(big_dict[gene][10] - (big_dict[gene][11] + big_dict[gene][12])), big_dict[gene][11] + big_dict[gene][12]]
		esp6500_disruptive_plus_dany_vs = [(big_dict[gene][10] - (big_dict[gene][11] + big_dict[gene][13])), big_dict[gene][11] + big_dict[gene][13]]
		esp6500_pc_vs = [(big_dict[gene][10] - big_dict[gene][14]), big_dict[gene][14]]
		for numbers in [cases_disruptive_plus_dall_vs, controls_disruptive_plus_dall_vs, esp6500_disruptive_plus_dall_vs]:
			for n in numbers:
				if n <0:
					print gene, numbers, n


		##fischers exact test
		disruptive_results = fischers_test(cases_disruptive_vs, controls_disruptive_vs, esp6500_disruptive_vs)
		disruptive_dall_results = fischers_test(cases_disruptive_plus_dall_vs, controls_disruptive_plus_dall_vs, esp6500_disruptive_plus_dall_vs) 
		disruptive_dany_results = fischers_test(cases_disruptive_plus_dany_vs, controls_disruptive_plus_dany_vs, esp6500_disruptive_plus_dany_vs) 
		pc_results = fischers_test(cases_pc_vs, controls_pc_vs, esp6500_pc_vs) 

		##add to dict
		big_dict[gene].extend(disruptive_results + disruptive_dall_results + disruptive_dany_results + pc_results)

	# for gene in big_dict:
	# 	print gene, big_dict[gene]

	##covert all values to str and then write to file
	with open(outfile, "w") as out_fh:
		head = ['gene', 'cases_exonic', 'cases_disruptive', 'cases_damaging_all', 'cases_damaging_any', 'cases_protein_changing', 'controls_exonic', 'controls_disruptive', 'controls_damaging_all', 'controls_damaging_any', 'controls_protein_changing', 'esp6500_exonic', 'esp6500_disruptive', 'esp6500_damaging_all', 'esp6500_damaging_any', 'esp6500_protein_changing', 'case_esp6500_disruptive_pvalue', 'control_esp6500_disruptive_pvalue', 'case_control_disruptive_pvalue', 'case_esp6500_disruptive_case_dir', 'control_esp6500_disruptive_case_dir', 'case_control_disruptive_case_dir', 'case_esp6500_dall_pvalue', 'control_esp6500_dall_pvalue', 'case_control_dall_pvalue', 'case_esp6500_dall_case_dir', 'control_esp6500_dall_case_dir', 'case_control_dall_case_dir', 'case_esp6500_dany_pvalue', 'control_esp6500_dany_pvalue', 'case_control_dany_pvalue', 'case_esp6500_dany_case_dir', 'control_esp6500_dany_case_dir', 'case_control_dany_case_dir', 'case_esp6500_pc_pvalue', 'control_esp6500_pc_pvalue', 'case_control_pc_pvalue', 'case_esp6500_pc_case_dir', 'control_esp6500_pc_case_dir', 'case_control_pc_case_dir', '\n']
		out_fh.write(delim.join(head))
		big_dict_str = {}
		for gene in big_dict:
			number_list = []
			for i in big_dict[gene]:
				i_str = str(i)
				number_list.append(i_str)
			big_dict_str[gene] = number_list
		for gene in big_dict_str:
			out_fh.write(delim.join([gene] + big_dict_str[gene] + ['\n']))



def rvis_lite_exac(freqs, outfile_suffix):
	for freq in freqs:
		outfile = 'css_rnaseq.' + str(freq) + outfile_suffix
		big_dict = get_var_counts_per_gene('cases.' + str(freq) + '.exonic.xls', case_ctl_rvis_col[4], 'case_control')
		cases_dis = get_var_counts_per_gene('cases.' + str(freq) + '.disruptive.xls', case_ctl_rvis_col[4], 'case_control')
		cases_dam_all = get_var_counts_per_gene('cases.' + str(freq) + '.damaging_all.xls', case_ctl_rvis_col[4], 'case_control')
		cases_dam_any= get_var_counts_per_gene('cases.' + str(freq) + '.damaging_any.xls', case_ctl_rvis_col[4], 'case_control')
		cases_pc = get_var_counts_per_gene('cases.' + str(freq) + '.protein_changing.xls', case_ctl_rvis_col[4], 'case_control')
		control_ex = get_var_counts_per_gene('controls.' + str(freq) + '.exonic.xls', case_ctl_rvis_col[4], 'case_control')
		control_dis = get_var_counts_per_gene('controls.' + str(freq) + '.disruptive.xls', case_ctl_rvis_col[4], 'case_control')
		control_dam_all = get_var_counts_per_gene('controls.' + str(freq) + '.damaging_all.xls', case_ctl_rvis_col[4], 'case_control')
		control_dam_any= get_var_counts_per_gene('controls.' + str(freq) + '.damaging_any.xls', case_ctl_rvis_col[4], 'case_control')
		control_pc = get_var_counts_per_gene('controls.' + str(freq) + '.protein_changing.xls', case_ctl_rvis_col[4], 'case_control')
		exac3_1115_ex = get_var_counts_per_gene('exac3_1115.' + str(freq) + '.exonic.xls', esp6500_rvis_col[4], 'exac3_1115')
		exac3_1115_dis = get_var_counts_per_gene('exac3_1115.' + str(freq) + '.disruptive.xls', esp6500_rvis_col[4], 'exac3_1115')
		exac3_1115_dam_all = get_var_counts_per_gene('exac3_1115.' + str(freq) + '.damaging_all.xls', esp6500_rvis_col[4], 'exac3_1115')
		exac3_1115_dam_any = get_var_counts_per_gene('exac3_1115.' + str(freq) + '.damaging_any.xls', esp6500_rvis_col[4], 'exac3_1115')
		exac3_1115_pc = get_var_counts_per_gene('exac3_1115.' + str(freq) + '.protein_changing.xls', esp6500_rvis_col[4], 'exac3_1115')
		##combine results from all dicts
		for gene in big_dict:
			for genedict in [cases_dis, cases_dam_all, cases_dam_any, cases_pc, control_ex, control_dis, control_dam_all, control_dam_any, control_pc, exac3_1115_ex, exac3_1115_dis, exac3_1115_dam_all, exac3_1115_dam_any, exac3_1115_pc]:
				if gene in genedict:
					big_dict[gene].extend(genedict[gene])
				else:
					big_dict[gene].extend([0])
		##do tests and add
		for gene in big_dict:
			cases_disruptive_vs = [big_dict[gene][0] - big_dict[gene][1], big_dict[gene][1]]
			cases_disruptive_plus_dall_vs = [(big_dict[gene][0] - (big_dict[gene][1] + big_dict[gene][2])), big_dict[gene][1] + big_dict[gene][2]]
			cases_disruptive_plus_dany_vs = [(big_dict[gene][0] - (big_dict[gene][1] + big_dict[gene][3])), big_dict[gene][1] + big_dict[gene][3]]
			cases_pc_vs = [(big_dict[gene][0] - big_dict[gene][4]), big_dict[gene][4]]

			controls_disruptive_vs = [(big_dict[gene][5] - big_dict[gene][6]), big_dict[gene][6]]
			controls_disruptive_plus_dall_vs = [(big_dict[gene][5] - (big_dict[gene][6] + big_dict[gene][7])), big_dict[gene][6] + big_dict[gene][7]]
			controls_disruptive_plus_dany_vs = [(big_dict[gene][5] - (big_dict[gene][6] + big_dict[gene][8])), big_dict[gene][6] + big_dict[gene][8]]
			controls_pc_vs = [(big_dict[gene][5] - big_dict[gene][9]), big_dict[gene][9]]

			exac3_1115_disruptive_vs = [(big_dict[gene][10] - big_dict[gene][11]), big_dict[gene][11]]
			exac3_1115_disruptive_plus_dall_vs = [(big_dict[gene][10] - (big_dict[gene][11] + big_dict[gene][12])), big_dict[gene][11] + big_dict[gene][12]]
			exac3_1115_disruptive_plus_dany_vs = [(big_dict[gene][10] - (big_dict[gene][11] + big_dict[gene][13])), big_dict[gene][11] + big_dict[gene][13]]
			exac3_1115_pc_vs = [(big_dict[gene][10] - big_dict[gene][14]), big_dict[gene][14]]
			for numbers in [cases_disruptive_plus_dall_vs, controls_disruptive_plus_dall_vs, exac3_1115_disruptive_plus_dall_vs]:
				for n in numbers:
					if n <0:
						print gene, numbers, n


			##fischers exact test
			disruptive_results = fischers_test(cases_disruptive_vs, controls_disruptive_vs, exac3_1115_disruptive_vs)
			disruptive_dall_results = fischers_test(cases_disruptive_plus_dall_vs, controls_disruptive_plus_dall_vs, exac3_1115_disruptive_plus_dall_vs) 
			disruptive_dany_results = fischers_test(cases_disruptive_plus_dany_vs, controls_disruptive_plus_dany_vs, exac3_1115_disruptive_plus_dany_vs) 
			pc_results = fischers_test(cases_pc_vs, controls_pc_vs, exac3_1115_pc_vs) 

			##add to dict
			big_dict[gene].extend(disruptive_results + disruptive_dall_results + disruptive_dany_results + pc_results)

		# for gene in big_dict:
		# 	print gene, big_dict[gene]

		##covert all values to str and then write to file
		with open(outfile, "w") as out_fh:
			head = ['gene', 'cases_exonic', 'cases_disruptive', 'cases_damaging_all', 'cases_damaging_any', 'cases_protein_changing', 'controls_exonic', 'controls_disruptive', 'controls_damaging_all', 'controls_damaging_any', 'controls_protein_changing', 'exac3_1115_exonic', 'exac3_1115_disruptive', 'exac3_1115_damaging_all', 'exac3_1115_damaging_any', 'exac3_1115_protein_changing', 'case_exac3_1115_disruptive_pvalue', 'control_exac3_1115_disruptive_pvalue', 'case_control_disruptive_pvalue', 'case_exac3_1115_disruptive_case_dir', 'control_exac3_1115_disruptive_case_dir', 'case_control_disruptive_case_dir', 'case_exac3_1115_dall_pvalue', 'control_exac3_1115_dall_pvalue', 'case_control_dall_pvalue', 'case_exac3_1115_dall_case_dir', 'control_exac3_1115_dall_case_dir', 'case_control_dall_case_dir', 'case_exac3_1115_dany_pvalue', 'control_exac3_1115_dany_pvalue', 'case_control_dany_pvalue', 'case_exac3_1115_dany_case_dir', 'control_exac3_1115_dany_case_dir', 'case_control_dany_case_dir', 'case_exac3_1115_pc_pvalue', 'control_exac3_1115_pc_pvalue', 'case_control_pc_pvalue', 'case_exac3_1115_pc_case_dir', 'control_exac3_1115_pc_case_dir', 'case_control_pc_case_dir', '\n']
			out_fh.write(delim.join(head))
			big_dict_str = {}
			for gene in big_dict:
				number_list = []
				for i in big_dict[gene]:
					i_str = str(i)
					number_list.append(i_str)
				big_dict_str[gene] = number_list
			for gene in big_dict_str:
				out_fh.write(delim.join([gene] + big_dict_str[gene] + ['\n']))


##prepare files

##remove not passed from vcf
# remove_not_passed_from_vcf(exac_vcf, exac_passed_vcf)

##annotate all exac variants
# annotate_vcf_file(exac_passed_vcf, exac_avinput, exac_prefix)
# annotate_vcf_file(exac_vcf, exac_avinput, exac_unfiltered_prefix)


##make tests i.e. rvis_lite and case control

##filter variants
# filter_ann_txt_files(exac_annotated, cases_passed_annotated, controls_passed_annotated, esp6500_rvis_col, case_ctl_rvis_col, frequencies)
filter_ann_txt_files_just_exac(exac_unfiltered_annotated, esp6500_rvis_col, frequencies)

##exac
##rvis_lite - done on local machine
# rvis_lite_exac(frequencies, rvis_lite_exac_out_file_suffix)






