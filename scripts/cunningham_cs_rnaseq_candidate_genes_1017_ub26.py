#!/usr/bin/env python
import filtering_annotated
import os
import subprocess
import glob
import re

##set input variables and parameters
delim = '\t'
##working directory
working_dir = '/data/atimms/cs_rnaseq_1116'
os.chdir(working_dir)


def convert_list_to_dict_with_zero(gene_list):
	##convert list to dict to track var_count
	gene_dict = {}
	for gene in gene_list:
		gene_dict[gene] = 0
	return gene_dict

def filter_var_by_genename(genelist, input_file,output_file, gene_col):
	##convert list to dict to track var_count
	genedict = convert_list_to_dict_with_zero(genelist)
# 	print genedict
# 	print len(genedict), len(genelist)
	with open(output_file, "w") as outfh, open(input_file, "U") as infh:
		line_count = 0
		for line in infh:
			line_count += 1
			if line_count == 1:
				outfh.write(line)
			else:
				line = line.strip('\n').strip('\r').rstrip().split(delim)
				genes = re.sub(r'\(.+?\)', '', line[gene_col])
				genes = genes.split(',')
				for gene in genes:
					if gene in genedict:
						outfh.write(delim.join(line) + '\n')
						genedict[gene] += 1
						break
		#print genedict
		print 'of the variants checked %i are in the genes specified'% (sum(genedict.values()))
		print '\n'
# 		for g in genedict:
# 			print 'in file %s we have %s variants in gene %s'%(output_file, genedict[g], g)


def combine_ann_txt(sample_list, prefix, suffix, outfile):
	with open(outfile, "w") as final_file:
		sample_count = 0
		for sample in sample_list:
			file = prefix + sample + suffix
			sample_count += 1
			with open(file, "r") as open_file:
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


##simple list maker -- give string and dictionary to get length from
def make_list(string, dict):
	l = []
	for i in range(len(dict)):
		l.append(string)
	return l

def filter_ann_txt_files(infile, cols_to_filter, freq_req, outfile_prefix, protein_changing_definitions, exonic_definitions, genename):
	##filter rnaseq data
	for freq in freq_req:
		##only 'rare' (<=1%)
		filtering_annotated.filter(working_dir, "and", infile, "1.temp", [cols_to_filter[0]], ['<='], [freq])
		##q>=30 and coverage >=5
		filtering_annotated.filter(working_dir, "and", "1.temp", "2.temp", cols_to_filter[1], ['>=', '>='], [30,5])
		##exonic_variants in refGene
		filtering_annotated.filter(working_dir, "or", "2.temp", outfile_prefix + '.' + str(freq) +  ".exonic_temp.xls", [cols_to_filter[2], cols_to_filter[2]], ['==','=='], [exonic_definitions[0],exonic_definitions[1]])
		##get all protein changing
		filtering_annotated.filter(working_dir, "or", outfile_prefix + '.' + str(freq) +  ".exonic_temp.xls", outfile_prefix + '.' + str(freq) +  ".protein_changing." + genename + ".xls", make_list(cols_to_filter[3], protein_changing_definitions), make_list('==', protein_changing_definitions), protein_changing_definitions)




def make_maf_dict_from_pc_vars(var_file, maf_col):
	maf_dict = {}
	with open(var_file, "r") as infh:
		line_count = 0
		for line in infh:
			line_count += 1
			if line_count > 1:
				line = line.split(delim)
				pro = line[0]
				maf = line[maf_col]
				# print pro, maf
				if pro in maf_dict:
					maf_dict[pro].append(maf)
				else:
					maf_dict[pro] = [maf]
	return maf_dict

def combine_all_data_make_final_file(pro_dict, file_suffix, outfile):
	with open(outfile, "w") as outfh:
		file_count = 0
		for pro in pro_dict:
			input_file = pro + file_suffix
			file_count += 1
			with open(input_file, "r") as infh:
				line_count = 0
				for line in infh:
					line = line.split(delim)
					line_count += 1
					if line_count == 1:
						if file_count == 1:
							outfh.write(delim.join(['proband', 'maf'] + line))
					else:
						outfh.write(delim.join([pro,'_'.join(pro_dict[pro])] + line))


def filter_ann_txt_files_just_maf(infile, cols_to_filter, freq_req, outfile_prefix, protein_changing_definitions, exonic_definitions, genename):
	##filter rnaseq data
	for freq in freq_req:
		##only 'rare' (<=1%)
		filtering_annotated.filter(working_dir, "and", infile, "1.temp", [cols_to_filter[0]], ['<='], [freq])
		##q>=30 and coverage >=5
		filtering_annotated.filter(working_dir, "and", "1.temp", outfile_prefix + '.' + str(freq) + '.' + genename + ".xls", cols_to_filter[1], ['>=', '>='], [30,5])


def summarize_final_file(in_file, out_file):
	gene_dict = {}
	with open(in_file, "r") as infh:
		line_count = 0
		for line in infh:
			line_count += 1
			if line_count > 1:
				line = line.rstrip().split(delim)
				gene = line[17]
				proband = line[0]
				var = '_'.join(line[2:7])
				if gene in gene_dict:
					if proband not in gene_dict[gene][0]:
						gene_dict[gene][0].append(proband)
					gene_dict[gene][1].append(var)
					if var not in gene_dict[gene][2]:
						gene_dict[gene][2].append(var)
				else:
					 gene_dict[gene] = [[proband], [var], [var]]

	with open(out_file, "w") as outfh:
		outfh.write(delim.join(['gene', 'proband#', 'probands', 'total_variant#', 'unique_variant#', '\n']))
		for g in gene_dict:
			outfh.write(delim.join([g, str(len(gene_dict[g][0])), '_'.join(gene_dict[g][0]), str(len(gene_dict[g][1])), str(len(gene_dict[g][2])), '\n']))






# genes_of_interest = ['PIEZO1', 'AXL', 'FLNA', 'FLNB', 'FLNC', 'PDGFRA']
genes_of_interest = ['FLNA'] 
cases_list = ['cohort1.cases.95445', 'cohort1.cases.95446', 'cohort1.cases.95447', 'cohort1.cases.95448', 'cohort1.cases.95450', 'cohort1.cases.95452', 'cohort1.cases.95453', 'cohort1.cases.95454', 'cohort1.cases.95456', 'cohort1.cases.95457', 'cohort1.cases.95458', 'cohort1.cases.95459', 'cohort1.cases.95460', 'cohort1.cases.95462', 'cohort1.cases.95463', 'cohort1.cases.95464', 'cohort1.cases.95466', 'cohort1.cases.95467', 'cohort1.cases.95468', 'cohort1.cases.95469', 'cohort1.cases.95470', 'cohort1.cases.95471', 'cohort1.cases.95472', 'cohort1.cases.95474', 'cohort1.cases.95475', 'cohort1.cases.95476', 'cohort1.cases.95477', 'cohort1.cases.95479', 'cohort1.cases.95481', 'cohort1.cases.95482', 'cohort1.cases.95484', 'cohort1.cases.95485', 'cohort1.cases.95486', 'cohort1.cases.95487', 'cohort1.cases.95488', 'cohort1.cases.95490', 'cohort1.cases.95491', 'cohort1.cases.95492', 'cohort1.cases.95493', 'cohort1.cases.95494', 'cohort1.cases.95495', 'cohort1.cases.95496', 'cohort1.cases.95498', 'cohort1.cases.95500', 'cohort1.cases.95501', 'cohort1.cases.95502', 'cohort1.cases.95505', 'cohort1.cases.95506', 'cohort1.cases.95507', 'cohort1.cases.95509', 'cohort1.cases.95510', 'cohort1.cases.95511', 'cohort1.cases.95512', 'cohort1.cases.95514', 'cohort1.cases.95515', 'cohort1.cases.95516', 'cohort1.cases.95518', 'cohort1.cases.95519', 'cohort1.cases.95520', 'cohort1.cases.95521', 'cohort1.cases.95522', 'cohort1.cases.95524', 'cohort1.cases.95525', 'cohort1.cases.95526', 'cohort1.cases.95527', 'cohort1.cases.95528', 'cohort1.cases.95529', 'cohort1.cases.95530', 'cohort1.cases.95531', 'cohort1.cases.95534', 'cohort1.cases.95535', 'cohort1.cases.95536', 'cohort1.cases.95537', 'cohort1.cases.95540', 'cohort1.cases.95542', 'cohort1.cases.95543', 'cohort1.cases.95544', 'cohort1.cases.95545', 'cohort1.cases.95546', 'cohort1.cases.95547', 'cohort1.cases.95549', 'cohort1.cases.95550', 'cohort1.cases.95551', 'cohort1.cases.95552', 'cohort1.cases.95555', 'cohort1.cases.95556', 'cohort1.cases.95557', 'cohort1.cases.95558', 'cohort1.cases.95559', 'cohort1.cases.95560', 'cohort1.cases.95562', 'cohort1.cases.95563', 'cohort1.cases.95564', 'cohort1.cases.95565', 'cohort1.cases.95566', 'cohort1.cases.95567', 'cohort1.cases.95569', 'cohort1.cases.95570', 'cohort1.cases.95571', 'cohort1.cases.95572', 'cohort1.cases.95573', 'cohort1.cases.95575', 'cohort1.cases.95576', 'cohort1.cases.95578', 'cohort1.cases.95579', 'cohort1.cases.95580', 'cohort1.cases.95581', 'cohort1.cases.95583', 'cohort1.cases.95584', 'cohort1.cases.95585', 'cohort1.cases.95586', 'cohort1.cases.95587', 'cohort1.cases.95588', 'cohort1.cases.95589', 'cohort1.cases.95591', 'cohort1.cases.95592', 'cohort1.cases.95593', 'cohort1.cases.95594', 'cohort1.cases.95595', 'cohort1.cases.95596', 'cohort1.cases.95597', 'cohort1.cases.95598', 'cohort1.cases.95599', 'cohort1.cases.95600', 'cohort1.cases.95601', 'cohort1.cases.95602', 'cohort1.cases.95603', 'cohort1.cases.95604', 'cohort1.cases.95605', 'cohort1.cases.95606', 'cohort1.cases.95609', 'cohort1.cases.95610', 'cohort1.cases.95611', 'cohort1.cases.95612', 'cohort1.cases.95614', 'cohort1.cases.95615', 'cohort1.cases.95616', 'cohort1.cases.95617', 'cohort1.cases.95619', 'cohort1.cases.95620', 'cohort1.cases.95622', 'cohort1.cases.95623', 'cohort1.cases.95624', 'cohort1.cases.95626', 'cohort1.cases.95627', 'cohort1.cases.95628', 'cohort1.cases.95629', 'cohort1.cases.95630', 'cohort1.cases.95631', 'cohort1.cases.95633', 'cohort1.cases.95634', 'cohort1.cases.95635', 'cohort1.cases.95636', 'cohort1.cases.95637', 'cohort1.cases.95638', 'cohort1.cases.95639', 'cohort1.cases.95640', 'cohort1.cases.95642', 'cohort1.cases.95644', 'cohort1.cases.95645', 'cohort1.cases.95646', 'cohort1.cases.95649', 'cohort1.cases.95650', 'cohort1.cases.95651', 'cohort1.cases.95652', 'cohort1.cases.95653', 'cohort1.cases.95654', 'cohort1.cases.95655', 'cohort1.cases.95656', 'cohort1.cases.95657', 'cohort1.cases.95658', 'cohort1.cases.95661', 'cohort1.cases.95662', 'cohort1.cases.95663', 'cohort1.cases.95665', 'cohort1.cases.95668', 'cohort1.cases.95669', 'cohort1.cases.95670', 'cohort1.cases.95671', 'cohort1.cases.95673', 'cohort1.cases.95675', 'cohort1.cases.95676', 'cohort1.cases.95677', 'cohort1.cases.95678', 'cohort1.cases.95679', 'cohort1.cases.95680', 'cohort1.cases.95681', 'cohort1.cases.95682', 'cohort1.cases.95683', 'cohort1.cases.95684', 'cohort1.cases.95685', 'cohort1.cases.95686', 'cohort1.cases.95688', 'cohort1.cases.95689', 'cohort1.cases.95690', 'cohort1.cases.95691', 'cohort1.cases.95693', 'cohort1.cases.95694', 'cohort1.cases.95696', 'cohort1.cases.95697', 'cohort1.cases.95698', 'cohort1.cases.95699', 'cohort1.cases.95700', 'cohort1.cases.95702', 'cohort1.cases.95703', 'cohort2.cases.163499', 'cohort2.cases.163502', 'cohort2.cases.163503', 'cohort2.cases.163504', 'cohort2.cases.163505', 'cohort2.cases.163506', 'cohort2.cases.163507', 'cohort2.cases.163508', 'cohort2.cases.163509', 'cohort2.cases.163510', 'cohort2.cases.163511', 'cohort2.cases.163513', 'cohort2.cases.163516', 'cohort2.cases.163517', 'cohort2.cases.163518', 'cohort2.cases.163519', 'cohort2.cases.163520', 'cohort2.cases.163521', 'cohort2.cases.163522', 'cohort2.cases.163523', 'cohort2.cases.163524', 'cohort2.cases.163525', 'cohort2.cases.163526', 'cohort2.cases.163527', 'cohort2.cases.163528', 'cohort2.cases.163529', 'cohort2.cases.163530', 'cohort2.cases.163531', 'cohort2.cases.163532', 'cohort2.cases.163533', 'cohort2.cases.163534', 'cohort2.cases.163535', 'cohort2.cases.163537', 'cohort2.cases.163538', 'cohort2.cases.163539', 'cohort2.cases.163540', 'cohort2.cases.163541', 'cohort2.cases.163543', 'cohort2.cases.163544', 'cohort2.cases.163545', 'cohort2.cases.163546', 'cohort2.cases.163547', 'cohort2.cases.163548', 'cohort2.cases.163549', 'cohort2.cases.163553', 'cohort2.cases.163555', 'cohort2.cases.163556', 'cohort2.cases.163557', 'cohort2.cases.163558', 'cohort2.cases.163559', 'cohort2.cases.163560', 'cohort2.cases.163561', 'cohort2.cases.163562', 'cohort2.cases.163563', 'cohort2.cases.163564', 'cohort2.cases.163565', 'cohort2.cases.163566', 'cohort2.cases.163567', 'cohort2.cases.163568', 'cohort2.cases.163569', 'cohort2.cases.163571', 'cohort2.cases.163572', 'cohort2.cases.163573', 'cohort2.cases.163575', 'cohort2.cases.163576', 'cohort2.cases.163577', 'cohort2.cases.163578', 'cohort2.cases.163579', 'cohort2.cases.163580', 'cohort2.cases.163581', 'cohort2.cases.163582', 'cohort2.cases.163584', 'cohort2.cases.163585', 'cohort2.cases.163586', 'cohort2.cases.163587', 'cohort2.cases.163588', 'cohort2.cases.163590', 'cohort2.cases.163591', 'cohort2.cases.163592', 'cohort2.cases.163593', 'cohort2.cases.163594', 'cohort2.cases.163595', 'cohort2.cases.163596', 'cohort2.cases.163597', 'cohort2.cases.163598', 'cohort2.cases.163601', 'cohort2.cases.163602', 'cohort2.cases.163603', 'cohort2.cases.163604', 'cohort2.cases.163605', 'cohort2.cases.163606', 'cohort2.cases.163607', 'cohort2.cases.163609', 'cohort2.cases.163610', 'cohort2.cases.163611', 'cohort2.cases.163614', 'cohort2.cases.163615', 'cohort2.cases.163616', 'cohort2.cases.163619', 'cohort2.cases.163620', 'cohort2.cases.163621', 'cohort2.cases.163622', 'cohort2.cases.163623', 'cohort2.cases.163624', 'cohort2.cases.163625', 'cohort2.cases.163626', 'cohort2.cases.163627', 'cohort2.cases.163628', 'cohort2.cases.163629', 'cohort2.cases.163630', 'cohort2.cases.163631', 'cohort2.cases.163632', 'cohort2.cases.163633', 'cohort2.cases.163634', 'cohort2.cases.163635', 'cohort2.cases.163636', 'cohort2.cases.163637', 'cohort2.cases.163639', 'cohort2.cases.163640', 'cohort2.cases.163641', 'cohort2.cases.163642', 'cohort2.cases.163643', 'cohort2.cases.163644', 'cohort2.cases.163645', 'cohort2.cases.163646', 'cohort2.cases.163647', 'cohort2.cases.163649', 'cohort2.cases.163651', 'cohort2.cases.163652', 'cohort2.cases.163653', 'cohort2.cases.163654', 'cohort2.cases.163655', 'cohort2.cases.163656', 'cohort2.cases.163657', 'cohort2.cases.163658', 'cohort2.cases.163659', 'cohort2.cases.163660', 'cohort2.cases.163661', 'cohort2.cases.163662', 'cohort2.cases.163663', 'cohort2.cases.163664', 'cohort2.cases.163665', 'cohort2.cases.163666', 'cohort2.cases.163667', 'cohort2.cases.163668', 'cohort2.cases.163669', 'cohort2.cases.163672', 'cohort2.cases.163674', 'cohort2.cases.163678', 'cohort2.cases.163679', 'cohort2.cases.163680', 'cohort2.cases.163681', 'cohort2.cases.163682', 'cohort2.cases.163683', 'cohort2.cases.163684', 'cohort2.cases.163685', 'cohort2.cases.163686', 'cohort2.cases.163687', 'cohort2.cases.163688', 'cohort2.cases.163690', 'cohort2.cases.163691', 'cohort2.cases.163692', 'cohort2.cases.163694', 'cohort2.cases.163695', 'cohort2.cases.163696', 'cohort2.cases.163697', 'cohort2.cases.163698', 'cohort2.cases.163699', 'cohort2.cases.163700', 'cohort2.cases.163701', 'cohort2.cases.163702', 'cohort2.cases.163703', 'cohort2.cases.163704', 'cohort2.cases.163706', 'cohort2.cases.163707', 'cohort2.cases.163709', 'cohort2.cases.163710', 'cohort2.cases.163711', 'cohort2.cases.163712', 'cohort2.cases.163713', 'cohort2.cases.163714', 'cohort2.cases.163717', 'cohort2.cases.163718', 'cohort2.cases.163719', 'cohort2.cases.163720', 'cohort2.cases.163721']
control_list = ['cohort1.ctls.95449', 'cohort1.ctls.95451', 'cohort1.ctls.95455', 'cohort1.ctls.95461', 'cohort1.ctls.95465', 'cohort1.ctls.95473', 'cohort1.ctls.95478', 'cohort1.ctls.95480', 'cohort1.ctls.95483', 'cohort1.ctls.95489', 'cohort1.ctls.95497', 'cohort1.ctls.95499', 'cohort1.ctls.95504', 'cohort1.ctls.95508', 'cohort1.ctls.95513', 'cohort1.ctls.95517', 'cohort1.ctls.95523', 'cohort1.ctls.95532', 'cohort1.ctls.95533', 'cohort1.ctls.95538', 'cohort1.ctls.95539', 'cohort1.ctls.95541', 'cohort1.ctls.95548', 'cohort1.ctls.95553', 'cohort1.ctls.95554', 'cohort1.ctls.95561', 'cohort1.ctls.95568', 'cohort1.ctls.95574', 'cohort1.ctls.95577', 'cohort1.ctls.95582', 'cohort1.ctls.95590', 'cohort1.ctls.95607', 'cohort1.ctls.95608', 'cohort1.ctls.95613', 'cohort1.ctls.95618', 'cohort1.ctls.95621', 'cohort1.ctls.95625', 'cohort1.ctls.95647', 'cohort1.ctls.95648', 'cohort1.ctls.95659', 'cohort1.ctls.95660', 'cohort1.ctls.95666', 'cohort1.ctls.95667', 'cohort1.ctls.95672', 'cohort1.ctls.95674', 'cohort1.ctls.95692', 'cohort1.ctls.95695', 'cohort1.ctls.95704', 'cohort2.ctls.163500', 'cohort2.ctls.163501', 'cohort2.ctls.163512', 'cohort2.ctls.163514', 'cohort2.ctls.163515', 'cohort2.ctls.163536', 'cohort2.ctls.163542', 'cohort2.ctls.163550', 'cohort2.ctls.163551', 'cohort2.ctls.163552', 'cohort2.ctls.163554', 'cohort2.ctls.163570', 'cohort2.ctls.163574', 'cohort2.ctls.163583', 'cohort2.ctls.163589', 'cohort2.ctls.163599', 'cohort2.ctls.163600', 'cohort2.ctls.163608', 'cohort2.ctls.163612', 'cohort2.ctls.163613', 'cohort2.ctls.163617', 'cohort2.ctls.163618', 'cohort2.ctls.163638', 'cohort2.ctls.163648', 'cohort2.ctls.163650', 'cohort2.ctls.163670', 'cohort2.ctls.163671', 'cohort2.ctls.163673', 'cohort2.ctls.163675', 'cohort2.ctls.163676', 'cohort2.ctls.163677', 'cohort2.ctls.163689', 'cohort2.ctls.163693', 'cohort2.ctls.163705', 'cohort2.ctls.163708', 'cohort2.ctls.163715', 'cohort2.ctls.163716']

mafs_required = [0.04]
##for filtering variants: refgene function, frequency dbs, refgene exonic function, pp2 x2/cadd/gerp, refgene geneame, quality and coverage
##using exac all
filtering_cols = [64, [8,9], 16, 19]
filtering_cols_ind = [63, [7,8], 15, 18]
pc_defs =['','.', 'stopgain', 'stoploss', 'nonsynonymous SNV','frameshift insertion', 'frameshift deletion', 'nonsynonymous_SNV','frameshift_insertion', 'frameshift_deletion']
exonic_defs = ['exonic', 'splicing']

##get all protein changing variants at a maf for all cases to be used later
# for case in cases_list:
# 	filter_ann_txt_files(case + '.annotated.txt', filtering_cols_ind, mafs_required, case, pc_defs, exonic_defs, 'all_genes')

##get all variants in set of genes
'''
for gene in genes_of_interest:
	##get all variants in gene
	for case in cases_list:
		filter_var_by_genename([gene], case + '.annotated.txt', case + '.' + gene + '.temp.xls', 15)
	##combine in one file
	combine_ann_txt(cases_list, '', '.' + gene + '.temp.xls', 'cases.' +  gene + '.all_vars.xls')
	filter_ann_txt_files('cases.' +  gene + '.all_vars.xls', filtering_cols, mafs_required, 'cases', pc_defs, exonic_defs, gene)
	for maf in mafs_required:
		##make dict for each proband with a gene var
		pc_vars = 'cases.' + str(maf) + ".protein_changing." + gene + ".xls"
		var_maf_dict = make_maf_dict_from_pc_vars(pc_vars, filtering_cols[0])
		# for p in var_maf_dict:
		# 	print p, var_maf_dict[p]
		##get other variants and make final file
		combine_all_data_make_final_file(var_maf_dict, '.' +  str(maf) + '.protein_changing.all_genes.xls', 'cases.' + str(maf) + '.protein_changing.' + gene + '_associated.xls')
'''

##get allvariants at a maf of 0.04 for all cases (FLNA associated)
##not used
'''
cases_list = ['cohort1.cases.95479', 'cohort2.cases.163641', 'cohort1.cases.95661', 'cohort2.cases.163603', 'cohort2.cases.163714', 'cohort2.cases.163606', 'cohort1.cases.95597', 'cohort1.cases.95545', 'cohort2.cases.163560', 'cohort1.cases.95624']
for case in cases_list:
	filter_ann_txt_files_just_maf(case + '.annotated.txt', filtering_cols_ind, mafs_required, case, pc_defs, exonic_defs, 'all_genes')
for freq in mafs_required:
	combine_ann_txt(cases_list, '',  '.' + str(freq) + ".all_genes.xls", 'FLNA_cases.all_genes.' +  str(freq) + '.all_vars.xls')
'''
##second way
genes_of_interest = ['FLNA'] 
cases_list = ['cohort1.cases.95479', 'cohort2.cases.163641', 'cohort1.cases.95661', 'cohort2.cases.163603', 'cohort2.cases.163714', 'cohort2.cases.163606', 'cohort1.cases.95597', 'cohort1.cases.95545', 'cohort2.cases.163560', 'cohort1.cases.95624']
control_list = ['cohort1.ctls.95449', 'cohort1.ctls.95451', 'cohort1.ctls.95455', 'cohort1.ctls.95461', 'cohort1.ctls.95465', 'cohort1.ctls.95473', 'cohort1.ctls.95478', 'cohort1.ctls.95480', 'cohort1.ctls.95483', 'cohort1.ctls.95489', 'cohort1.ctls.95497', 'cohort1.ctls.95499', 'cohort1.ctls.95504', 'cohort1.ctls.95508', 'cohort1.ctls.95513', 'cohort1.ctls.95517', 'cohort1.ctls.95523', 'cohort1.ctls.95532', 'cohort1.ctls.95533', 'cohort1.ctls.95538', 'cohort1.ctls.95539', 'cohort1.ctls.95541', 'cohort1.ctls.95548', 'cohort1.ctls.95553', 'cohort1.ctls.95554', 'cohort1.ctls.95561', 'cohort1.ctls.95568', 'cohort1.ctls.95574', 'cohort1.ctls.95577', 'cohort1.ctls.95582', 'cohort1.ctls.95590', 'cohort1.ctls.95607', 'cohort1.ctls.95608', 'cohort1.ctls.95613', 'cohort1.ctls.95618', 'cohort1.ctls.95621', 'cohort1.ctls.95625', 'cohort1.ctls.95647', 'cohort1.ctls.95648', 'cohort1.ctls.95659', 'cohort1.ctls.95660', 'cohort1.ctls.95666', 'cohort1.ctls.95667', 'cohort1.ctls.95672', 'cohort1.ctls.95674', 'cohort1.ctls.95692', 'cohort1.ctls.95695', 'cohort1.ctls.95704', 'cohort2.ctls.163500', 'cohort2.ctls.163501', 'cohort2.ctls.163512', 'cohort2.ctls.163514', 'cohort2.ctls.163515', 'cohort2.ctls.163536', 'cohort2.ctls.163542', 'cohort2.ctls.163550', 'cohort2.ctls.163551', 'cohort2.ctls.163552', 'cohort2.ctls.163554', 'cohort2.ctls.163570', 'cohort2.ctls.163574', 'cohort2.ctls.163583', 'cohort2.ctls.163589', 'cohort2.ctls.163599', 'cohort2.ctls.163600', 'cohort2.ctls.163608', 'cohort2.ctls.163612', 'cohort2.ctls.163613', 'cohort2.ctls.163617', 'cohort2.ctls.163618', 'cohort2.ctls.163638', 'cohort2.ctls.163648', 'cohort2.ctls.163650', 'cohort2.ctls.163670', 'cohort2.ctls.163671', 'cohort2.ctls.163673', 'cohort2.ctls.163675', 'cohort2.ctls.163676', 'cohort2.ctls.163677', 'cohort2.ctls.163689', 'cohort2.ctls.163693', 'cohort2.ctls.163705', 'cohort2.ctls.163708', 'cohort2.ctls.163715', 'cohort2.ctls.163716']
for gene in genes_of_interest:
	##get all variants in gene
	for case in cases_list:
		filter_var_by_genename([gene], case + '.annotated.txt', case + '.' + gene + '.temp.xls', 15)
	##combine in one file
	combine_ann_txt(cases_list, '', '.' + gene + '.temp.xls', 'cases.' +  gene + '.all_vars.xls')
	filter_ann_txt_files('cases.' +  gene + '.all_vars.xls', filtering_cols, mafs_required, 'cases', pc_defs, exonic_defs, gene)
	for maf in mafs_required:
		##make dict for each proband with a gene var
		pc_vars = 'cases.' + str(maf) + ".protein_changing." + gene + ".xls"
		var_maf_dict = make_maf_dict_from_pc_vars(pc_vars, filtering_cols[0] - 1)
		# for p in var_maf_dict:
		# 	print p, var_maf_dict[p]
		##get other variants and make final file
		combine_all_data_make_final_file(var_maf_dict, '.' +  str(maf) + '.protein_changing.all_genes.xls', 'cases_1117.' + str(maf) + '.protein_changing.' + gene + '_associated.xls')
		##summarize final file
		summarize_final_file('cases_1117.' + str(maf) + '.protein_changing.' + gene + '_associated.xls', 'cases_1117.' + str(maf) + '.protein_changing.' + gene + '_associated.summary.xls')





