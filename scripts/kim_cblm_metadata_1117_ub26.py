#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import filtering_annotated
import shutil

##set input variables and parameters
delim = '\t'


'''
##make dict from varaint files
var_dict = {}
##data upto 0317
with open('all_exome_data_std_pipeline.0317.xls', 'r') as vars_fh:
	lc = 0
	for line in vars_fh:
		lc += 1
		if lc >1:
			line = line.rstrip().split(delim)
			fam_id = line[20]
			analysis = line[24]
			gene_id = line[5]
			if gene_id in htseq_dict:
				is_expressed = htseq_dict[gene_id][0]
			else:
				is_expressed = 'yes'
			cadd = line[11]
			cadd_passed = 'no'
			if cadd == 'None':
				cadd_passed = 'yes'
			elif float(cadd) >= 10:
				cadd_passed = 'yes'
			##filter cadd >=10 and is expressed
			if cadd_passed == 'yes' and is_expressed == 'yes':
				if fam_id in var_dict:
					if analysis == 'de_novo' or analysis == 'x_linked_dn' or analysis == 'x_linked_de_novo':
						if line[0] == 'chrX':
							var_dict[fam_id][4] += 1
						else:
							var_dict[fam_id][0] += 1
					elif analysis == 'autosomal_recessive' or analysis == 'auto_rec':
						var_dict[fam_id][1] += 1
					elif analysis == 'compound_het' or analysis == 'comp_hets':
						var_dict[fam_id][5].append(gene_id)
					elif analysis == 'x_linked' or analysis == 'x_linked_recessive':
						var_dict[fam_id][3] += 1
					# else:
					# 	print 'analysis %s not counted'%analysis
				else:
					var_dict[fam_id] = [0,0,0,0,0,[]]
					if analysis == 'de_novo'or analysis == 'x_linked_dn' or analysis == 'x_linked_de_novo':
						if line[0] == 'chrX':
							var_dict[fam_id][4] += 1
						else:
							var_dict[fam_id][0] += 1
					elif analysis == 'autosomal_recessive' or analysis == 'auto_rec':
						var_dict[fam_id][1] += 1
					elif analysis == 'compound_het' or analysis == 'comp_hets':
						var_dict[fam_id][5].append(gene_id)
					elif analysis == 'x_linked' or analysis == 'x_linked_recessive':
						var_dict[fam_id][3] += 1
					# else:
					# 	print 'analysis %s not counted'%analysis
##add new data (slightly different format)
for var_data in ['LR14-071.std_analysis.xls', 'LR14-221.std_analysis.xls', 'LR16-079.std_analysis.xls', 'LR16-451.std_analysis.xls']:
	with open(var_data, 'r') as vars_fh:
		lc =0
		for line in vars_fh:
			lc += 1
			if lc >1:
				line = line.rstrip().split(delim)
				fam_id = line[28]
				analysis = line[33]
				gene_id = line[5]
				if gene_id in htseq_dict:
					is_expressed = htseq_dict[gene_id][0]
				else:
					is_expressed = 'yes'
				cadd = line[12]
				cadd_passed = 'no'
				if cadd == 'None':
					cadd_passed = 'yes'
				elif float(cadd) >= 10:
					cadd_passed = 'yes'
				##filter cadd >=10 and is expressed
				if cadd_passed == 'yes' and is_expressed == 'yes':
					if fam_id in var_dict:
						if analysis == 'de_novo' or analysis == 'x_linked_dn' or analysis == 'x_linked_de_novo':
							if line[0] == 'chrX':
								var_dict[fam_id][4] += 1
							else:
								var_dict[fam_id][0] += 1
						elif analysis == 'autosomal_recessive' or analysis == 'auto_rec':
							var_dict[fam_id][1] += 1
						elif analysis == 'compound_het' or analysis == 'comp_hets':
							var_dict[fam_id][5].append(gene_id)
						elif analysis == 'x_linked' or analysis == 'x_linked_recessive':
							var_dict[fam_id][3] += 1
						# else:
						# 	print 'analysis %s not counted'%analysis
					else:
						var_dict[fam_id] = [0,0,0,0,0,[]]
						if analysis == 'de_novo'or analysis == 'x_linked_dn' or analysis == 'x_linked_de_novo':
							if line[0] == 'chrX':
								var_dict[fam_id][4] += 1
							else:
								var_dict[fam_id][0] += 1
						elif analysis == 'autosomal_recessive' or analysis == 'auto_rec':
							var_dict[fam_id][1] += 1
						elif analysis == 'compound_het' or analysis == 'comp_hets':
							var_dict[fam_id][5].append(gene_id)
						elif analysis == 'x_linked' or analysis == 'x_linked_recessive':
							var_dict[fam_id][3] += 1
						# else:
						# 	print 'analysis %s not counted'%analysis
'''

##methods
def combine_all_denovo_snps(in_files, out_file):
	with open(out_file, "w") as outfh:
		# header = []
		# outfh.write(delim.join(header))
		for in_file in in_files:
			with open(in_file, "r") as infh:
				lc = 0
				for line in infh:
					lc +=1
					line = line.rstrip().split(delim)
					if lc > 1:
						if in_file.split('.')[0] == 'all_exome_data_std_pipeline':
							fam_id = line[20]
							analysis = line[24]
							gene_id = line[5]
						else:
							fam_id = line[28]
							analysis = line[33]
							gene_id = line[5]
						var = line[:5]
						ref = line[3]
						alt = line[4]
						if analysis == 'de_novo'or analysis == 'x_linked_dn' or analysis == 'x_linked_de_novo':
							if len(ref) == 1 and len(alt) ==1:
								outfh.write(delim.join(var + [gene_id, fam_id, analysis, '\n']))

def make_var_table(ped_dict, infile, outfile):
	count_dict = {}
	for group in ped_dict:
		peds = ped_dict[group]
		with open(infile, "r") as infh:
			count_dict[group] = [0,0,0,0,0,0,]
			for line in infh:
				line = line.rstrip().split(delim)
				ped = line[6]
				if ped in peds:
					snp_type = line[3] + line[4]
					# print ped, snp_type
					if snp_type == 'AG' or snp_type =='TC':
						count_dict[group][0] += 1
					elif snp_type == 'CT' or snp_type =='GA':
						count_dict[group][1] += 1
					elif snp_type == 'CA' or snp_type =='GT':
						count_dict[group][2] += 1
					elif snp_type == 'CG' or snp_type =='GC':
						count_dict[group][3] += 1
					elif snp_type == 'AC' or snp_type =='TG':
						count_dict[group][4] += 1
					elif snp_type == 'AT' or snp_type =='TA':
						count_dict[group][5] += 1				
					else:
						print 'snp type %s not recognized'%snp_type
	with open(outfile, "w") as outfh:
		header = ['Group', 'AG_TC', 'CT_GA', 'CA_GT', 'CG_GC', 'AC_TG', 'AT_TA', '\n']
		outfh.write(delim.join(header))
		for g in count_dict:
			print g, count_dict[g]
			counts = [str(i) for i in count_dict[g]]
			outfh.write(delim.join([g] + counts + ['\n']))



##get denovo variants for groups of samples
##cblh vs dwm and non gx vs gx
##count transistions
##make table
project_name = 'cblm_metadata_1117'
dwm_peds = ['LR04-186', 'LR05-203', 'LR05-354', 'LR12-115', 'LR12-434', 'LR12-443', 'LR03-274', 'LR06-085', 'LR03-077', 'LR03-223', 'LR04-017', 'LR04-084', 'LR05-118', 'LR01-079', '3C-4', 'LR06-157', 'LR11-152', 'LR06-278', 'LR08-002', 'LR03-278', 'LR03-298', 'LR12-313a2', 'LR04-106', 'LR13-315', 'LR03-332', 'LR04-233', 'LR04-399_2', 'LR10-222', 'LR09-023']
cblm_peds = ['LR05-160', 'LR10-102', 'LR11-033', 'LR03-305', 'LR08-323', 'LR08-390', 'LR13-002', 'LR13-153', 'LR03-206', 'LR04-020', 'LR04-371', 'LR09-280', 'LR03-169', 'LR04-208', 'LR05-396', 'LR14-221', 'LR16-451', 'LR12-439_2', 'LR10-199', 'LR04-376', 'LR03-130', 'LR12-316', 'LR06-207', 'LR10-026', 'LR10-228a1', 'LR11-042', 'LR08-396', 'LR11-169', 'LR12-032', 'LR12-463', 'LR03-055', 'LR05-398', 'LR14-071', 'LR16-079', 'LR03-120', 'LR05-120', 'LR10-243', 'LR10-230', 'LR12-464', 'LR13-037', 'LR13-085', 'LR13-199', 'LR13-200', 'LR09-227', 'LR11-241', 'LR05-007']
gx_ped = ['LR04-186', 'LR05-203', 'LR05-354', 'LR12-115', 'LR12-434', 'LR12-443', 'LR03-274', 'LR06-085', 'LR03-077', 'LR03-223', 'LR04-017', 'LR04-084', 'LR05-118', 'LR01-079', '3C-4', 'LR06-157', 'LR11-152', 'LR06-278', 'LR08-002', 'LR05-160', 'LR10-102', 'LR11-033', 'LR03-305', 'LR08-323', 'LR08-390', 'LR13-002', 'LR13-153', 'LR03-206', 'LR04-020', 'LR04-371', 'LR09-280', 'LR03-169', 'LR04-208', 'LR05-396', 'LR14-221', 'LR16-451', 'LR12-439_2', 'LR10-199', 'LR04-376', 'LR03-130', 'LR12-316', 'LR06-207', 'LR10-026', 'LR10-228a1', 'LR11-042']
nongx_ped = ['LR03-278', 'LR03-298', 'LR12-313a2', 'LR08-396', 'LR11-169', 'LR12-032', 'LR12-463', 'LR03-055', 'LR05-398', 'LR14-071', 'LR16-079', 'LR03-120', 'LR05-120', 'LR10-243']
var_files = ['all_exome_data_std_pipeline.0317.xls', 'LR14-071.std_analysis.xls', 'LR14-221.std_analysis.xls', 'LR16-079.std_analysis.xls', 'LR16-451.std_analysis.xls']
ped_dict = {'dwm':dwm_peds, 'cblm':cblm_peds, 'gx':gx_ped, 'nongx': nongx_ped}

denovo_vars = project_name + '.denovo_vars.xls'
final_table = project_name + '.variant_counts.xls'

##combine all denovo vars
combine_all_denovo_snps(var_files, denovo_vars)
make_var_table(ped_dict, denovo_vars, final_table)

