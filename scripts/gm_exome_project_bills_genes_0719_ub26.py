#!/usr/bin/env python
import sys
import os

##parameters
delim = '\t'

##working dir
working_dir = '/data/atimms/gm_exome_project_0719'
os.chdir(working_dir)

##methods
def gene_dict_from_file(infile):
	gene_dict = {}
	with open(infile, "r") as infh:
		line_count = 0
		for line in infh:
			line_count += 1
			if line_count > 1:
				line = line.rstrip().split(delim)
				gene = line[0]
				info = line[1:]
				if gene in gene_dict:
					gene_dict[gene].append(info)
				else:
					gene_dict[gene] = [info]
	return(gene_dict)

def ped_dict_from_file(infile):
	ped_dict = {}
	with open(infile, "r") as infh:
		line_count = 0
		for line in infh:
			line_count += 1
			if line_count > 1:
				line = line.rstrip().split(delim)
				ped = line[0]
				info = line[1:]
				if ped in ped_dict:
					print('ped %s seen multiple times'%ped)
				else:
					ped_dict[ped] = info
	return(ped_dict)

def get_vars_from_bills_list(ped_file, bills_genelist, infile_suffix):
	gene_dict = gene_dict_from_file(bills_genelist)
	ped_dict = ped_dict_from_file(ped_file)
	print(len(gene_dict), len(ped_dict))
	for pedigree in ped_dict:
		pedigree_info = ped_dict[pedigree]
		ped_type = pedigree_info[1]
		var_file = pedigree + infile_suffix
		var_file_bill = var_file.rsplit('.', 1)[0] + '.bl.xls'
		with open(var_file, "r") as infh, open(var_file_bill, "w") as outfh:
			line_count = 0
			for line in infh:
				line_count += 1
				if line_count == 1:
					outfh.write(line)
				else:
					line = line.split(delim)
					gene = line[5]
					if gene in gene_dict:
						outfh.write(delim.join(line))

def get_counts_from_var_file(var_file):
	with open(var_file, "r") as infh:
		line_count = 0
		for line in infh:
			line_count += 1
			if line_count == 1:
				outfh.write(line)
			else:
				line = line.split(delim)
				gene = line[5]
				if gene in gene_dict:
					outfh.write(delim.join(line))


def get_solved_numbers(ped_dict, outfile):
	counts_dict = {}
	for ped in ped_dict:
		dx_type = '_'.join([ped_dict[ped][2], ped_dict[ped][1]])
		status = ped_dict[ped][8]
		# print(ped, dx_type, status)
		if dx_type in counts_dict:
			counts_dict[dx_type].append(status)
		else:
			counts_dict[dx_type] = [status]
	all_values = sum(counts_dict.values(), [])
	groups = list(set(all_values))
	groups.sort()
	with open(outfile, "w") as outfh:
		header = ['DxGroup1', 'type', 'total peds'] + groups
		outfh.write(delim.join(header) + '\n')
		for dt in counts_dict:
			line_out = dt.split('_')
			line_out.append(str(len(counts_dict[dt])))
			print(dt, counts_dict[dt])
			for g in groups:
				count_g = counts_dict[dt].count(g)
				line_out.append(str(count_g))
			outfh.write(delim.join(line_out) + '\n')
			print(line_out)

def get_counts_from_var_file(in_file):
	anal_types = []
	ch_genes = []
	dn_vars = []
	##add var types to dict for each ped
	with open(in_file, "r") as infh:
		line_count = 0
		for line in infh:
			line_count += 1
			if line_count > 1:
				line = line.rstrip().split(delim)
				biotype = line[7]
				rmsk = line[29]
				segdup = line[30]
				gene = line[5]
				var = '_'.join(line[:5])
				##filter vars
				if biotype == 'protein_coding' and rmsk == 'None' and segdup == '0':
					anal_type = line[-1]
					##only count comp hets one per gene
					if anal_type == 'comp_hets':
						if gene not in ch_genes:
							anal_types.append(anal_type)
							ch_genes.append(gene)
					##no dup de novo vars
					elif anal_type == 'de_novo' or anal_type == 'x_linked_de_novo':
						if var not in dn_vars:
							anal_types.append(anal_type)
							dn_vars.append(var)
					else:
						anal_types.append(anal_type)
	return(anal_types)

def get_header_from_file(infile):
	with open(infile, "r") as infh:
		line_count = 0
		for line in infh:
			line_count += 1
			if line_count == 1:
				header = line.rstrip().split(delim)
	return(header)

def summarize_var_counts(infile, outfile):
	counts_dict, ped_count_dict = {}, {}
	with open(infile, "r") as infh:
		line_count = 0
		for line in infh:
			line_count += 1
			line = line.rstrip().split(delim)
			if line_count == 1:
				header = ['DxGroup1', 'type', 'total peds'] + line[12:]
			else:
				dx_type = '_'.join([line[3], line[2]])
				counts = line[12:]
				counts = [int(i) for i in counts]
				if dx_type in counts_dict:
					ped_count_dict[dx_type] += 1
					# print(dx_type, counts_dict[dx_type], counts)
					counts_dict[dx_type] = [counts_dict[dx_type][i] + counts[i] for i in range(len(counts_dict[dx_type]))]
					# print(dx_type, counts_dict[dx_type], counts)
				else:
					counts_dict[dx_type] = counts
					ped_count_dict[dx_type] = 1
	with open(outfile, "w") as outfh:
		outfh.write(delim.join(header) + '\n')
		for dt in counts_dict:
			line_out = dt.split('_')
			ped_count = ped_count_dict[dt]
			line_out.append(str(ped_count))
			print(line_out)
			counts = counts_dict[dt]
			print(counts)
			average_counts = [i / float(ped_count) for i in counts]
			print(average_counts)
			average_counts = [str(round(i, 2)) for i in average_counts]
			print(average_counts)
			line_out.extend(average_counts)
			outfh.write(delim.join(line_out) + '\n')
			print(line_out)


	# for ped in ped_dict:
	# 	dx_type = '_'.join([ped_dict[ped][2], ped_dict[ped][1]])
	# 	status = ped_dict[ped][8]
	# 	# print(ped, dx_type, status)
	# 	if dx_type in counts_dict:
	# 		counts_dict[dx_type].append(status)
	# 	else:
	# 		counts_dict[dx_type] = [status]
	# all_values = sum(counts_dict.values(), [])
	# groups = list(set(all_values))
	# groups.sort()
	# with open(outfile, "w") as outfh:
	# 	header = ['DxGroup1', 'type', 'total peds'] + groups
	# 	outfh.write(delim.join(header) + '\n')
	# 	for dt in counts_dict:
	# 		line_out = dt.split('_')
	# 		line_out.append(str(len(counts_dict[dt])))
	# 		print(dt, counts_dict[dt])
	# 		for g in groups:
	# 			count_g = counts_dict[dt].count(g)
	# 			line_out.append(str(count_g))
	# 		outfh.write(delim.join(line_out) + '\n')
	# 		print(line_out)


def get_var_numbers(ped_dict, ped_outfile, sum_outfile, varfile_suffix, start_header):
	counts_dict, info_dict = {}, {}
	##populate dict = [[ped_info], [var_types]]
	for pedigree in ped_dict:
		pedigree_info = ped_dict[pedigree]
		var_file = pedigree + varfile_suffix
		anal_counts = get_counts_from_var_file(var_file)
		# print(pedigree, len(anal_counts), anal_counts)
		if pedigree in counts_dict:
			print('pedigree %s seen multiple times!!', pedigree)
		else:
			##add ped info and analysis counts to file
			counts_dict[pedigree] = anal_counts
			info_dict[pedigree] = pedigree_info
	with open(ped_outfile, "w") as poutfh:
		##get list of all analysis types
		all_values = sum(counts_dict.values(), [])
		# print(all_values)
		groups = list(set(all_values))
		groups.sort()
		print(groups)
		header = start_header + groups
		poutfh.write(delim.join(header) + '\n')
		for p in info_dict:
			line_out = [p] + info_dict[p]
			for g in groups:
				count_g = counts_dict[p].count(g)
				line_out.append(str(count_g))
			poutfh.write(delim.join(line_out) + '\n')
			print(line_out)
	##take xounts by ped file and summrize by group
	summarize_var_counts(ped_outfile, sum_outfile)

def summarize_data(ped_file, bills_genelist, varfile_suffix, out_prefix):
	##bills list so gene:info
	gene_dict = gene_dict_from_file(bills_genelist)
	##ped info so pedid:info
	ped_dict = ped_dict_from_file(ped_file)
	print(len(gene_dict), len(ped_dict))
	##get ped numbers by dx and type
	'''
	solved_outfile = out_prefix + '.solved_counts.xls'
	get_solved_numbers(ped_dict, solved_outfile)
	'''
	##get var numbers i.e. b
	var_by_ped_file = out_prefix + '.var_counts.by_ped.xls'
	var_by_group_file = out_prefix + '.var_counts.summary.xls'
	pi_header = get_header_from_file(ped_file)
	get_var_numbers(ped_dict, var_by_ped_file, var_by_group_file, varfile_suffix, pi_header)













##run methods
ped_info = 'ped_info_0719.txt'
bills_genes = 'bills_genes.xls'
std_anal_suffix = '.std_analysis.xls'
results_prefix =  'gm_exome_project_0719'
##filter var file by presence in bills list
# get_vars_from_bills_list(ped_info, bills_genes, std_anal_suffix)

##summarize all data
summarize_data(ped_info, bills_genes, std_anal_suffix, results_prefix)




