#!/usr/bin/env python
import os
import subprocess
import glob
import shutil

##set input variables and parameters
delim = '\t'
working_dir = '/data/atimms/microtia_exomes/batch1-5_0117'

##programs and ref files
vt = '/data/atimms/gemini/vt/vt'
gemini = '/data/atimms/gemini/bin/gemini'
fasta = '/data/atimms/references/human_g1k_v37.fasta'
snpeff_jar = '/data/atimms/gemini/snpEff/snpEff.jar'
bcftools_12 = '/home/atimms/programs/bcftools-1.2/bcftools'
# dbdb_data = '/data/atimms/references/gemini/dbdb.gene_association.txt'
gemini_ref_dir = '/data/atimms/references/gemini/'
rvis_data = gemini_ref_dir + 'RVIS_ExAC_4KW.txt'
gdi_data = gemini_ref_dir + 'GDI_full_10282015.txt'
mgi_data = gemini_ref_dir + 'mgi.abnormal_outer_ear_morhology.txt'
hpo_data = gemini_ref_dir + 'genes_for_HP_0000356.csv'
sbse_rnaseq_data = gemini_ref_dir + 'esra.gene_exp.diff'
hmx1_rnaseq_data = gemini_ref_dir + 'jess.gene_exp.diff'
human_ear_rnaseq_data = gemini_ref_dir + 'human_ear.gene_exp.diff'
human_syndrome_data = gemini_ref_dir + 'microtia_human_syndromes.txt'
hoxa2_data = gemini_ref_dir + 'hoxa2.expression_data.txt'
ba_exp_data = gemini_ref_dir + 'ba_expression_papers_0816.txt'

##all pedigrees, all parents coded as unaffected
ped_dn_rec = 'microtia_exomes.dn_recessive.0117.ped'
##all peds, all unaffected coded as unknown
ped_dom = 'microtia_exomes.dominant.0117.ped'
##all peds, all carriers and ptags coded as unknown
ped_real = 'microtia_exomes.real.0117.ped'

##annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/data/atimms/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene']
av_operation = ['-operation', 'g']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-arg', '-splicing 10 ']
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'


def load_single_vcf_into_gemini(ped, db_name, input_vcf, ped_file, anal_type):
	##dcompress, change a gatk header thing, and decompse and normalize (should be redundant)
	temp_vcf = 'temp1.vcf'
	normalized_vcf = ped + '.' + anal_type + '.int.norm.vcf'
	# ped_file = ped + '.ped'
	with open (temp_vcf, 'w') as tvcf_fh:
		zless_vcf = subprocess.Popen(['zless', input_vcf], stdout=subprocess.PIPE)
		sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
		vt_decompose = subprocess.Popen([vt, 'decompose', '-s', '-'], stdin=sed_vcf.stdout, stdout=subprocess.PIPE)
		vt_normalize = subprocess.Popen([vt, 'normalize', '-r', fasta, '-'], stdin=vt_decompose.stdout, stdout=tvcf_fh)
		vt_normalize.wait()
	##annotate with snpeff and compress and index
	with open (normalized_vcf, 'w') as nvcf_fh:
		snpeff_vcf = subprocess.Popen(['java', '-Xmx20G', '-jar', snpeff_jar, 'GRCh37.75', '-v', '-formatEff', '-classic', temp_vcf], stdout=nvcf_fh)
		snpeff_vcf.wait()
		bgzip_vcf = subprocess.Popen(['bgzip', normalized_vcf])
		bgzip_vcf.wait()
		tabix_vcf = subprocess.Popen(['tabix', '-p', 'vcf', normalized_vcf + '.gz'])
		tabix_vcf.wait()
	##add file to gemini, requires ped file
	gemini_load = subprocess.Popen([gemini, 'load', '--cores', '10', '-t', 'snpEff', '-p', ped_file, '-v', normalized_vcf + '.gz', db_name])
	gemini_load.wait()
	##remove intermediate files
	os.remove(temp_vcf)

def cp_gemini_db(in_db, out_db):
	cp_command = subprocess.Popen(['cp', in_db, out_db])
	cp_command.wait()


def gemini_comp_het(ped, gem_db_name, max_aaf, coverage, anal_type):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.' + anal_type + '.comp_hets.xls'
	with open (outfile, 'w') as out_fh:
		##compound het: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'comp_hets', '--max-priority', '2', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter',  "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf) , '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

def gemini_recessive(ped, gem_db_name, max_aaf, coverage, anal_type):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.' + anal_type + '.autosomal_recessive.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
# 		gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--columns', 'chrom,start,end,ref,alt,gene,biotype,impact,aa_change,gerp_bp_score,cadd_scaled,max_aaf_all', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '-d', '10', gem_db_name], stdout=out_fh)

		gemini_query.wait()
	return outfile

def gemini_denovo(ped, gem_db_name, max_aaf, coverage, anal_type):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.' + anal_type + '.de_novo.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'de_novo', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

##not used atm
def gemini_dominant(ped, gem_db_name, max_aaf, coverage, anal_type):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.' + anal_type + '.autosomal_dominant.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'autosomal_dominant', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

def gemini_inherited(ped, gem_db_name, max_aaf, coverage, anal_type):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.' + anal_type + '.inher.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'autosomal_dominant', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile


def gemini_xlinked(ped, gem_db_name, max_aaf, coverage, anal_type):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.' + anal_type + '.x_linked_recessive.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'x_linked_recessive', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

def gemini_xlinked_de_novo(ped, gem_db_name, max_aaf, coverage, anal_type):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.' + anal_type + '.x_linked_de_novo.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'x_linked_de_novo', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

def gemini_de_novo_syn(ped, gem_db_name, max_aaf, coverage, anal_type):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.' + anal_type + '.de_novo_syn.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'de_novo', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity == 'LOW' and is_coding == 1 and gerp_bp_score >= 2 AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

def gemini_recessive_syn(ped, gem_db_name, max_aaf, coverage, anal_type):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.' + anal_type + '.auto_rec_syn.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity == 'LOW' and is_coding == 1 and gerp_bp_score >= 2 AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

def make_var_dict_from_gemini_result(gemini_file):
	var_dict = {}
	line_count = 0
	with open(gemini_file, 'r') as in_fh:
		for line in in_fh:
			line = line.split(delim)
			line_count += 1
			if line_count > 1:
				var = '_'.join(line[:3] + line[6:8])
				# print var
				var_dict[var] = 1
	return var_dict

def gemini_potential_cnv(ped, gem_db_name, max_aaf, coverage, anal_type):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.' + anal_type + '.potential_cnv.xls'
	all_auto_rec_file = 'all_auto_rec.temp.xls'
	correct_auto_rec_file = 'correct_auto_rec.temp.xls'
	with open (correct_auto_rec_file, 'w') as car_fh:
		gemini_query_a = subprocess.Popen([gemini, 'autosomal_recessive', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=car_fh)
		gemini_query_a.wait()
	with open (all_auto_rec_file, 'w') as aar_fh:
		gemini_query_b = subprocess.Popen([gemini, 'autosomal_recessive', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=aar_fh)
		gemini_query_b.wait()
	correct_cnv_var_dict = make_var_dict_from_gemini_result(correct_auto_rec_file)
	with open (outfile, 'w') as out_fh:
		with open(all_auto_rec_file, 'r') as all_ar_fh:
			line_count = 0
			for line in all_ar_fh:
				line = line.split(delim)
				line_count += 1
				if line_count == 1:
					out_fh.write(delim.join(line))
				else:
					var = '_'.join(line[:3] + line[6:8])
					chrom  = line[0]
					if chrom != 'chrX':
						if var not in correct_cnv_var_dict:
							out_fh.write(delim.join(line))
	##remove intermediate files
	os.remove(all_auto_rec_file)
	os.remove(correct_auto_rec_file)
	return outfile

##make dict using exac rvis file
def make_dict_from_exac_intolerence(intol_file):
	intol_dict = {}
	with open(intol_file, "U") as intol:
		line_count = 0
		for line in intol:
			line = line.strip('\n').split(delim)
			line_count += 1
			gene_list = line[0:5]
			genes = list(set(gene_list))
			# print genes, len(genes)
			rvis = line[5]
			percentile = line[6]
			# print genes, rvis, percentile
			for gene in genes:
				if gene not in intol_dict:
					intol_dict[gene] = [rvis, percentile]
					# print gene, rvis, percentile
				else:
					print 'gene %s seen multiple times'% gene
	return intol_dict

##make dict using gdi file
def make_dict_from_gdi(gdi_file):
	gdi_dict = {}
	with open(gdi_file, "U") as gdi_fh:
		for line in gdi_fh:
			line = line.strip('\n').split(delim)
			# print line
			gene = line[0]
			gdi = line[1]
			gdi_phred = line[2]
			if gene not in gdi_dict:
				gdi_dict[gene] = [gdi, gdi_phred]
				# print gene, gdi, gdi_phred
			else:
				print 'gene %s seen multiple times'% gene
	return gdi_dict

def add_rvis_and_GDI(in_file, out_file, dicts_to_add, gene_col, position_to_insert):
	with open(in_file, "U") as inf, open(out_file, "w") as outf:
		line_count, genes_not_found = 0, 0
		for line in inf:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				extra_header = ['RVIS exac', 'Percentile exac', 'GDI', 'GDI phred', 'mgi phenotype', 'hpo_syndrome', 'sbse_wt_fpkm', 'sbse_log2_fc', 'sbse_q_value', 
						'hmx1_wt_fpkm', 'hmx1_log2_fc', 'hmx1_q_value', 'ear_57d_fpkm', 'ear_59d_fpkm', 'human_syndrome', 'hoxa2_wt1', 'hoxa2_wt3', 'ba_expression_papers_0816']
				outf.write(delim.join(line[:position_to_insert] + extra_header + line[position_to_insert:] + ['\n']))
			else:
				gene = line[gene_col]
				# print gene
				##get values for all gene names and print minimum value
				extra_stuff = []
				for data_dict in dicts_to_add:
					if gene in data_dict:
						extra_stuff.extend(data_dict[gene])
					else:
						extra_stuff.extend(['na']*len(data_dict[data_dict.keys()[0]]))

				# print gene, extra_stuff
				outf.write(delim.join(line[:position_to_insert] + extra_stuff + line[position_to_insert:] + ['\n']))

##make dict using mgi file
def make_dict_from_mgi(mgi_file):
	mgi_dict = {}
	mgi_dict2 = {}
	with open(mgi_file, "U") as mgi_fh:
		for line in mgi_fh:
			line = line.strip('\n').split(delim)
			# print line
			gene = line[2]
			pheno = line[3]
			if gene not in mgi_dict:
				mgi_dict[gene] = [pheno]
				# print gene, inheritance, pheno, syndrome, loe
			else:
				mgi_dict[gene].append(pheno)
	for g in mgi_dict:
		mgi_dict2[g] = [', '.join(mgi_dict[g])]
	# for g in mgi_dict:
	# 	print g, mgi_dict[g], mgi_dict2[g]
	return mgi_dict2

def make_dict_from_hpo(hpo_file):
	hpo_dict = {}
	with open(hpo_file, "U") as hpo_fh:
		line_count = 0
		for line in hpo_fh:
			line_count += 1
			if line_count > 2:
				line = line.strip('\n').split(',',1)

				gene = line[0].strip('"').split(' ')[0]
				syndrome = line[1].strip('"')
				# print gene, syndrome
				if gene not in hpo_dict:
					hpo_dict[gene] = [syndrome]
					# print gene, gdi, gdi_phred
				else:
					print 'gene %s seen multiple times'% gene
	return hpo_dict

def make_dict_from_cufflinks_rnaseq(gene_diff_file, cols_req):
	rnaseq_dict = {}
	with open(gene_diff_file, "U") as gdf_fh:
		for line in gdf_fh:
			line = line.strip('\n').split(delim)
			# print line
			gene = line[2]
			data = []
			for col in cols_req:
				d = line[col]
				data.append(d)
			# print gene, data
			if gene != '-':
				if gene not in rnaseq_dict:
					rnaseq_dict[gene] = data
					# print gene, gdi, gdi_phred
				else:
					# print 'gene %s seen multiple times'% gene
					pass
		return rnaseq_dict

def make_dict_from_human_syndrome(human_syndrome_file):
	hum_syn_dict = {}
	with open(human_syndrome_file, "U") as hs_fh:
		line_count = 0
		for line in hs_fh:
			line_count += 1
			if line_count > 1:
				line = line.strip('\n').split(delim)
				gene = line[0]
				syndrome = line[1]
				# print gene, syndrome
				if gene not in hum_syn_dict:
					hum_syn_dict[gene] = [syndrome]
				else:
					# print 'gene %s seen multiple times'% gene
					pass
	return hum_syn_dict

def make_dict_from_hoxa2(hoxa2_file):
	hoxa2_dict = {}
	with open(hoxa2_file, "U") as hoxa2_fh:
		line_count = 0
		for line in hoxa2_fh:
			line_count += 1
			if line_count > 1:
				line = line.strip('\n').split(delim)
				gene = line[1]
				data = line[2:4]
				# print gene, data
				if gene not in hoxa2_dict:
					hoxa2_dict[gene] = data
				else:
					# print 'gene %s seen multiple times'% gene
					pass
	return hoxa2_dict

def make_dict_from_ba_expression_papers(ba_exp_file):
	ba_dict = {}
	ba_dict2 = {}
	with open(ba_exp_file, "U") as ba_fh:
		for line in ba_fh:
			line = line.strip('\n').split(delim)
			# print line
			gene = line[0]
			paper = line[1]
			if gene not in ba_dict:
				ba_dict[gene] = [paper]
				# print gene, inheritance, pheno, syndrome, loe
			else:
				ba_dict[gene].append(paper)
	for g in ba_dict:
		ba_dict2[g] = [', '.join(ba_dict[g])]
	# for g in mgi_dict:
	# 	print g, mgi_dict[g], mgi_dict2[g]
	return ba_dict2

def add_columns_to_file(in_file, out_file, dicts_to_add, gene_col, position_to_insert):
	with open(in_file, "U") as inf, open(out_file, "w") as outf:
		line_count, genes_not_found = 0, 0
		for line in inf:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				extra_header = ['RVIS exac', 'Percentile exac', 'GDI', 'GDI phred', 'mgi phenotype', 'hpo_syndrome', 'sbse_wt_fpkm', 'sbse_log2_fc', 'sbse_q_value', 
						'hmx1_wt_fpkm', 'hmx1_log2_fc', 'hmx1_q_value', 'ear_57d_fpkm', 'ear_59d_fpkm', 'human_syndrome', 'hoxa2_wt1', 'hoxa2_wt3', 'ba_expression_papers_0816']
				outf.write(delim.join(line[:position_to_insert] + extra_header + line[position_to_insert:] + ['\n']))
			else:
				gene = line[gene_col]
				# print gene
				##get values for all gene names and print minimum value
				extra_stuff = []
				for data_dict in dicts_to_add:
					if gene in data_dict:
						extra_stuff.extend(data_dict[gene])
					else:
						extra_stuff.extend(['na']*len(data_dict[data_dict.keys()[0]]))

				# print gene, extra_stuff
				outf.write(delim.join(line[:position_to_insert] + extra_stuff + line[position_to_insert:] + ['\n']))


def combine_gemini_results(ped, file_list, gene_col, position_to_insert, position_to_finish, add_annovar_pos, out_suffix):
	##file names
	outfile = ped + out_suffix
	temp_file = 'combined.temp.xls'
	temp2_file = 'combined.temp2.xls'
	print 'combining files:', file_list
	##combine all files, making sure col lengths are all the same
	with open(temp_file, "w") as temp_fh:
		file_count, total_line_count = 0, 0
		for filename in file_list:
			# print filename
			analysis_type = filename.split('.')[-2]
			# file_count += 1
			# print analysis_type, file_count
			##add header from first file
			with open(filename, "r") as in_fh:
				line_count = 0
				for line in in_fh:
					line_count +=1
					total_line_count += 1
					line = line.rstrip().split(delim)
					if line_count == 1:
						file_count += 1
						if file_count == 1:
							header = delim.join(line[:position_to_insert] + line[position_to_insert:position_to_finish] + ['analysis', '\n'])
							temp_fh.write(header)
					else:
						line_out = delim.join(line[:position_to_insert] + line[position_to_insert:position_to_finish] + [analysis_type, '\n'])
						temp_fh.write(line_out)
	if total_line_count >=2:
		##add dbdb and mgi data
		##add rvis, gdi, mgi data
		rvis_data_dict = make_dict_from_exac_intolerence(rvis_data)
		gdi_data_dict = make_dict_from_gdi(gdi_data)
		mgi_data_dict = make_dict_from_mgi(mgi_data)
		hpo_data_dict = make_dict_from_hpo(hpo_data)
		sbse_rnaseq_dict = make_dict_from_cufflinks_rnaseq(sbse_rnaseq_data, [8,9,11])
		hmx1_rnaseq_dict = make_dict_from_cufflinks_rnaseq(hmx1_rnaseq_data, [8,9,11])
		human_ear_rnaseq_dict = make_dict_from_cufflinks_rnaseq(human_ear_rnaseq_data, [7,8])
		human_syndrome_dict = make_dict_from_human_syndrome(human_syndrome_data)
		hoxa2_dict = make_dict_from_hoxa2(hoxa2_data)
		ba_paper_dict = make_dict_from_ba_expression_papers(ba_exp_data)
		add_columns_to_file(temp_file, temp2_file, [rvis_data_dict, gdi_data_dict, mgi_data_dict, hpo_data_dict, sbse_rnaseq_dict, hmx1_rnaseq_dict, human_ear_rnaseq_dict, human_syndrome_dict, hoxa2_dict, ba_paper_dict], gene_col, position_to_insert)
		# add_columns_to_file(temp_file, temp2_file, dbdb_data_dict, mgi_data_dict, omim_data_dict, gene_col, position_to_insert)
		add_annovar_ts_info(ped, temp2_file, outfile, add_annovar_pos)
		##remove intermediate files
		os.remove(temp_file)
		os.remove(temp2_file)

def add_annovar_ts_info(pedigree, in_file, out_file, pos_to_insert):
	##files
	normalized_vcf = pedigree + '.int.norm.vcf.gz'
	region_file = pedigree + '_temp.regions'
	temp_vcf = pedigree + 'temp0.vcf.gz'
	av_file = pedigree + '.avinput'
	multianno = pedigree + '.hg19_multianno.txt'
	## make_list_of_regions(var_file, region_file):
	with open(in_file, 'r') as var_fh, open(region_file, 'w') as reg_fh:
		line_count = 0
		for line in var_fh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count >1:
				chrom = line[0][3:]
				# start = int(line[1]) - 2
				# end = int(line[2]) + 2
				start = line[1]
				end = line[2]
				reg_fh.write(delim.join([chrom, str(start), str(end), '\n']))
				# print chrom, start_end
	## get_var_from_regions_in_vcf(in_vcf, regions, out_vcf):
	bcftools_filter = subprocess.Popen([bcftools_12, 'view', '-R', region_file, '-o', temp_vcf, '-O', 'z', normalized_vcf])
	bcftools_filter.wait()
	##convert vcf file to individual avinput file
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', temp_vcf, '-allsample', '--withfreq', '-outfile', av_file])
	con_ann.wait()
	##run_table_annovar(avinput, av_prefix):
	command = [table_annovar] + av_buildver + [av_file] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', pedigree]
	annovar = subprocess.Popen(command)
	annovar.wait()
	## add_annovar_info(var_file, multianno, outfile, pos_to_insert):
	##make dict with annovar data
	ann_dict = {}
	with open(multianno, 'r') as multi_fh:
		for line in multi_fh:
			line = line.rstrip().split(delim)
			chrom = line[0]
			pos = line[2]
			ch_pos = '_'.join([chrom,pos])
			ref = line[3]
			alt = line[4]
			ts_info = line[9]
			if ts_info == '.':
				ts_info = line[7]
			# print chrom, pos, ts_info
			ann_dict[ch_pos] = ts_info

	with open(in_file, 'r') as in_fh, open(out_file, 'w') as out_fh:
		line_count = 0
		for line2 in in_fh:
			line2 = line2.split(delim)
			line_count += 1
			if line_count == 1:
				out_fh.write(delim.join(line2[:pos_to_insert] + ['annovar_ts_info'] + line2[pos_to_insert:]))
			else:
				chrom2 = line2[0][3:]
				pos2 = line2[2]
				ch_pos2 = '_'.join([chrom2,pos2])
				if ch_pos2 in ann_dict:
					# print line2[:5],ch_pos2, ann_dict[ch_pos2]
					out_fh.write(delim.join(line2[:pos_to_insert] + [ann_dict[ch_pos2]] + line2[pos_to_insert:]))
				else:
					print 'not found:', line2[:5], ch_pos2

def add_new_ped_file_to_gemini_db(pedfile, geminidb):
	gemini_amend = subprocess.Popen([gemini, 'amend', '--sample', pedfile, geminidb])
	gemini_amend.wait()

##run all other methods
def modified_gemini_protocol(working_directory, pedigree, ped_type):
	os.chdir(working_directory)
	for analysis in ['freebayes', 'gatk', 'intesect']:
		gemini_db = pedigree +  '.' +  analysis + '_gemini.db'
		if analysis == 'intesect':
			in_vcf = pedigree + '.intersected_vcfs/0002.vcf.gz'
		elif analysis == 'gatk':
			in_vcf = pedigree + '.intersected_vcfs/0000.vcf'
		elif analysis == 'freebayes':
			in_vcf = pedigree + '.intersected_vcfs/0001.vcf'

		freq_req = 0.01
		freq_req_recessive = 0.05
		coverage_req = 5
		##load vcf file into db
		load_single_vcf_into_gemini(pedigree, gemini_db, in_vcf, ped_dn_rec, analysis)
		##run queries and get files to combine
		files_to_combine = []
		if ped_type == 'singleton' or ped_type == 'duo':
			files_to_combine.append(gemini_comp_het(pedigree, gemini_db, freq_req_recessive, coverage_req,analysis))
			files_to_combine.append(gemini_recessive(pedigree, gemini_db, freq_req_recessive, coverage_req,analysis))
			# gemini_dominant(pedigree, gemini_db, freq_req, coverage_req)
		elif ped_type == 'trio':
			files_to_combine.append(gemini_comp_het(pedigree, gemini_db, freq_req_recessive, coverage_req,analysis))
			files_to_combine.append(gemini_recessive(pedigree, gemini_db, freq_req_recessive, coverage_req,analysis))
			files_to_combine.append(gemini_denovo(pedigree, gemini_db, freq_req, coverage_req,analysis))
			files_to_combine.append(gemini_xlinked(pedigree, gemini_db, freq_req_recessive, coverage_req,analysis))
			files_to_combine.append(gemini_xlinked_de_novo(pedigree, gemini_db, freq_req, coverage_req,analysis))
			files_to_combine.append(gemini_de_novo_syn(pedigree, gemini_db, freq_req, coverage_req,analysis))
			files_to_combine.append(gemini_recessive_syn(pedigree, gemini_db, freq_req_recessive, coverage_req,analysis))
			files_to_combine.append(gemini_potential_cnv(pedigree, gemini_db, freq_req_recessive, coverage_req,analysis))
			# gemini_dominant(pedigree, gemini_db, freq_req, coverage_req)
		else:
			print 'ped type %s not recognized'%ped_type
		##combine results
		combine_gemini_results(pedigree, files_to_combine, 5, 13, 22, 6, '.' + analysis + '_analysis.xls')
		##get auto dom vars regardless of family history
		add_new_ped_file_to_gemini_db(ped_dom, gemini_db)
		gemini_dominant(pedigree, gemini_db, freq_req, coverage_req,analysis)
		combine_gemini_results(pedigree, [pedigree + '.cov' + str(coverage_req) + '.maf'+ str(freq_req) + '.' + analysis + '.autosomal_dominant.xls'], 5, 13, 22, 6, '.' + analysis + '.auto_dom.xls')
		##get inherited vars - should be same as auto dom unless is family history
		add_new_ped_file_to_gemini_db(ped_real, gemini_db)
		gemini_inherited(pedigree, gemini_db, freq_req, coverage_req,analysis)
		combine_gemini_results(pedigree, [pedigree + '.cov' + str(coverage_req) + '.maf'+ str(freq_req) + '.' + analysis +  '.inher.xls'], 5, 13, 22, 6, '.' + analysis + '.inherited.xls')

def combine_files_and_filter_by_gene_list(outfile_prefix, infile_suffix, genelist, genename_col):
	outfile = outfile_prefix + infile_suffix.rsplit('.', 1)[0] + '.txt' 
	infiles = glob.glob('*' + infile_suffix)
	# print outfile
	with open(outfile, 'w') as out_fh: 
		file_count = 0
		for infile in infiles:
			# print infile
			file_count += 1
			line_count = 0
			with open(infile, 'r') as in_fh:
				for line in in_fh:
					line = line.rstrip().split(delim)
					line_count += 1
					if line_count == 1:
						if file_count == 1:
							out_fh.write(delim.join(line + ['\n']))
					else:
						genename = line[genename_col]
						# print analysis_type
						if genename in genelist:
							out_fh.write(delim.join(line + ['\n']))



##run methods

##parameters
rho_genes = ['A2M', 'ABR', 'AKAP13', 'ARAP1', 'ARAP2', 'ARAP3', 'ARHGAP1', 'ARHGAP10', 'ARHGAP11A', 'ARHGAP11B', 'ARHGAP12', 'ARHGAP15', 'ARHGAP17', 'ARHGAP18', 
		'ARHGAP19', 'ARHGAP20', 'ARHGAP21', 'ARHGAP22', 'ARHGAP23', 'ARHGAP24', 'ARHGAP25', 'ARHGAP26', 'ARHGAP27', 'ARHGAP28', 'ARHGAP29', 'ARHGAP30', 'ARHGAP31', 
		'ARHGAP32', 'ARHGAP33', 'ARHGAP35', 'ARHGAP36', 'ARHGAP39', 'ARHGAP4', 'ARHGAP40', 'ARHGAP42', 'ARHGAP44', 'ARHGAP5', 'ARHGAP6', 'ARHGAP8', 'ARHGAP9', 
		'ARHGDIA', 'ARHGDIB', 'ARHGDIG', 'ARHGEF1', 'ARHGEF10', 'ARHGEF10L', 'ARHGEF11', 'ARHGEF12', 'ARHGEF15', 'ARHGEF16', 'ARHGEF17', 'ARHGEF18', 'ARHGEF19', 
		'ARHGEF2', 'ARHGEF26', 'ARHGEF3', 'ARHGEF33', 'ARHGEF35', 'ARHGEF37', 'ARHGEF38', 'ARHGEF39', 'ARHGEF4', 'ARHGEF40', 'ARHGEF5', 'ARHGEF6', 'ARHGEF7', 
		'ARHGEF9', 'BCR', 'CDC42', 'CHN1', 'CHN2', 'DEPDC1B', 'DEPDC7', 'DLC1', 'ECT2', 'FAM13A', 'FAM13B', 'FGD1', 'FGD2', 'FGD3', 'FGD4', 'GDI1', 'GDI2', 
		'GDP', 'GMIP', 'GNA13', 'GTP', 'HMHA1', 'ITSN1', 'KALRN', 'MCF2L', 'MYO9A', 'MYO9B', 'NET1', 'NGEF', 'OBSCN', 'OCRL', 'OPHN1', 'PIK3R2', 'PLEKHG2', 
		'PLEKHG5', 'PREX1', 'RAC1', 'RAC2', 'RAC3', 'RACGAP1', 'RALBP1', 'RASGRF2', 'RHOA', 'RHOB', 'RhoBTB', 'RHOC', 'RHOD', 'RHOF', 'RHOG', 'RHOH', 'RHOJ', 
		'RHOQ', 'RHOT1', 'RHOT2', 'RHOU', 'RHOV', 'SOS1', 'SOS2', 'SRGAP1', 'SRGAP2', 'SRGAP3', 'STARD13', 'STARD8', 'SYDE1', 'SYDE2', 'TAGAP', 'TIAM1', 'TIAM2', 
		'TRIO', 'TRIP10', 'VAV1', 'VAV2', 'VAV3', 'ARHGEF28', 'ARHGEF25', 'ARHGEF8', 'ARHGEF29', 'DOCK2', 'DOCK3', 'DOCK180', 'DOCK10', 'ALS2', 'ARHGEF34P', 
		'DNMBP', 'ECT2L', 'FARP1', 'FARP2', 'FGD5', 'FGD6', 'ITSN2', 'MCF2', 'MCF2L2', 'PLEKHG1', 'PLEKHG3', 'PLEKHG4', 'PREX2', 'RASGRF1', 'SPATA13', 'ARHGAP16P', 
		'ARHGAP45', 'SH3BP1', 'RHOBTB1', 'RHOBTB2', 'RND1', 'RND2', 'RND3']
trio_peds = ['1010002', '1010004', '1010005', '1010007', '1030001', '1010008', '1010010', '1010013', '1010019', '1010020', '1010021', '1010022', '1010023', '1010024', 
		'1010025', '1010026', '1010028', '1030002', '1030006', '1060004', '1070001', '1070002', '1070003', '1070011', '1010031', '1010032', '1030013', '1030014', 
		'1030015', '1060015', '1060020', '1070012', '1070013', '1070015', '1030011', '1060023', '1060029', 'CFM-MOS-03', '1060018', '1060012', '4010002']
duo_peds = ['1010003', '1010006', '1030007']
infile_suffixes = ['.gatk.inherited.xls', '.freebayes.inherited.xls', '.intesect.inherited.xls', '.gatk.auto_dom.xls', '.freebayes.auto_dom.xls', '.intesect.auto_dom.xls', 
			'.freebayes_analysis.xls', '.gatk_analysis.xls', '.intesect_analysis.xls']
combined_name = 'combined.rho_genes'
gene_name_col = 5
##run gemini on all pedigrees
# for ped in trio_peds:
# 	modified_gemini_protocol(working_dir, ped, 'trio')
# for ped in duo_peds:
# 	modified_gemini_protocol(working_dir, ped, 'duo')

##combine results files and filter for rho genes
for infile_suffix in infile_suffixes:
		combine_files_and_filter_by_gene_list(combined_name, infile_suffix, rho_genes, gene_name_col)

