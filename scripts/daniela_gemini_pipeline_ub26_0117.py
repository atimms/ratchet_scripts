#!/tools/BioBuilds-2015.04/bin/python
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




def load_single_vcf_into_gemini(ped, db_name, input_vcf, ped_file):
	##dcompress, change a gatk header thing, and decompse and normalize (should be redundant)
	temp_vcf = 'temp1.vcf'
	normalized_vcf = ped + '.int.norm.vcf'
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


def gemini_comp_het(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.comp_hets.xls'
	with open (outfile, 'w') as out_fh:
		##compound het: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'comp_hets', '--max-priority', '2', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter',  "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf) , '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

def gemini_recessive(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.autosomal_recessive.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
# 		gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--columns', 'chrom,start,end,ref,alt,gene,biotype,impact,aa_change,gerp_bp_score,cadd_scaled,max_aaf_all', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '-d', '10', gem_db_name], stdout=out_fh)

		gemini_query.wait()
	return outfile

def gemini_denovo(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.de_novo.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'de_novo', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

##not used atm
def gemini_dominant(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.autosomal_dominant.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'autosomal_dominant', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

def gemini_inherited(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.inher.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'autosomal_dominant', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile


def gemini_xlinked(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.x_linked_recessive.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'x_linked_recessive', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

def gemini_xlinked_de_novo(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.x_linked_de_novo.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'x_linked_de_novo', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

def gemini_de_novo_syn(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.de_novo_syn.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'de_novo', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity == 'LOW' and is_coding == 1 and gerp_bp_score >= 2 AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

def gemini_recessive_syn(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.auto_rec_syn.xls'
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

def gemini_potential_cnv(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.potential_cnv.xls'
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
def standard_gemini_protocol(working_directory, pedigree, ped_type):
	os.chdir(working_directory)
	gemini_db = pedigree + '.gemini.db'
	in_vcf = pedigree + '.intersected_vcfs/0002.vcf.gz'
	freq_req = 0.01
	freq_req_recessive = 0.05
	coverage_req = 10
	##load vcf file into db
	load_single_vcf_into_gemini(pedigree, gemini_db, in_vcf, ped_dn_rec)
	##run queries and get files to combine
	files_to_combine = []
	if ped_type == 'singleton' or ped_type == 'duo':
		files_to_combine.append(gemini_comp_het(pedigree, gemini_db, freq_req_recessive, coverage_req))
		files_to_combine.append(gemini_recessive(pedigree, gemini_db, freq_req_recessive, coverage_req))
		# gemini_dominant(pedigree, gemini_db, freq_req, coverage_req)
	elif ped_type == 'trio':
		files_to_combine.append(gemini_comp_het(pedigree, gemini_db, freq_req_recessive, coverage_req))
		files_to_combine.append(gemini_recessive(pedigree, gemini_db, freq_req_recessive, coverage_req))
		files_to_combine.append(gemini_denovo(pedigree, gemini_db, freq_req, coverage_req))
		files_to_combine.append(gemini_xlinked(pedigree, gemini_db, freq_req_recessive, coverage_req))
		files_to_combine.append(gemini_xlinked_de_novo(pedigree, gemini_db, freq_req, coverage_req))
		files_to_combine.append(gemini_de_novo_syn(pedigree, gemini_db, freq_req, coverage_req))
		files_to_combine.append(gemini_recessive_syn(pedigree, gemini_db, freq_req_recessive, coverage_req))
		files_to_combine.append(gemini_potential_cnv(pedigree, gemini_db, freq_req_recessive, coverage_req))
		# gemini_dominant(pedigree, gemini_db, freq_req, coverage_req)
	else:
		print 'ped type %s not recognized'%ped_type
	##combine results
	combine_gemini_results(pedigree, files_to_combine, 5, 13, 22, 6, '.std_analysis.xls')
	##get auto dom vars regardless of family history
	add_new_ped_file_to_gemini_db(ped_dom, gemini_db)
	gemini_dominant(pedigree, gemini_db, freq_req, coverage_req)
	combine_gemini_results(pedigree, [pedigree + '.cov' + str(coverage_req) + '.maf'+ str(freq_req) + '.autosomal_dominant.xls'], 5, 13, 22, 6, '.auto_dom.xls')
	##get inherited vars - should be same as auto dom unless is family history
	add_new_ped_file_to_gemini_db(ped_real, gemini_db)
	gemini_inherited(pedigree, gemini_db, freq_req, coverage_req)
	combine_gemini_results(pedigree, [pedigree + '.cov' + str(coverage_req) + '.maf'+ str(freq_req) + '.inher.xls'], 5, 13, 22, 6, '.inherited.xls')


#call methods
##tests
# standard_gemini_protocol(working_dir, 'CFM-MOS-03', 'trio')
# standard_gemini_protocol(working_dir, '4010002', 'trio')
# standard_gemini_protocol(working_dir, '1060018', 'trio')


##batch1 -- ready
# standard_gemini_protocol(working_dir, '1010002', 'trio')
# standard_gemini_protocol(working_dir, '1010003', 'duo')
# standard_gemini_protocol(working_dir, '1010004', 'trio')
# standard_gemini_protocol(working_dir, '1010005', 'trio')
# standard_gemini_protocol(working_dir, '1010006', 'duo')
# standard_gemini_protocol(working_dir, '1010007', 'trio')
# standard_gemini_protocol(working_dir, '1030001', 'trio')

##batch2 -- ready
# standard_gemini_protocol(working_dir, '1010008', 'trio')
# standard_gemini_protocol(working_dir, '1010010', 'trio')
# standard_gemini_protocol(working_dir, '1010013', 'trio')
# standard_gemini_protocol(working_dir, '1010019', 'trio')
# standard_gemini_protocol(working_dir, '1010020', 'trio')
# standard_gemini_protocol(working_dir, '1010021', 'trio')
# standard_gemini_protocol(working_dir, '1010022', 'trio')
# standard_gemini_protocol(working_dir, '1010023', 'trio')
# standard_gemini_protocol(working_dir, '1010024', 'trio')
# standard_gemini_protocol(working_dir, '1010025', 'trio')
# standard_gemini_protocol(working_dir, '1010026', 'trio')
# standard_gemini_protocol(working_dir, '1010028', 'trio')

##batch3 -- ready
# standard_gemini_protocol(working_dir, '1030002', 'trio')
# standard_gemini_protocol(working_dir, '1030006', 'trio')
# standard_gemini_protocol(working_dir, '1030007', 'duo')
# standard_gemini_protocol(working_dir, '1060004', 'trio')
# standard_gemini_protocol(working_dir, '1070001', 'trio')
# standard_gemini_protocol(working_dir, '1070002', 'trio')
# standard_gemini_protocol(working_dir, '1070003', 'trio')
# standard_gemini_protocol(working_dir, '1070011', 'trio')

##batch4 -- ready
# standard_gemini_protocol(working_dir, '1010031', 'trio')
# standard_gemini_protocol(working_dir, '1010032', 'trio')
# standard_gemini_protocol(working_dir, '1030013', 'trio')
# standard_gemini_protocol(working_dir, '1030014', 'trio')
# standard_gemini_protocol(working_dir, '1030015', 'trio')
# standard_gemini_protocol(working_dir, '1060015', 'trio')
# standard_gemini_protocol(working_dir, '1060020', 'trio')
# standard_gemini_protocol(working_dir, '1070012', 'trio')
# standard_gemini_protocol(working_dir, '1070013', 'trio')
# standard_gemini_protocol(working_dir, '1070015', 'trio')


##batch5 -- ready
# standard_gemini_protocol(working_dir, '1030011', 'trio')
# standard_gemini_protocol(working_dir, '1060023', 'trio')
# standard_gemini_protocol(working_dir, '1060029', 'trio')
# standard_gemini_protocol(working_dir, '4010002','trio' )
# standard_gemini_protocol(working_dir, 'CFM-MOS-03', 'trio')
# standard_gemini_protocol(working_dir, '1060018', 'trio')
# standard_gemini_protocol(working_dir, '1060012', 'trio')



def combine_std_variant_type_and_filter_by_gene_list(outfile_prefix, infile_suffix, analysis_name, analysis_types, genelist_name, genelist, analysis_col, genename_col):
	outfile = outfile_prefix + analysis_name + '.' + genelist_name + '.xls'
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
						analysis_type = line[analysis_col]
						genename = line[genename_col]
						# print analysis_type
						if analysis_type in analysis_types and genename in genelist:
							out_fh.write(delim.join(line + ['\n']))

def combine_files_and_filter_by_gene_list(outfile_prefix, infile_suffix, genelist_name, genelist, genename_col):
	analysis_name = infile_suffix.split('.')[1]
	outfile = outfile_prefix + analysis_name + '.' + genelist_name + '.xls'
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

def combine_var_results_files(outfile_prefix, infile_suffix):
	outfile = outfile_prefix + infile_suffix.split('.')[1] + '.txt'
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
						out_fh.write(delim.join(line + ['\n']))




##combine different variants and compare against gene lists

##for genets genes
combined_name = 'combined.'
known_microtia_genes = ['GMNN', 'FGF10', 'FGF3', 'FGFR3', 'CDC6', 'FGFR2', 'CDT1', 'POLR1C', 'POLR1A', 'POLR1D', 'FAT4', 'TCOF1', 'RPS28', 'DHODH', 'RPS26', 'DCHS1', 'GRIP1', 'FRAS1', 'EYA1', 'SIX1', 'TSR2', 'SF3B4', 'HOXD', 'CHD7', 'GLI3', 'SALL4', 'EFTUD2', 'KMT2D', 'CDC45', 'HOXA2', 'SIX5', 'EDNRA', 'ORC1', 'PHF9', 'PLCB4', 'ORC6', 'HOXA1', 'FREM2', 'CDC18L', 'TWIST2', 'MGORS6', 'POMT1', 'HMX1', 'GDF6', 'ORC4', 'SALL1', 'SEMA3E', 'GNAI3', 'TFAP2A']
genets_genes = ['MCM6', 'HSPA8', 'NCK1', 'CD3EAP', 'FRS2', 'ORC4L', 'FASN', 'PNN', 'WEE1', 'GSTK1', 'MNAT1', 'HOXB2', 'FRS3', 'FGF23', 'STAT1', 'ORC2L', 'RRN3', 'LIX1L', 'CCNH', 'EEF1A1', 'GNA14', 'FGF5', 'FGF4', 'UBTF', 'RPLP1', 'EPHA4', 'GNA11', 'NDRG1', 'MYH9', 'EIF3B', 'EIF1AX', 'FGF13', 'POLR1E', 'FGF18', 'FGF6', 'ZYX', 'FGF19', 'TTF1', 'FGF17', 'RPL21', 'FGF22', 'FGF9', 'HTRA1', 'POLR1B', 'FGF8', 'ZNRD1', 'POLR3G', 'ZFPM1', 'EYA3', 'MCC', 'ORC5L', 'ARRB2', 'PIK3R2', 'GTF2H2', 'ATP5F1', 'PTRF', 'CBL', 'POLR3A', 'ERBB3', 'IKBKG', 'ORC6L', 'FGF7', 'GTF2H3', 'MAN2C1', 'TAF1A', 'EIF3E', 'POLR3B', 'TAF1B', 'GNA15', 'ANAPC1', 'LMNA', 'PHLDA3', 'TCF3', 'DFFA', 'RBM10', 'GNAQ', 'RGS12', 'RPL26L1', 'WBP4', 'CDC14A', 'PYGL', 'IKBKB', 'NFKB1', 'BARX1', 'CANX', 'FGFR1', 'ORC1L', 'RPL26', 'TRAP1', 'ASL', 'CUL1', 'EIF3G', 'DDX6', 'CDK7', 'EDN3', 'SIX3', 'CNOT1', 'SMARCD1', 'ERCC3', 'RBL2', 'RPL9P7', 'EIF4A2', 'POLR2G', 'ERCC2', 'RPS5', 'ETF1', 'MCM8', 'RPL9', 'PPFIA1', 'POLR2K', 'RPL24', 'TAF1C', 'NFKB2', 'PLXND1', 'RPL19', 'EIF3J', 'EIF5B', 'CSTF1', 'SOS1', 'ORC3L', 'CEBPD', 'HCFC1', 'SNORA7A', 'MAP3K14', 'RPL37A', 'HMGB1', 'EIF3I', 'CDH2', 'CDC2L5', 'FGF20', 'DGKZ', 'GTF2H1', 'FGF1', 'SALL3', 'CD59', 'RERE', 'YSK4', 'RPL30', 'EIF3H', 'SMAD3', 'LSM2', 'VARS2', 'RPL15', 'RAB1A', 'POLR3H', 'EDN1', 'WDR18', 'COPG', 'POLR2E', 'ARF1', 'TTN', 'RPL14', 'RASAL2', 'PDGFB', 'GRM6', 'CLIC1', 'MKI67', 'PPA1', 'RBM8A', 'RPL8', 'TUBGCP2', 'NR3C1', 'RPS8', 'RPS29', 'RPS21', 'DACH1', 'REL', 'FLG2', 'IGF2BP1', 'RPS25', 'POLR3F', 'DDX1', 'RPL7', 'DST', 'NOLC1', 'CEP70', 'EIF4B', 'STK19', 'GRB2', 'GOLM1', 'UBE1C', 'SNRPA', 'NRG1', 'LBR', 'WDR77', 'RPL36A', 'TERF1', 'NSUN2', 'PTPRS', 'PABPC3', 'CKAP5', 'POLR3K', 'SMU1', 'UBE2N', 'EP300', 'ICK', 'PCBP2', 'HNRPH2', 'BRCA1', 'SFRS3', 'RPL35', 'CENPF', 'HNRNPR', 'EPN2', 'CFL2', 'ITGB1', 'ARHGEF6', 'POLR3D', 'NCBP1', 'SRRM1', 'ACTR8', 'FGF14', 'BARD1', 'PKLR', 'EIF3D', 'FGF21', 'ATP6V1B2', 'CUTL1', 'SKI', 'CDKN2A', 'GTF2F2', 'FLII', 'WWP2', 'RPL13', 'FGF11', 'BTRC', 'CDC7', 'TBL2', 'EIF2S1', 'MYC', 'EIF3A', 'IGF2BP3', 'RFC4', 'POLR3C', 'HNRNPA0', 'SUMO1P3', 'RPL38', 'RPS14', 'RPS20', 'PEPD', 'MOBKL1B', 'SNRPD1', 'CBS', 'FOXA2', 'ZNF207', 'LIPT2', 'CDK2', 'RPL35A', 'FGF12', 'ITGA5', 'YWHAH', 'HIST1H1C', 'CDKN1A', 'FLNA', 'DDX23', 'EIF3F', 'RPL31', 'NFKBIB', 'MAGEA11', 'ARRB1', 'CKS2', 'UBE4B', 'ASNA1', 'RPL34', 'CREBBP', 'DCP1A', 'PPIF', 'CHUK', 'RPL12', 'LRP11', 'PIK3C3', 'RPL3L', 'NDUFS1', 'ME1', 'RPL41', 'ATP6V1A', 'RPL18', 'MOV10', 'PLS3', 'RBX1', 'KRTAP4-12', 'PDGFRA', 'HTATIP2', 'UBE3C', 'MDFI', 'ARPC3', 'HNRNPA1', 'GNB4', 'NCBP2', 'SNRPA1', 'INS', 'DGKQ', 'RPL10A', 'RPL23', 'RPS4X', 'RPS4Y1', 'SRRM2', 'DCP2', 'NUDT21', 'PTP4A3', 'RPS26P10', 'SFRS1', 'NPM1', 'CAPN1', 'EGR2', 'HIST2H2BE', 'C12orf4', 'CIR1', 'TMPO', 'GSPT2', 'MCM10', 'RPL11', 'NHP2L1', 'DKC1', 'DHX9', 'PCF11', 'IRF2BP1', 'HDAC7A', 'PHGDH', 'RPS3A', 'GTF2F1', 'RPL14P1', 'LAMC1', 'RBM5', 'SYBL1', 'RPS15', 'DLX4', 'SUPT4H1', 'RPL28', 'SPRY2', 'TBCA', 'PABPC1', 'RAB5A', 'CYP4F2', 'RBM39', 'GNA13', 'HNRPUL1', 'CCAR1', 'NONO', 'VASP', 'PRDM16', 'SETX', 'NDUFA10', 'RPL37', 'RPL4', 'TTC4', 'CBLB', 'FGFR1OP2', 'SMAD2', 'PTPRF', 'ACTA1', 'SFRS11', 'EIF1B', 'EIF2B4', 'RPS13', 'NFX1', 'FGFR4', 'JUN', 'EFNB1', 'NOLA1', 'SMNDC1', 'POLR2I', 'RPLP2', 'IPO13', 'KYNU', 'EEF1A2', 'DNM1', 'STK4', 'MAPK8', 'PRKACG', 'HIP1', 'EIF4H', 'CPSF2', 'SFRS2', 'TXNL4A', 'BUD13', 'EIF3K', 'RPL29', 'CTU2', 'KRT8', 'RPS9', 'TIMM13', 'RNPS1', 'MCM4', 'TLX3', 'RPS16', 'BCLAF1', 'SF3B1', 'RPL22', 'DYNC1H1', 'FYN', 'RPA4', 'RPLP0', 'RPS23', 'RPL5', 'POLR2J', 'GNAZ', 'MCM3', 'EPHA7', 'AKAP8', 'EIF2B2', 'EIF2S3', 'EIF2B1', 'RPS26P3', 'BMPR1A', 'EGFR', 'FEN1', 'MPRIP', 'SMC2', 'GRAMD1A', 'SFRS7', 'SARS', 'TEX10', 'RPS6KA2', 'CHMP4B', 'EIF3C', 'GCDH', 'EPHB2', 'DDX3X', 'FBXW11', 'QKI', 'P4HB', 'RPL7A', 'CUL3', 'PPP2R3A', 'SSBP1', 'PRKY', 'UBE2G2', 'IQGAP1', 'PRPF4', 'CAMKK1', 'RPL10', 'RPS26P6', 'RBMX2', 'PDGFRB', 'PPP2R2A', 'HNRPM', 'PDIA6', 'RBPJ', 'IDH3G', 'LSM1', 'ETFA', 'EIF4E', 'HNRPD', 'RASA1', 'FGF2', 'IGBP1', 'OAT', 'EIF2B3', 'PABPN1', 'IKZF5', 'DMD', 'CPSF7', 'RPS19', 'CORO1B', 'RELA', 'PPM1A', 'RPSAP58', 'POLR2D', 'DDX20', 'LRP4', 'SNRPE', 'CCNA2', 'TBCE', 'POLR2H', 'CSTF3', 'RPS6', 'MLLT4', 'PRKAG2', 'RPS7', 'SRR', 'TCP1', 'C6orf211', 'KIF26A', 'RFC5', 'ITPR1', 'LRP2', 'EIF2S2', 'CAPZA2', 'LAMA5', 'UPF3B', 'POLR2L', 'CSTF2', 'PCBP1', 'HNRNPC', 'DEK', 'PPFIA3', 'U2AF1', 'RFX2', 'UBE2R2', 'MRPL27', 'CNOT2', 'CLP1', 'SCN1A', 'COMMD10', 'PC', 'NOL5A', 'WWP1', 'MAP1B', 'SFRS12', 'MYO1C', 'FAM62C', 'CERK', 'EIF3CL', 'STAT3', 'PDCD2', 'CEBPZ', 'SNRPD2', 'RBM17', 'NAF1', 'LIMA1', 'SFRS9', 'FUS', 'UFD1L', 'AGTR2', 'SIP1', 'CCDC5', 'POLR2F', 'SNRPF', 'SNORA67', 'RPS27A', 'DDX24', 'USP15', 'AP3D1', 'RPS12', 'GNAO1', 'SEC22A', 'CAMKK2', 'IARS', 'THRAP3', 'CIAO1', 'VPS28', 'KLF9', 'BMPR2', 'LEO1', 'DBN1', 'NHLH1', 'RPS10', 'HSF1', 'ARR3', 'SOS2', 'TCERG1', 'CAD', 'MCF2L', 'CTNNB1', 'ATP5A1', 'EDG1', 'MAF1', 'RPL36', 'ATIC', 'CDC40', 'PFDN6', 'TUBA4A', 'HADHB', 'ASCC3L1', 'MYO1B', 'EIF2S2P4', 'KIAA1377', 'PRKCA', 'TPM2', 'NIPSNAP1', 'U2AF2', 'PPM1B', 'SC5DL', 'RSL1D1', 'YWHAZ', 'SF3B2', 'POLR2A', 'HNRNPL', 'ADORA1', 'GRPEL1', 'TBP', 'SNRPB2', 'CPSF1', 'KRI1', 'GCN1L1', 'SORBS1', 'RPS24', 'TPM1', 'TOP1', 'MDH2', 'VCP', 'SH3KBP1', 'IPO4', 'TOX', 'SERPINH1', 'SMAD5', 'ATP5H', 'TUBB2A', 'CHERP', 'PLSCR1', 'RTF1', 'HSPA1L', 'DDX19B', 'G3BP1', 'WDR92', 'XRCC5', 'PRMT1', 'CD48', 'ARF4', 'FBXW8', 'NDUFS2', 'RPL23A', 'APC', 'TRAF6', 'RPL27A', 'GRIP2', 'TPTEP1', 'EP400', 'SLC25A3', 'ELAVL1', 'PRLHR', 'MCM2', 'SF4', 'C4orf27', 'RPN2', 'CSNK1D', 'MAGEB2', 'GTF3C4', 'MAGOH', 'ATXN1', 'WBP11', 'ARFGEF2', 'GPR88', 'SMN2', 'PRKACA', 'PCGF2', 'DDX31', 'FAU', 'OPN4', 'PPP2R1A', 'PICALM', 'GRIK2', 'SCFD1', 'DNAJC8', 'DNAJB6', 'HIBCH', 'WDR57', 'RBMX', 'GRIN1', 'EIF2AK4', 'NFASC', 'CCT6A', 'PPP1R12A', 'SF3A1', 'TXN', 'GSTM3', 'RCC2', 'HOXD10', 'HTR1A', 'XPO1', 'MDN1', 'RRAGB', 'DHX38', 'DAPK3', 'EIF4A3', 'HMGB2', 'PPP2R5C', 'TRIM28', 'GNG3', 'CHEK1', 'FECH', 'MAPKAP1', 'CITED4', 'HNRPF', 'GNB2L1', 'ACTA2', 'IDH3A', 'SNRPN', 'MYO9A', 'EIF2A', 'TP53', 'AHSA1', 'NSFL1C', 'RPS17', 'ACOT8', 'SHOC2', 'CUL2', 'EGF', 'DDX17', 'FXR1', 'CLSPN', 'RANBP9', 'NOLA3', 'CYP3A5', 'RPL6', 'ERH', 'ATP5C1', 'RPL27AP', 'PRPF6', 'HIST1H2AE', 'PSMC5', 'WDR74', 'TOX3', 'ACTC1', 'PHC1', 'PLCG1', 'HMG1L1', 'PLCE1', 'CAPZB', 'CYR61', 'YWHAG', 'NTHL1', 'MAP3K7', 'CCR5', 'KRR1', 'CTNNA1', 'DIAPH1', 'UTX', 'PAPOLA', 'CCNB2', 'SMC4', 'H1FX', 'RGS2', 'IQGAP2', 'CTTN', 'ACTB', 'TUBB2C', 'DICER1', 'YWHAB', 'PPP6C', 'UTP3', 'YBX1', 'MYL6', 'CUL5', 'IQGAP3', 'LYAR', 'YWHAQ', 'MEIS1', 'LONP1', 'DDX55', 'GMPPB', 'ACACB', 'COPA', 'CCNE1', 'PCSK5', 'WIPF2', 'SCMH1', 'RPL39', 'DDX27', 'GARS', 'PRKAB1', 'RPL18A', 'ATP6V1E1', 'PPIA', 'TUBB2B', 'SOD1', 'HNRPA3', 'IGHMBP2', 'ARPC5', 'EIF4G1', 'DDX3Y', 'METTL7A', 'PDCD11', 'GNAI2', 'NDUFA8', 'RAC1', 'SDHA', 'ATP2A2', 'ULK3', 'CALD1', 'LYN', 'DNM1L', 'SPIN1', 'ASS1', 'MSH3', 'RFC3', 'OPRD1', 'APLP1', 'GSC', 'ARHGAP32', 'RPS2', 'VIM', 'NOS1AP', 'GABPA', 'ARHGAP17', 'COPB2', 'KRT18', 'SFRS5', 'NUP133', 'PSPC1', 'UBA52', 'CCDC94', 'USO1', 'ROS1', 'ATAD3A', 'HIST1H2AI', 'ZDHHC17', 'HDLBP', 'HNRNPA2B1', 'PGD', 'GNAI1', 'VDAC3', 'RARS2', 'RAF1', 'SURF6', 'SPECC1', 'ISCA2', 'VPS13A', 'SAPS3', 'AARS', 'RIOK2', 'ABLIM1', 'UQCC', 'RAP1A', 'BMPR1B', 'GDI1', 'WDR26', 'UBE2M', 'C1orf25', 'ROCK1', 'H2AFX', 'SAG', 'ACTR2', 'RPP30', 'YKT6', 'PARP1', 'PIAS1', 'PCSK7', 'IL8RB', 'C9orf114', 'SS18L1', 'SDC3', 'CHMP6', 'SLMAP', 'TRIP6', 'TMOD3', 'EWSR1', 'SUMO3', 'HNRPK', 'PIK3C2A', 'TRPV4', 'LAS1L', 'RPL3', 'DGKG', 'MAPK10', 'DCD', 'MYBBP1A', 'PAPD5', 'MOGS', 'SAMM50', 'MED6', 'RC3H1', 'GRAMD1B', 'PELP1', 'STX12', 'MCM7', 'SFRS6', 'PTBP1', 'HIST3H2A', 'CALM2', 'PICK1', 'EFNB3', 'SDC2', 'AP2M1', 'SLC25A1', 'XPR1', 'MED18', 'PTGES2', 'AKAP8L', 'SF3B3', 'ROBO2', 'TRRAP', 'POLR2J3', 'HIST1H2BO', 'BUB3', 'ARHGEF11', 'PSMC6', 'P4HA1', 'EIF2B5', 'HIST1H2BM', 'EEF2', 'DDX52', 'ATP5B', 'HIST1H1A', 'RPS3', 'MRPS15', 'C15orf15', 'ATF3', 'NDUFA5', 'MYO9B', 'PRPF4B', 'MOCS3', 'ANP32B', 'PSMA7', 'DNAJA2', 'FH', 'HELLS', 'DDB1', 'MAPK9', 'MED9', 'OGT', 'TRIP10', 'LMO7', 'WRN', 'KCTD1', 'AP2A1', 'CDK4', 'RPS6KA3', 'NCOA1', 'CDC2', 'GDI2', 'PTPRU', 'CPNE1', 'BXDC5', 'HIST1H3B', 'WASF2', 'SRPK1', 'SFRS15', 'SNRP70', 'SLC25A5', 'PPP1R9B', 'TOP2A', 'RGS14', 'ENAH', 'CCT4', 'HIST1H2AG', 'RPL27', 'XRCC3', 'DLSTP', 'SMAD1', 'UBC', 'SYNPO', 'POLR2C', 'PKMYT1', 'SLC5A5', 'ATP5O', 'BTF3L4', 'SARM1', 'HIGD1A', 'HNRNPU', 'HIST1H3J', 'RPL13AP25', 'BDP1', 'RBPMS', 'ADCY2', 'MRPS34', 'IMMT', 'POLR2B', 'RSF1', 'PPP1CB', 'CCDC55', 'CAMK1D', 'CUL4B', 'BAT3', 'CSNK1A1L', 'HIST1H3D', 'CRKRS', 'TTLL3', 'RAB11FIP5', 'ATR', 'GRN', 'DIAPH2', 'CAMK2D', 'AP2B1', 'PDCD6', 'HIST1H2AM', 'MCRS1', 'PPIAP22', 'PTCD3', 'PRMT5', 'DDX42', 'MON2', 'ADCY5', 'MRPS16', 'RPL15P3', 'CTPS', 'PPP2R1B', 'ZNF121', 'DDX50', 'CPA1', 'SNRPD3', 'SDAD1', 'IDH3B', 'PCID2', 'PRPF8', 'H3F3AP6', 'DLD', 'COX6B1', 'RNF130', 'SNIP1', 'NDUFS6', 'TUBA1B', 'GNL3', 'WDR59', 'HIST1H2AD', 'HUWE1', 'ADCY1', 'KIF2C', 'RANGAP1', 'SNRPG', 'RGS5', 'SUPT5H', 'ZNF501', 'EDG2', 'YES1', 'CTSK', 'ANKRD25', 'CHIC2', 'SNRPB', 'ABCF3', 'KLC1', 'PRKX', 'CLTB', 'SYT1', 'ZIC3', 'PLEC1', 'TOP3A', 'NCKAP1L', 'RANBP5', 'CCT8', 'UTP11L', 'ARPC1A', 'VAV3', 'SLC25A4', 'GRWD1', 'STK3', 'MMS19', 'CITED2', 'GTF3C2', 'RBM34', 'SMARCC2', 'XBP1', 'SAE1', 'AHCTF1P', 'SFRS4', 'PHF5A', 'CHD1', 'GNL1', 'DHX15', 'RPS18', 'RNH1', 'EFNB2', 'AOF1', 'MED16', 'RADIL', 'SUFU', 'TH1L', 'SACM1L', 'HIST1H3C', 'RRP12', 'SFRS17A', 'MED25', 'OLA1', 'PCSK6', 'GTF3C3', 'PLA1A', 'PEX5', 'RGS7', 'MYL12A', 'RACGAP1', 'SRSF10', 'ILF3', 'CUL4A', 'MCM5', 'CWF19L1', 'DYNLL1', 'HIST1H3A', 'PSMB2', 'ZMAT2', 'SF3A3', 'ANAPC2', 'PIAS2', 'MAGI1', 'FXR2', 'HOXB7', 'RGS3', 'RIPK3', 'MAP3K1', 'LRP1B', 'POU2F1', 'PWP1', 'MKI67IP', 'UGDH', 'MTIF2', 'RAB10', 'IRS4', 'PLK1', 'NPLOC4', 'TUBB', 'NOTCH1', 'DLG4', 'HOXC8', 'SOBP', 'SF3B5', 'EXOSC10', 'SPDYA', 'SRF', 'UBE2H', 'HIST1H2AK', 'LIAS', 'ACTN4', 'G3BP2', 'GSK3B', 'S100A9', 'MED31', 'LARS', 'ESR1', 'FLNB', 'NUMBL', 'KPNA5', 'USP24', 'TCF7L2', 'GTF3C1', 'CCDC17', 'PSPHL', 'ABCB7', 'CAMK2G', 'KLHDC3', 'PTPN11', 'GNAS', 'CALM1', 'MYL12B', 'FBL', 'USP52', 'MSH2', 'ISCA1', 'HSP90AB1', 'CPNE8', 'NAV2', 'C20orf191', 'CYP3A43', 'CCDC97', 'CALM3', 'PSMD6', 'EBNA1BP2', 'FAM115A', 'UPF1', 'HIST1H3E', 'HCLS1', 'DOM3Z', 'PRKAA2', 'HIST2H2BA', 'IMPDH1', 'SH3D19', 'NETO2', 'MAPRE1', 'RBM4', 'COX5A', 'PLOD3', 'TBL1XR1', 'ALKBH8', 'LSM4', 'TSR1', 'ATP6V0D1', 'HERC5', 'MATR3', 'FANCD2', 'PELO', 'PAH', 'ADCK1', 'SLC25A13', 'HSPA6', 'APPBP1', 'AP2S1', 'BCKDHB', 'RPA2', 'RUVBL1', 'DGKE', 'MED13', 'HIST1H2BD', 'RPS15A', 'UBR5', 'SRP9', 'GBF1', 'NANOGP8', 'BOP1', 'SUCLG1', 'RRP9', 'C1QBP', 'PTK2B', 'PFKP', 'SCYL2', 'HSPA9', 'NOC2L', 'NAP1L1', 'FKBP1A', 'MRPS23', 'SUMO4', 'PHC1P1', 'RNGTT', 'SHC1', 'HIST1H2BJ', 'NGB', 'ZIC1', 'RGS10', 'NADSYN1', 'SART3', 'SRP14', 'CAMK2A', 'ILF2', 'SMC1A', 'GUK1', 'NVL', 'GUSBL1', 'MAP3K7IP1', 'DRD2', 'LY6G5B', 'IKBKAP', 'MRPS27', 'RNF149', 'GPSM2', 'HIST1H2BL', 'DPH1', 'KPNA2', 'UNC93A', 'MRM1', 'ZNF622', 'XPO5', 'BMS1', 'PES1', 'SMARCB1', 'LPXN', 'RPN1', 'CSE1L', 'MYO5C', 'ACTG2', 'HMGCR', 'MYH10', 'KARS', 'LSM6', 'PDLIM7', 'MED30', 'PMM2', 'GUSBP9', 'AASS', 'SAFB2', 'SMEK1', 'PAK7', 'SLC25A11', 'CSPG4', 'NALCN', 'FSCN1', 'GRK6', 'RB1', 'VPS24', 'HIST1H4H', 'C20orf160', 'NSA2', 'NANOG', 'IRX2', 'CAP2', 'PSMB9', 'MYH1', 'PSMD11', 'KIAA0372', 'CCT3', 'ABL1', 'ACTN1', 'MED15', 'GRM3', 'IMPDH2', 'SPARCL1', 'SMG7', 'HIST1H4L', 'DYNLL2', 'WDR12', 'PIH1D1', 'SYNE1', 'RBM9', 'GTPBP4', 'HSPA1B', 'PPP2R3B', 'EIF2C2', 'OR10B1P', 'ACTR10', 'LRPPRC', 'CAPRIN1', 'KHSRP', 'SFPQ', 'NIP7', 'DDX56', 'TIAF1', 'DHX30', 'HIST1H4F', 'DSCC1', 'P2RY13', 'TPM4', 'CETN3', 'AP2A2', 'LARP7', 'NOL10', 'NUMA1', 'LMNB1', 'CSNK2B', 'SNX13', 'GNL2', 'PPP3CC', 'SHARPIN', 'TTLL9', 'SHANK1', 'C5AR1', 'MED11', 'NEDD8', 'DCAKD', 'PCNA', 'HSPH1', 'MARS', 'MRTO4', 'AFF4', 'BUD31', 'CASP6', 'EVL', 'TIA1', 'HIST1H4D', 'FOXA1', 'CPT1B', 'NLE1', 'MORF4L1', 'VPS25', 'EDN2', 'PPP2CB', 'PCMT1', 'ARHGAP15', 'RPL13A', 'EMG1', 'TARDBP', 'PTCH2', 'POMT2', 'TUBB4', 'SMARCAD1', 'DHCR7', 'SKP2', 'STOML2', 'HAP1', 'HIST1H4K', 'SMARCA2', 'INA', 'KRAS', 'AP3B2', 'PPME1', 'HIST3H2BB', 'RNASEH2B', 'PKM2', 'MED8', 'PPP1CC', 'LARP2', 'PRKACB', 'TUBA3D', 'BAG2', 'DDX43', 'SERBP1', 'MED21', 'SPTAN1', 'ANXA2P3', 'HBZ', 'IDE', 'PFAS', 'RQCD1', 'SIN3A', 'CCNC', 'THOC4', 'DHX36', 'FMNL1', 'RDX', 'MIOS', 'CCDC85B', 'HIST1H4I', 'HSPB1', 'FAM55B', 'ZNF593', 'PPP1CA', 'RRP1B', 'GUF1', 'MED4', 'ACSL4', 'XPOT', 'ERCC4', 'INCENP', 'SKIV2L', 'GSTP1', 'GSG2', 'HIST1H4B', 'PSMB3', 'VAV2', 'MED27', 'MED19', 'HIST1H4E', 'UBE2I', 'KIAA0664', 'CEPT1', 'PIK3R4', 'RPS27', 'SCIN', 'MYOCD', 'HIST1H4J', 'C1orf156', 'CFL1', 'NEXN', 'CLTA', 'TAF5', 'CD2BP2', 'NEK9', 'DDX54', 'DUSP4', 'KRT17', 'APBB1', 'HIST1H2AH', 'H2BFS', 'U2SURP', 'MAP3K3', 'PPP1R10', 'FTSJ3', 'YWHAE', 'GRIK3', 'KIAA1967', 'DDOST', 'TFAP2D', 'ZNF326', 'UBB', 'RPS11', 'PRDX1', 'TPT1', 'TUBA1A', 'NDUFS3', 'SSRP1', 'ARPC5L', 'TUBA3C', 'FZR1', 'ACACA', 'C9orf140', 'PFKFB3', 'MTNR1B', 'NME2', 'POR', 'SHMT2', 'GCC1', 'SMC6', 'GEMIN4', 'GSN', 'EIF6', 'AHCY', 'NAT10', 'BXDC1', 'PPP2CA', 'BRAF', 'ATP1A1', 'NCOA3', 'HIST2H4A', 'DDX18', 'CKAP4', 'GBA2', 'ARD1A', 'ANAPC5', 'BYSL', 'TAOK2', 'HDAC2', 'HIST1H3H', 'HIST1H2BH', 'CORO1C', 'GLI2', 'RGS4', 'CNGA3', 'PGRMC1', 'SPARC', 'DPY30', 'PSMD14', 'PPM1G', 'PSMA6', 'NCK2', 'GMFB', 'PNO1', 'MED14', 'KCTD3', 'SPECC1L', 'TUBB1', 'PGM1', 'CSNK1A1', 'POLE', 'RPA1', 'PSMD4', 'SEH1L', 'VRK1', 'SSR4', 'HIST2H2BF', 'C17orf79', 'NEDD4L', 'NUP153', 'POLD1', 'DPM1', 'GRIA3', 'GTF2I', 'STK38L', 'TLE1', 'CACNG7', 'ATF2', 'RBBP4', 'ETV7', 'PPIG', 'BACH2', 'SLC25A6', 'TUFM', 'HEATR1', 'EIF4ENIF1', 'DHX29', 'PLAA', 'PSMA5', 'EIF5', 'MYO5B', 'PPP3CB', 'PROX1', 'VPS35', 'PRDM1', 'MED20', 'AMPH', 'PHB2', 'MED24', 'TMEM33', 'KIAA0368', 'HARS', 'MYO6', 'EPB41L3', 'KHDRBS1', 'HNRPH1', 'GRIPAP1', 'PAWR', 'ACO2', 'RGS19', 'HDAC6', 'S100A6', 'ACTBL2', 'ATAD3B', 'CDC2L6', 'GADD45G', 'PRKD2', 'ARHGEF2', 'RGS16', 'CBR1', 'SPTBN1', 'TUBA1C', 'POLA2', 'HTATIP', 'HIST1H4A', 'SLC47A1', 'PSKH1', 'ARPC2', 'UACA', 'C6orf157', 'CDK3', 'UTP20', 'IDH1', 'SMAP1', 'HMX3', 'BRSK1', 'CSNK1E', 'MYO3B', 'POU2F2', 'TWISTNB', 'XRN1', 'PLRG1', 'CLPB', 'DUSP21', 'AES', 'HSPA5', 'RPAP2', 'SNURF', 'HOXA7', 'XRCC6', 'DHX16', 'SEC61A1', 'CCT5', 'IRAK1', 'PTPRD', 'RFC2', 'DNM2', 'SLC7A9', 'CCDC75', 'SPTLC1', 'MYH14', 'FBXL13', 'MYLK', 'LHX8', 'LDLRAP1', 'ZRANB2', 'LUZP1', 'AR', 'DDX5', 'SMARCA4', 'VPS52', 'CTNND1', 'CSNK2A2', 'GNA12', 'AGGF1', 'MED23', 'TRAPPC6A', 'POTEKP', 'CDC45L', 'HIST1H2BK', 'MRPL15', 'HOXA11', 'DBF4', 'PARD6G', 'SUPT16H', 'AMN1', 'SPIN2B', 'SAR1B', 'FBXW7', 'ZNF650', 'UBR2', 'OSBPL5', 'CDC5L', 'DNAH3', 'CASP3', 'P4HA2', 'HSP90AA2', 'DDX49', 'CHTF18', 'CDC123', 'WDR68', 'PDE8A', 'EYA4', 'WAS', 'SH2B1', 'INOC1', 'GRIA2', 'HIST2H4B', 'TUBA8', 'SMAD4', 'TBCK', 'IGF2BP2', 'PA2G4', 'PPP3R2', 'SART1', 'ARHGEF4', 'GRIK1', 'MVD', 'GANAB', 'PTPN14', 'CACNA1D', 'PAK1IP1', 'MSH6', 'HOXA10', 'CSNK2A1P', 'DUSP14', 'SNORD23', 'XAB1', 'ABLIM3', 'CHD4', 'IPO9', 'CNOT4', 'SMARCC1', 'PHB', 'NUCB1', 'PSME4', 'CPSF3', 'DNAJA3', 'PBRM1', 'MED26', 'MED12', 'CHD9', 'GSTM1', 'EGR1', 'ARL1', 'UTP15', 'PDE4D', 'KIAA0020', 'ARIH2', 'VTI1A', 'PSMA4', 'GART', 'GNB1', 'HIST1H2BB', 'TUBG1', 'CLTC', 'ZNF281', 'DHX8', 'CASP1', 'PSMD7', 'SORT1', 'VAPB', 'PRPF3', 'SF1', 'NOL8', 'RPS27L', 'NOC4L', 'SLC4A5', 'EEF1G', 'DDX46', 'TUBB3', 'RCL1', 'C11orf41', 'SPIN3', 'PSMD8', 'DYM', 'CALR', 'HIPK3', 'TTC27', 'CNR1', 'CACNG5', 'RBBP7', 'HIST2H2AB', 'DDX39', 'EFEMP2', 'ASF1B', 'PPT1', 'GNB2', 'PHKA2', 'ANXA6', 'VPS4B', 'SYNCRIP', 'MED29', 'PAK4', 'DDX21', 'SLC25A12', 'PRPF19', 'SF3A2', 'BXDC2', 'ALDH1B1', 'EPRS', 'ITPR3', 'C9orf41', 'PSMC2', 'MRPL16', 'TTL', 'BMP2K', 'IMP3', 'HIST1H2AC', 'SLIT2', 'CSN2', 'TLE4', 'DIS3', 'PDPK1', 'GRIA4', 'PABPC4', 'CHN1', 'RAD54L2', 'SQRDL', 'FBXO44', 'KIF3A', 'COG4', 'LUC7L2', 'STXBP5', 'MAP2K2', 'MAT1A', 'SMARCA5', 'SQSTM1', 'ZIC2', 'NOL6', 'HIST4H4', 'ZBTB43', 'HTATSF1', 'ZBTB9', 'MYO1A', 'LRRC59', 'FOXA3', 'STK36', 'PPP3CA', 'ALDH18A1', 'MAX', 'SEC31A', 'MRPL17', 'NOL14', 'FREQ', 'PRDX2', 'MEX3D', 'ADORA3', 'WDR46', 'GRAP', 'C20orf4', 'SUPT6H', 'CCT2', 'ACTR3', 'CCNA1', 'KIAA1033', 'RICTOR', 'HRNR', 'HDAC1', 'SPAST', 'CPSF6', 'RPA3', 'HIST1H2BC', 'SHROOM3', 'TMEM16K', 'SRPK2', 'MED22', 'H3F3B', 'NCOR1', 'CDC23', 'ERO1LB', 'PATL1', 'GAS7', 'AURKA', 'SKP1A', 'ITM2C', 'TCP11L2', 'PSMD1', 'C13orf34', 'BSN', 'DPH2', 'BOLA2B', 'CWC22', 'KPNA4', 'LIPA', 'DGKH', 'CLNS1A', 'DHX37', 'C6orf47', 'HSPD1', 'C21orf77', 'CAMK2B', 'KPNA6', 'FNTB', 'LRRFIP2', 'PSMD3', 'PSMB5', 'SORCS2', 'C11orf84', 'PSMB6', 'MED10', 'IHPK1', 'MAEA', 'SMEK2', 'GPSM3', 'DEPDC5', 'PWP2', 'HIST1H2BE', 'EYA2', 'CNP', 'SEC23A', 'KCNAB1', 'LUC7L', 'TAF1', 'PI4KB', 'RNMT', 'SNW1', 'PNPLA7', 'CSDA', 'ACTR3B', 'MRPS22', 'EIF2C1', 'LUC7L3', 'YOD1', 'HOXC9', 'TREH', 'STIP1', 'PRPF40A', 'DUT', 'GMPS', 'PQBP1', 'CIRH1A', 'KTN1', 'UTP18', 'RANBP2', 'GFPT2', 'IPO7', 'SNRPEL1', 'PRKCB1', 'ZBTB7A', 'MED17', 'TBL3', 'CAB39', 'RRM2', 'HIST2H2AA4', 'XYLB', 'STXBP2', 'XAB2', 'MED13L', 'RAI14', 'SLC36A1', 'NACAP1', 'PRKCG', 'CYFIP1', 'SH3BP5L', 'ACIN1', 'NOL1', 'PI4KA', 'HIST1H2BG', 'TAF7L', 'KRT7', 'EIF4G3', 'RBM19', 'IMP4', 'MED28', 'HIST1H1D', 'ATP5J2', 'PSMD10', 'CENTD2', 'VAC14', 'MTNR1A', 'PPID', 'SLC13A4', 'GRIA1', 'WDSOF1', 'PSMC3', 'BRD2', 'LSM5', 'PSMC1', 'RGS18', 'TARS', 'WBSCR22', 'DLAT', 'NEK2', 'RASA4', 'BCR', 'SKIV2L2', 'PSMD12', 'EPS15', 'BTAF1', 'CAPZA1', 'PRKDC', 'SIX2', 'HLTF', 'USP14', 'ABCE1', 'RUVBL2', 'RBM28', 'HIST1H2BI', 'CRNKL1', 'ZCCHC10', 'CAP1', 'RASA4B', 'ARHGEF12', 'WDR1', 'ARHGAP21', 'ARFIP1', 'MED1', 'GAPDH', 'R3HDM1', 'EIF5AP1', 'MED7', 'NKTR', 'DHRS2', 'PSMB4', 'SVIL', 'PHLPPL', 'HRB', 'DNAJA1', 'ACTL6A', 'EEF1AL3', 'HAT1', 'DDX47', 'ESF1', 'C3orf26', 'SHROOM4', 'DHDH', 'TWIST1', 'SAPS1', 'PSMD13', 'LSM7', 'GRTP1', 'STK38', 'SLC13A1', 'WDR48', 'PPP5C', 'WDR36', 'PIP5K1A', 'RAD50', 'PSMA1', 'WDR3', 'DDX10', 'PSMB1', 'PSMA3', 'PIM2', 'KPNB1', 'DIMT1L', 'PRPF31', 'MTHFD1', 'AATF', 'TNPO1', 'ACTG1', 'NOB1', 'NET1', 'BRD3', 'NR2F1', 'ELL2', 'MRPS5', 'CCDC85C', 'FUSIP1', 'RIC8A', 'SMARCD3', 'PSMD2', 'MPHOSPH10', 'IRF3', 'BAT1', 'DNAJB1', 'PSMC4', 'NCL', 'FMR1', 'RAB43', 'PPAN-P2RY11', 'HIST1H2BF', 'PPP2R2C', 'RPSAP15', 'CHGB', 'ATG9A', 'NOP58', 'TSHR', 'SH3YL1', 'PSMA2', 'ADCY3', 'KIFC1', 'CDK8', 'CHKB', 'RGL3', 'ADCY6', 'MAPK1', 'TFIP11', 'DNALI1']
all_genes = known_microtia_genes + genets_genes
analysis_type_col = 41
gene_name_col = 5
std_analysis_types_dict = {'denovo':['de_novo', 'x_linked_de_novo'], 'recessive':['autosomal_recessive', 'comp_hets'], 'xlinked': ['x_linked_recessive', 'x_linked_de_novo']}
# gene_lists_dict = {'known_microtia_genes': known_microtia_genes, 'genets_genes': genets_genes, 'known_and_genets_genes': all_genes}
gene_lists_dict = {'known_microtia_genes': known_microtia_genes, 'known_and_genets_genes': all_genes}
infile_suffixes = ['.inherited.xls', '.auto_dom.xls']
##for the variants in the std analysis files
# for analysis_type in std_analysis_types_dict:
# 	for genelist in gene_lists_dict:
# 		combine_std_variant_type_and_filter_by_gene_list(combined_name, '.std_analysis.xls', analysis_type, std_analysis_types_dict[analysis_type], genelist, gene_lists_dict[genelist], analysis_type_col, gene_name_col)
##for auto dom and inherited
# for infile_suffix in infile_suffixes:
# 	for genelist in gene_lists_dict:
# 		combine_files_and_filter_by_gene_list(combined_name, infile_suffix, genelist, gene_lists_dict[genelist], gene_name_col)


##for ARHGEF6 related
# infile_suffixes = ['.inherited.xls', '.auto_dom.xls', '.std_analysis.xls']
# genelist = ['PARVB', 'GIT1', 'GIT2', 'PAK1', 'PAK2', 'PAK3', 'ARHGEF7', 'CDC42', 'RAC1', 'RHOA', 'ARHGEF6']
# for infile_suffix in infile_suffixes:
# 	combine_files_and_filter_by_gene_list(combined_name, infile_suffix, 'ARHGEF6_related', genelist, gene_name_col)

##combine all files regardless of analysis or gene
infile_suffixes = ['.inherited.xls', '.auto_dom.xls', '.std_analysis.xls']
for infile_suffix in infile_suffixes:
	combine_var_results_files(combined_name, infile_suffix)












