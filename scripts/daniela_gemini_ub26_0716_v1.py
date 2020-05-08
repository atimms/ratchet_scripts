#!/tools/BioBuilds-2015.04/bin/python
import os
import subprocess
import glob
import shutil

##set input variables and parameters
delim = '\t'
##working directory
working_dir = '/data/atimms/microtia_exomes/for_gemini'
# gemini_temp_dir = '/home/atimms/gemini_temp'
os.chdir(working_dir)

##programs and ref files
vt = '/data/atimms/gemini/vt/vt'
gemini = '/data/atimms/gemini/bin/gemini'
fasta = '/data/atimms/references/human_g1k_v37.fasta'
snpeff_jar = '/data/atimms/gemini/snpEff/snpEff.jar'
bcftools = '/home/atimms/programs/bcftools-1.3/bcftools'
bgzip = '/tools/BioBuilds-2015.04/bin/bgzip'
tabix = '/tools/BioBuilds-2015.04/bin/tabix'

##all pedigrees, all parents coded as unaffected
ped_dn_rec = 'microtia_exomes.dn_recessive.0616.ped'
prefix_dn_rec = 'microtia_exomes.dn_recessive.0616'
gemini_db_dn_rec = 'microtia_exomes.dn_recessive.0616.db'
##all peds, all unaffected coded as unknown
ped_dom = 'microtia_exomes.dominant.0616.ped'
prefix_dom = 'microtia_exomes.dominant.0616'
gemini_db_dom = 'microtia_exomes.dominant.0616.db'
##all peds, all carriers and ptags coded as unknown
ped_real = 'microtia_exomes.real.0616.ped'
prefix_real = 'microtia_exomes.real.0616'
gemini_db_real = 'microtia_exomes.real.0616.db'

query_suffix_dict = {'denovo_a': '.de_novo.c10.pass.not_low.aaf.xls', 'denovo_b': '.de_novo.c10.not_low.aaf.xls', 'denovo_c': '.de_novo.c10.pass.not_low.xls', 
		'denovo_d': '.de_novo.c20.pass.not_low.aaf.xls', 'denovo_e': '.de_novo.c10.pass.not_low.aaf.allow_unaff.xls', 'denovo_f': '.de_novo.not_low.xls', 
		'denovo_g': '.de_novo.c5.pass.not_low.aaf.xls','comp_het_a':'.comp_het.c10.pass.not_low.aaf.xls', 'auto_rec_a': '.auto_rec.c10.pass.not_low.aaf.xls', 
		'auto_rec_b': '.auto_rec.c10.not_low.aaf.xls', 'auto_rec_c': '.auto_rec.c10.pass.not_low.aaf10.xls', 'auto_rec_d': '.auto_rec.c20.pass.not_low.aaf.xls', 
		'auto_rec_e': '.auto_rec.c10.pass.not_low.aaf.allow_unaff.xls', 'auto_rec_f': '.auto_rec.not_low.aaf10.xls', 'auto_rec_g': '.auto_rec.c5.pass.not_low.aaf.xls', 
		'auto_rec_h': '.auto_rec.c10.pass.not_low.aaf5.xls','comp_het_b':'.comp_het.c10.pass.not_low.aaf5.xls', 'auto_dom_a':'.auto_dom.c10.pass.not_low.aaf.xls', 
		'auto_dom_b':'.auto_dom.c20.pass.not_low.aaf.xls', 'xlinked_a':'.x_linked.c10.pass.not_low.xls', 'pm_a':'.aff_ge0.35.parent_le0.35.c10.pass.not_low.aaf.xls', 
		'pm_b':'.aff_ge0.4.parent_le0.2.c10.pass.not_low.aaf.xls', 'pm_c':'.aff_ge0.35.parent_le0.35.c20.pass.not_low.aaf.xls', 'pm_d':'.aff_ge0.4.parent_le0.2.c20.pass.not_low.aaf.xls', 
		'xlinked_dn_a':'.xlinked_dn.c10.pass.not_low.aaf.xls', 'denovo_syn_a': '.de_novo.c10.pass.synonymous.g2.aaf1.xls', 'auto_rec_syn_a': '.auto_rec.c10.pass.synonymous.g2.aaf1.xls', 
		'auto_rec_syn_b': '.auto_rec.c10.pass.synonymous.g2.aaf5.xls', 'potential_cnv_a': '.potential_cnv.c10.pass.not_low.aaf1.xls', 'potential_cnv_b': '.potential_cnv.c10.pass.not_low.aaf5.xls', 
		'xlinked_rec_a': '.xlinked_rec.c10.pass.not_low.aaf1.xls', 'xlinked_rec_b': '.xlinked_rec.c10.pass.not_low.aaf5.xls', 'xlinked_dom_a': '.xlinked_dom.c10.pass.not_low.aaf1.xls',
		'xlinked_de_novo_a': '.xlinked_de_novo.c10.pass.not_low.aaf1.xls',  'auto_dom_c':'.auto_dom.c10.pass.not_low.not_lenient.aaf.xls'}

def prep_vcf_for_gemini(peds, prefix, ped_file, db_name):
	##collect all vcf file names
	all_vcfs = []
	for ped in peds:
		input_vcf = ped + '.intersected_vcfs/0002.vcf.gz'
		bcf_index = subprocess.Popen([bcftools, 'index', input_vcf])
		bcf_index.wait()
		all_vcfs.append(input_vcf)
	##combine vcfs
	input2_vcf = 'temp0.vcf.gz'
	bcf_merge = subprocess.Popen([bcftools, 'merge'] + all_vcfs + ['-O', 'z', '-o', input2_vcf, '-m', 'none'])
	bcf_merge.wait()
	##vcf names
	normalized_vcf = prefix + '.int.norm.vcf'
	temp_vcf = 'temp1.vcf'
	##dcompress, change a gatk header thing, and decompse and normalize (should be redundant)
	with open (temp_vcf, 'w') as tvcf_fh:
		zless_vcf = subprocess.Popen(['zless', input2_vcf], stdout=subprocess.PIPE)
		sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
		vt_decompose = subprocess.Popen([vt, 'decompose', '-s', '-'], stdin=sed_vcf.stdout, stdout=subprocess.PIPE)
		vt_normalize = subprocess.Popen([vt, 'normalize', '-r', fasta, '-'], stdin=vt_decompose.stdout, stdout=tvcf_fh)
		vt_normalize.wait()
	##annotate with snpeff and compress and index
	with open (normalized_vcf, 'w') as nvcf_fh:
		snpeff_vcf = subprocess.Popen(['java', '-Xmx20G', '-jar', snpeff_jar, 'GRCh37.75', '-v', '-formatEff', '-classic', temp_vcf], stdout=nvcf_fh)
		snpeff_vcf.wait()
		bgzip_vcf = subprocess.Popen([bgzip, normalized_vcf])
		bgzip_vcf.wait()
		tabix_vcf = subprocess.Popen([tabix, '-p', 'vcf', normalized_vcf + '.gz'])
		tabix_vcf.wait()
	#add file to gemini, requires ped file
	gemini_load = subprocess.Popen([gemini, 'load', '--cores', '10', '-t', 'snpEff', '-p', ped_file, '-v', normalized_vcf + '.gz', db_name])
	gemini_load.wait()
	##remove intermediate files
	os.remove(input2_vcf)
	os.remove(temp_vcf)

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


# gemini de_novo --filter "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01" -d 10 /data/atimms/gemini/dobyns_exome_gemini.db
def query_gemini(query_name, prefix, db_name):
	outfile = prefix + query_suffix_dict[query_name]
	with open (outfile, 'w') as out_fh:
		##de novo: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		if query_name == 'denovo_a':
			gemini_query = subprocess.Popen([gemini, 'de_novo', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '-d', '10', db_name], stdout=out_fh)
			gemini_query.wait()
		##de novo: impact med or high, aaf <=1% and coverage >=10 in all family members
		elif query_name == 'denovo_b':
			gemini_query = subprocess.Popen([gemini, 'de_novo', '--filter', "impact_severity != 'LOW' AND max_aaf_all <= 0.01", '-d', '10', db_name], stdout=out_fh)
			gemini_query.wait()
		##de novo: passed by gatk, impact med or high, and coverage >=10 in all family members
		elif query_name == 'denovo_c':
			gemini_query = subprocess.Popen([gemini, 'de_novo', '--filter', "filter IS NULL AND impact_severity != 'LOW'", '-d', '10', db_name], stdout=out_fh)
			gemini_query.wait()
		##de novo: passed by gatk, impact med or high, aaf <=1% and coverage >=20 in all family members
		elif query_name == 'denovo_d':
			gemini_query = subprocess.Popen([gemini, 'de_novo', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '-d', '20', db_name], stdout=out_fh)
			gemini_query.wait()
		##de novo: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members and allow variant to be in unaffected
		elif query_name == 'denovo_e':
			gemini_query = subprocess.Popen([gemini, 'de_novo', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '--allow-unaffected', '-d', '10', db_name], stdout=out_fh)
			gemini_query.wait()
		##de novo: passed by gatk, impact med or high
		elif query_name == 'denovo_f':
			gemini_query = subprocess.Popen([gemini, 'de_novo', '--filter', "impact_severity != 'LOW'", db_name], stdout=out_fh)
			gemini_query.wait()
		##de novo: passed by gatk, impact med or high, aaf <=1% and coverage >=5 in all family members
		elif query_name == 'denovo_g':
			gemini_query = subprocess.Popen([gemini, 'de_novo', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '-d', '5', db_name], stdout=out_fh)
			gemini_query.wait()
		##compound het: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		elif query_name == 'comp_het_a':
			gemini_query = subprocess.Popen([gemini, 'comp_hets', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '-d', '10', db_name], stdout=out_fh)
			gemini_query.wait()
		##compound het: passed by gatk, impact med or high, aaf <=5% and coverage >=10 in all family members
		elif query_name == 'comp_het_b':
			gemini_query = subprocess.Popen([gemini, 'comp_hets', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.05", '-d', '10', db_name], stdout=out_fh)
			gemini_query.wait()
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		elif query_name == 'auto_rec_a':
			gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '-d', '10', db_name], stdout=out_fh)
			gemini_query.wait()
		##auto recessive: impact med or high, aaf <=1% and coverage >=10 in all family members
		elif query_name == 'auto_rec_b':
			gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--filter', "impact_severity != 'LOW' AND max_aaf_all <= 0.01", '-d', '10', db_name], stdout=out_fh)
			gemini_query.wait()
		##auto recessive: passed by gatk, impact med or high, aaf <=10% and coverage >=10 in all family members
		elif query_name == 'auto_rec_c':
			gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.1", '-d', '10', db_name], stdout=out_fh)
			gemini_query.wait()
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=20 in all family members
		elif query_name == 'auto_rec_d':
			gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '-d', '20', db_name], stdout=out_fh)
			gemini_query.wait()
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members and allow variant to be in unaffected
		elif query_name == 'auto_rec_e':
			gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '--allow-unaffected', '-d', '10', db_name], stdout=out_fh)
			gemini_query.wait()
		##auto recessive: passed by gatk, impact med or high, aaf <=10% 
		elif query_name == 'auto_rec_f':
			gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--filter', "impact_severity != 'LOW' AND max_aaf_all <= 0.1", db_name], stdout=out_fh)
			gemini_query.wait()
		##dauto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=5 in all family members
		elif query_name == 'auto_rec_g':
			gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '-d', '5', db_name], stdout=out_fh)
			gemini_query.wait()
		##dauto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=5 in all family members
		elif query_name == 'auto_rec_h':
			gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.05", '-d', '10', db_name], stdout=out_fh)
			gemini_query.wait()
		##autosomal dominant, as precurser to parental mosaic, nit working yet (need one parent to be classified as affected??)
		elif query_name == 'auto_dom_a':
			gemini_query = subprocess.Popen([gemini, 'autosomal_dominant', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '--lenient', '--allow-unaffected', '-d', '10', db_name], stdout=out_fh)
			gemini_query.wait()
		elif query_name == 'auto_dom_b':
			gemini_query = subprocess.Popen([gemini, 'autosomal_dominant', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '--lenient', '--allow-unaffected', '-d', '20', db_name], stdout=out_fh)
			gemini_query.wait()
		elif query_name == 'auto_dom_c':
			gemini_query = subprocess.Popen([gemini, 'autosomal_dominant', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '--allow-unaffected', '-d', '10', db_name], stdout=out_fh)
			gemini_query.wait()

		##de novo synonymous: passed by gatk, synonymous, aaf <=1% and coverage >=10 in all family members
		elif query_name == 'denovo_syn_a':
			gemini_query = subprocess.Popen([gemini, 'de_novo', '--filter', "filter IS NULL AND impact_severity == 'LOW' AND max_aaf_all <= 0.01 and is_coding == 1 and gerp_bp_score >= 2", '-d', '10', db_name], stdout=out_fh)
			gemini_query.wait()
		##auto rec synonymous: passed by gatk, synonymous, aaf <=1% and coverage >=10 in all family members
		elif query_name == 'auto_rec_syn_a':
			gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--filter', "filter IS NULL AND impact_severity == 'LOW' AND max_aaf_all <= 0.01 and is_coding == 1 and gerp_bp_score >= 2", '-d', '10', db_name], stdout=out_fh)
			gemini_query.wait()
		##auto rec synonymous: passed by gatk, synonymous, aaf <=1% and coverage >=10 in all family members
		elif query_name == 'auto_rec_syn_b':
			gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--filter', "filter IS NULL AND impact_severity == 'LOW' AND max_aaf_all <= 0.05 and is_coding == 1 and gerp_bp_score >= 2", '-d', '10', db_name], stdout=out_fh)
			gemini_query.wait()
		##potential cnv vars
		elif query_name == 'potential_cnv_a':
			all_auto_rec_file = 'all_auto_rec.temp.xls'
			correct_auto_rec_file = 'correct_auto_rec.temp.xls'
			with open (correct_auto_rec_file, 'w') as car_fh:
				gemini_query_a = subprocess.Popen([gemini, 'autosomal_recessive', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '-d', '10', db_name], stdout=car_fh)
				gemini_query_a.wait()
			with open (all_auto_rec_file, 'w') as aar_fh:
				gemini_query_b = subprocess.Popen([gemini, 'autosomal_recessive', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '--lenient', '-d', '10', db_name], stdout=aar_fh)
				gemini_query_b.wait()
			correct_cnv_var_dict = make_var_dict_from_gemini_result(correct_auto_rec_file)
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
		elif query_name == 'potential_cnv_b':
			all_auto_rec_file = 'all_auto_rec.temp.xls'
			correct_auto_rec_file = 'correct_auto_rec.temp.xls'
			with open (correct_auto_rec_file, 'w') as car_fh:
				gemini_query_a = subprocess.Popen([gemini, 'autosomal_recessive', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.05", '-d', '10', db_name], stdout=car_fh)
				gemini_query_a.wait()
			with open (all_auto_rec_file, 'w') as aar_fh:
				gemini_query_b = subprocess.Popen([gemini, 'autosomal_recessive', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.05", '--lenient', '-d', '10', db_name], stdout=aar_fh)
				gemini_query_b.wait()
			correct_cnv_var_dict = make_var_dict_from_gemini_result(correct_auto_rec_file)
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

		##xlined recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		elif query_name == 'xlinked_rec_a':
			gemini_query = subprocess.Popen([gemini, 'x_linked_recessive', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '-d', '10', db_name], stdout=out_fh)
			gemini_query.wait()
		##xlined recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		elif query_name == 'xlinked_rec_b':
			gemini_query = subprocess.Popen([gemini, 'x_linked_recessive', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.05", '-d', '10', db_name], stdout=out_fh)
			gemini_query.wait()
		##xlined dominant: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		elif query_name == 'xlinked_dom_a':
			gemini_query = subprocess.Popen([gemini, 'x_linked_dominant', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '-d', '10', db_name], stdout=out_fh)
			gemini_query.wait()
		##xlined de novo: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		elif query_name == 'xlinked_de_novo_a':
			gemini_query = subprocess.Popen([gemini, 'x_linked_de_novo', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '-d', '10', db_name], stdout=out_fh)
			gemini_query.wait()
		else:
			print 'query name %s not recognized'%query_name

def add_new_ped_file_to_gemini_db(pedfile, geminidb):
	gemini_amend = subprocess.Popen([gemini, 'amend', '--sample', pedfile, geminidb])
	gemini_amend.wait()


def get_parental_mosaic_from_auto_dom(query_name, prefix, db_name, min_aaf_affected, max_aaf_unknown):
	infile = prefix + query_suffix_dict[query_name]
	outfile = prefix + '.aff_ge' + str(min_aaf_affected) + '.parent_le' + str(max_aaf_unknown) + '.' + query_suffix_dict[query_name].split(".", 2)[2]
	# infile = 'temp.xls'
	print infile, outfile
	with open(infile, 'r') as in_fh, open(outfile, 'w') as out_fh:
		line_count, final_count = 0, 0
		for line in in_fh:
			line_count += 1
			if line_count == 1:
				out_fh.write(line)
			else:
				aff_aafs, unknown_aafs = [], []
				line = line.split(delim)
				var_id = line[4]
				family_members = line[151].split(',')
				# print family_members
				for family_member in family_members:
					family_member = family_member.split('(')[1]
					ind_id = family_member.split(';')[0]
					aff_status = family_member.split(';')[1]
					gemini_query = subprocess.Popen([gemini, 'query', '-q', "gt_depths." + ind_id + ", gt_ref_depths." + ind_id + ", gt_alt_depths." + ind_id + "from variants where variant_id=" + var_id , db_name], stdout=subprocess.PIPE)
					gemini_query.wait()
					ind_covs, err = gemini_query.communicate()
					ind_covs = ind_covs.rstrip().split(delim)
					if ind_covs[0] == '0':
						aaf = 0.0
					else:
						aaf = int(ind_covs[2])/ float(ind_covs[0])
					# print ind_id, aff_status, var_id, ind_covs, aaf
					if aff_status == 'affected':
						aff_aafs.append(aaf)
					elif aff_status == 'unknown':
						unknown_aafs.append(aaf)
					else:
						print 'affection status %s not recognised'%aff_status
				# print aff_aafs, unknown_aafs
				##if  minimum aaf in affected >= req (consitutional) and max aff in unknown is less than required (mosaic), both parents are not hom ref i.e max >0.01 aaf, and both are not mosaic i.e. min < 0.01
				if min(aff_aafs) >= min_aaf_affected and max(unknown_aafs) <= max_aaf_unknown and max(unknown_aafs) >= 0.01 and min(unknown_aafs) < 0.01:
					print var_id, family_members, aff_aafs, unknown_aafs
					out_fh.write(delim.join(line))
					final_count += 1
	print '%s autosomal dom variants checked and %s determined to be parental mosaic'%(line_count -1, final_count)


def cp_gemini_db(in_db, out_db):
	cp_command = subprocess.Popen(['cp', in_db, out_db])
	cp_command.wait()




##for trios and multiple affects sibs 
##parameters
pedigrees = ['1010002', '1010003', '1010004', '1010005', '1010006', '1010007', '1010008', '1010010', '1010013', '1010019', 
		'1010020', '1010021', '1010022', '1010023', '1010024', '1010025', '1010026', '1010028', '1030001', '1030002', '1030006', 
		'1030007', '1060004', '1060015', '1070001', '1070002', '1070003', '1070011', '1010031', '1010032', '1030013', '1030014', 
		'1030015', '1060018', '1060020', '1070012', '1070013', '1070015']




##run methods
##add vcfs to db
# prep_vcf_for_gemini(pedigrees, prefix_dn_rec, ped_dn_rec, gemini_db_dn_rec)

##run analysis
##denovo - use denovo_a and  denovo_g i.e. 10 and 5x coverage
# query_gemini('denovo_a', prefix_dn_rec, gemini_db_dn_rec)
# query_gemini('denovo_g', prefix_dn_rec, gemini_db_dn_rec)
# ##autosomal recessive - use 5% aap ie auto_rec_h and 1% i.e. a also g 1% and 5x coverage
# query_gemini('auto_rec_h', prefix_dn_rec, gemini_db_dn_rec)
# query_gemini('auto_rec_a', prefix_dn_rec, gemini_db_dn_rec)
# query_gemini('auto_rec_g', prefix_dn_rec, gemini_db_dn_rec)
# ##compound het - use 5% aap i.e. comp_het_b and 1% ie a
# query_gemini('comp_het_b', prefix_dn_rec, gemini_db_dn_rec)
# query_gemini('comp_het_a', prefix_dn_rec, gemini_db_dn_rec)
# ##get conserved synonymous mutations: denovo and ar
# query_gemini('denovo_syn_a', prefix_dn_rec, gemini_db_dn_rec)
# query_gemini('auto_rec_syn_a', prefix_dn_rec, gemini_db_dn_rec)
# query_gemini('auto_rec_syn_b', prefix_dn_rec, gemini_db_dn_rec)

# ##potential cnvs
# query_gemini('potential_cnv_a', prefix_dn_rec, gemini_db_dn_rec)
# query_gemini('potential_cnv_b', prefix_dn_rec, gemini_db_dn_rec)

# ##new xlinked tests
# query_gemini('xlinked_rec_a', prefix_dn_rec, gemini_db_dn_rec)
# query_gemini('xlinked_rec_b', prefix_dn_rec, gemini_db_dn_rec)
# query_gemini('xlinked_de_novo_a', prefix_dn_rec, gemini_db_dn_rec)


##make alternative db with parents coded as unknown and run auto_dominant
# cp_gemini_db(gemini_db_dn_rec, gemini_db_dom)
# ##add ped file
# add_new_ped_file_to_gemini_db(ped_dom, gemini_db_dom)
# ##autosomal dominant use 10x coverage 
# query_gemini('auto_dom_a', prefix_dom, gemini_db_dom)

# ##make alternative db with parents coded with best guess 
# cp_gemini_db(gemini_db_dn_rec, gemini_db_real)
# ##add ped file
add_new_ped_file_to_gemini_db(ped_real, gemini_db_real)
# ##autosomal dominant use 10x coverage 
query_gemini('auto_dom_c', prefix_real, gemini_db_real)




##compare numbers etc locally using script 

