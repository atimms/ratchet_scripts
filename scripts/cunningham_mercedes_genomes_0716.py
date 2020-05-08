#!/tools/BioBuilds-2015.04/bin/python
import os
import subprocess
import glob
import shutil

##set input variables and parameters
delim = '\t'
##working directory
working_dir = '/data/atimms/mercedes_genomes'
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
prefix_dn_rec = 'mercedes_genomes.0716'
ped_dn_rec = prefix_dn_rec + '.ped'
gemini_db_dn_rec = prefix_dn_rec + '.db'
vcf_files = ['cunningham_uwcmg_blls_m_1.HF.final.vcf.gz']
##paramters for queries
freq_req = 0.01
freq_req_recessive = 0.01
coverage_req = 10

def load_single_vcf_into_gemini(vcfs, prefix, ped_file, db_name):
	##dcompress, change a gatk header thing, and decompse and normalize (should be redundant)
	# input_vcf = vcfs[0]
	temp_vcf = 'temp1.vcf'
	normalized_vcf = prefix + '.int.norm.vcf'
	# with open (temp_vcf, 'w') as tvcf_fh:
	# 	zless_vcf = subprocess.Popen(['zless', input_vcf], stdout=subprocess.PIPE)
	# 	sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
	# 	vt_decompose = subprocess.Popen([vt, 'decompose', '-s', '-'], stdin=sed_vcf.stdout, stdout=subprocess.PIPE)
	# 	vt_normalize = subprocess.Popen([vt, 'normalize', '-r', fasta, '-'], stdin=vt_decompose.stdout, stdout=tvcf_fh)
	# 	vt_normalize.wait()
	# ##annotate with snpeff and compress and index
	# with open (normalized_vcf, 'w') as nvcf_fh:
	# 	snpeff_vcf = subprocess.Popen(['java', '-Xmx20G', '-jar', snpeff_jar, 'GRCh37.75', '-v', '-formatEff', '-classic', temp_vcf], stdout=nvcf_fh)
	# 	snpeff_vcf.wait()
	# 	bgzip_vcf = subprocess.Popen(['bgzip', normalized_vcf])
	# 	bgzip_vcf.wait()
	# 	tabix_vcf = subprocess.Popen(['tabix', '-p', 'vcf', normalized_vcf + '.gz'])
	# 	tabix_vcf.wait()
	##add file to gemini, requires ped file
	# gemini_load = subprocess.Popen([gemini, 'load', '--cores', '10', '-t', 'snpEff', '-p', ped_file, '-v', normalized_vcf + '.gz', db_name])
	##if need to change tmp directory i.e. it's filling up /tmp
	gemini_load = subprocess.Popen([gemini, 'load', '--cores', '18', '-t', 'snpEff', '-p', ped_file, '-v', normalized_vcf + '.gz', '--tempdir', working_dir, db_name])
	gemini_load.wait()
	os.remove(temp_vcf)

def gemini_comp_het(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.comp_hets.xls'
	with open (outfile, 'w') as out_fh:
		##compound het: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'comp_hets', '--max-priority', '2', '--columns', 'chrom,start,end,ref,alt,gene,biotype,impact,aa_change,gerp_bp_score,cadd_scaled,max_aaf_all', '--filter',  "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf) , '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

def gemini_recessive(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.autosomal_recessive.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--columns', 'chrom,start,end,ref,alt,gene,biotype,impact,aa_change,gerp_bp_score,cadd_scaled,max_aaf_all', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

def gemini_recessive_inc_noncoding(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.autosomal_recessive.inc_noncoding.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, aaf <=1% and coverage >=10 in all family members
		##need which cols i.e. exonic and encodde etc
		# gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--columns', 'chrom,start,end,ref,alt,gene,biotype,impact,aa_change,gerp_bp_score,cadd_scaled,max_aaf_all', '--filter', "filter IS NULL AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--filter', "filter IS NULL AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)

		gemini_query.wait()
	return outfile

def gemini_denovo(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.de_novo.xls'
	with open (outfile, 'w') as out_fh:
		##de novo: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'de_novo', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,biotype,impact,aa_change,gerp_bp_score,cadd_scaled,max_aaf_all', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

def gemini_denovo_inc_noncoding(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.de_novo.inc_noncoding.xls'
	with open (outfile, 'w') as out_fh:
		##de novo: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		# gemini_query = subprocess.Popen([gemini, 'de_novo', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,biotype,impact,aa_change,gerp_bp_score,cadd_scaled,max_aaf_all', '--filter', "filter IS NULL AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query = subprocess.Popen([gemini, 'de_novo', '--filter', "filter IS NULL AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile


def combine_gemini_results(ped, file_list, gene_col, position_to_insert, ped_lists):
	##file names
	outfile = ped + '.std_analysis.xls'
	temp_file = 'combined.temp.xls'
	print 'combining files:', file_list
	##combine all files, making sure col lengths are all the same
	with open(temp_file, "w") as temp_fh:
		file_count = 0
		for filename in file_list:
			analysis_type = filename.split('.')[-2]
			file_count += 1
			# print analysis_type, file_count
			##add header from first file
			with open(filename, "r") as in_fh:
				line_count = 0
				for line in in_fh:
					line_count +=1
					line = line.rstrip().split(delim)
					if line_count == 1:
						if file_count == 1:
							header = delim.join(line[:12] + line[13:17] + ['analysis', '\n'])
							temp_fh.write(header)
					else:
						if analysis_type[0] == 'x':
							# print analysis_type
							samples = ','.join(ped_lists[0])
							line_out = delim.join(line + [ped, 'na', 'na', samples, analysis_type, '\n'])
							temp_fh.write(line_out)
						else:
							line_out = delim.join(line[:12] + line[13:17] + [analysis_type, '\n'])
							temp_fh.write(line_out)





##run methods
# ##add vcfs to db
# load_single_vcf_into_gemini(vcf_files, prefix_dn_rec, ped_dn_rec, gemini_db_dn_rec)

##run queries and get files to combine
# files_to_combine = []
# files_to_combine.append(gemini_comp_het(prefix_dn_rec, gemini_db_dn_rec, freq_req_recessive, coverage_req))
# files_to_combine.append(gemini_recessive(prefix_dn_rec, gemini_db_dn_rec, freq_req_recessive, coverage_req))
# files_to_combine.append(gemini_denovo(prefix_dn_rec, gemini_db_dn_rec, freq_req, coverage_req))
# print files_to_combine
##combine results
# combine_gemini_results(pedigree, files_to_combine, 5, 12, x_linked_list)

##get non coding
gemini_denovo_inc_noncoding(prefix_dn_rec, gemini_db_dn_rec, freq_req, coverage_req)
gemini_recessive_inc_noncoding(prefix_dn_rec, gemini_db_dn_rec, freq_req_recessive, coverage_req)
##compare numbers etc locally using script 

