#!/usr/bin/env python
import os
import subprocess
import glob
import shutil


##set input variables and parameters
delim = '\t'
working_dir = '/data/atimms/microtia_exomes/daniela_add_annovar_0118'

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
av_protocol = ['-protocol', 'refGene,kaviar_20150923,gnomad_exome,dbnsfp33a,revel']
av_operation = ['-operation', 'g,f,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.']
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
ann_var = '/home/atimms/programs/annovar/annotate_variation.pl'

##methods

##download clinvar db
def download_annovar(file_to_download):
	annovar_dl = subprocess.Popen([ann_var, '-buildver', av_genome, '-downdb', file_to_download, av_ref_dir[0], '-webfrom', 'annovar'])
	annovar_dl.wait()

def combine_vcf_files(gz_vcfs, out_vcf):
	for vcf in gz_vcfs:
		# bgzip_vcf = subprocess.Popen(['bgzip', vcf])
		# bgzip_vcf.wait()
		bcf_index = subprocess.Popen(['bcftools', 'index', vcf])
		bcf_index.wait()
		# gz_vcfs.append(vcf + '.gz')
		# gz_vcfs.append(vcf)
	print gz_vcfs
	bcf_merge = subprocess.Popen(['bcftools', 'merge'] + gz_vcfs + ['-O', 'v', '-o', out_vcf, '-m', 'none'])
	bcf_merge.wait()

def add_annovar_info(normalized_vcf, in_file, out_file, pos_to_insert):
	##files
	region_file = 'temp.regions'
	temp_vcf = 'temp0.vcf.gz'
	av_file = 'temp.avinput'
	multianno = 'temp.hg19_multianno.txt'
	'''
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
	bgzip_vcf = subprocess.Popen(['bgzip', normalized_vcf])
	bgzip_vcf.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', normalized_vcf + '.gz'])
	bcf_index.wait()
	## get_var_from_regions_in_vcf(in_vcf, regions, out_vcf):
	bcftools_filter = subprocess.Popen([bcftools_12, 'view', '-R', region_file, '-o', temp_vcf, '-O', 'z', normalized_vcf + '.gz'])
	bcftools_filter.wait()
	##convert vcf file to individual avinput file
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', temp_vcf, '-allsample', '--withfreq', '-outfile', av_file])
	con_ann.wait()
	##run_table_annovar(avinput, av_prefix):
	command = [table_annovar] + av_buildver + [av_file] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', 'temp']
	annovar = subprocess.Popen(command)
	annovar.wait()
	'''
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
			av_info = line[10:89]
			# print chrom, pos, av_info
			ann_dict[ch_pos] = av_info

	with open(in_file, 'r') as in_fh, open(out_file, 'w') as out_fh:
		line_count = 0
		for line2 in in_fh:
			line2 = line2.split(delim)
			line_count += 1
			if line_count == 1:
				out_fh.write(delim.join(line2[:pos_to_insert] + ['Kaviar_AF', 'Kaviar_AC', 'Kaviar_AN', 'gnomAD_exome_ALL', 'gnomAD_exome_AFR', 'gnomAD_exome_AMR', 'gnomAD_exome_ASJ', 'gnomAD_exome_EAS', 'gnomAD_exome_FIN', 'gnomAD_exome_NFE', 'gnomAD_exome_OTH', 'gnomAD_exome_SAS', 'SIFT_score', 'SIFT_converted_rankscore', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_rankscore', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score', 'Polyphen2_HVAR_rankscore', 'Polyphen2_HVAR_pred', 'LRT_score', 'LRT_converted_rankscore', 'LRT_pred', 'MutationTaster_score', 'MutationTaster_converted_rankscore', 'MutationTaster_pred', 'MutationAssessor_score', 'MutationAssessor_score_rankscore', 'MutationAssessor_pred', 'FATHMM_score', 'FATHMM_converted_rankscore', 'FATHMM_pred', 'PROVEAN_score', 'PROVEAN_converted_rankscore', 'PROVEAN_pred', 'VEST3_score', 'VEST3_rankscore', 'MetaSVM_score', 'MetaSVM_rankscore', 'MetaSVM_pred', 'MetaLR_score', 'MetaLR_rankscore', 'MetaLR_pred', 'M-CAP_score', 'M-CAP_rankscore', 'M-CAP_pred', 'CADD_raw', 'CADD_raw_rankscore', 'CADD_phred', 'DANN_score', 'DANN_rankscore', 'fathmm-MKL_coding_score', 'fathmm-MKL_coding_rankscore', 'fathmm-MKL_coding_pred', 'Eigen_coding_or_noncoding', 'Eigen-raw', 'Eigen-PC-raw', 'GenoCanyon_score', 'GenoCanyon_score_rankscore', 'integrated_fitCons_score', 'integrated_fitCons_score_rankscore', 'integrated_confidence_value', 'GERP++_RS', 'GERP++_RS_rankscore', 'phyloP100way_vertebrate', 'phyloP100way_vertebrate_rankscore', 'phyloP20way_mammalian', 'phyloP20way_mammalian_rankscore', 'phastCons100way_vertebrate', 'phastCons100way_vertebrate_rankscore', 'phastCons20way_mammalian', 'phastCons20way_mammalian_rankscore', 'SiPhy_29way_logOdds', 'SiPhy_29way_logOdds_rankscore', 'Interpro_domain', 'GTEx_V6_gene', 'GTEx_V6_tissue', 'REVEL'] + line2[pos_to_insert:]))
			else:
				chrom2 = line2[0][3:]
				pos2 = line2[2]
				ch_pos2 = '_'.join([chrom2,pos2])
				if ch_pos2 in ann_dict:
					# print line2[:5],ch_pos2, ann_dict[ch_pos2]
					out_fh.write(delim.join(line2[:pos_to_insert] + ann_dict[ch_pos2] + line2[pos_to_insert:]))
				else:
					print 'not found:', line2[:5], ch_pos2
##run methods

##get extra dbs for annovar
# for db_file in ['gnomad_exome', 'dbnsfp33a', 'revel']:
# 	download_annovar(db_file)

##parameters
# intersected_vcfs = glob.glob('/data/atimms/microtia_exomes/batch1-5_0117/*intesect.int.norm.vcf.gz')
intersected_vcfs = glob.glob('/data/atimms/microtia_exomes/batch1-5_0117/*.intersected_vcfs/0002.vcf.gz')
combined_vcf = 'all_samples.intersect.0118.vcf'
danielas_infile = 'daniela_combined_var_file.txt'
danielas_outfile = 'daniela_combined_var_file.info_added.xls'
##combine individual files
# print intersected_vcfs, combined_vcf
# combine_vcf_files(intersected_vcfs, combined_vcf)

##annotate the vars and combine into final file
add_annovar_info(combined_vcf, danielas_infile, danielas_outfile, 33)





