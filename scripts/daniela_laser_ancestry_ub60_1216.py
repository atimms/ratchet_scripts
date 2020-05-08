#!/usr/bin/env python
import os
import subprocess


##parameters
delim = '\t'
working_dir = '/data/atimms/microtia_laser_1216'
os.chdir(working_dir)

##programs and program files
laser = '/home/atimms/programs/LASER-2.03/laser'
pileup2seq = '/home/atimms/programs/LASER-2.03/pileup2seq/pileup2seq.py'
hgdp_site = '/data/atimms/microtia_laser_1216/HGDP/HGDP_938.site'
hgdp_bed = 'HGDP_938.bed'
hgdp_geno = '/data/atimms/microtia_laser_1216/HGDP/HGDP_938.geno'
hgdp_coord = '/data/atimms/microtia_laser_1216/HGDP/HGDP_938.RefPC.coord'
fasta = '/data/atimms/references/human_g1k_v37.fasta'
b1_bams = ['101000201.bwa_gatk.bam', '101000202.bwa_gatk.bam', '101000203.bwa_gatk.bam', '101000301.bwa_gatk.bam', '101000302.bwa_gatk.bam', '101000401.bwa_gatk.bam', '101000402.bwa_gatk.bam', '101000403.bwa_gatk.bam', '101000501.bwa_gatk.bam', '101000502.bwa_gatk.bam', '101000503.bwa_gatk.bam', '101000601.bwa_gatk.bam', '101000602.bwa_gatk.bam', '101000701.bwa_gatk.bam', '101000702.bwa_gatk.bam', '101000703.bwa_gatk.bam', '103000101.bwa_gatk.bam', '103000102.bwa_gatk.bam', '103000103.bwa_gatk.bam']
b2_bams = ['101000801.bwa_gatk.bam', '101000802.bwa_gatk.bam', '101000803.bwa_gatk.bam', '101001001.bwa_gatk.bam', '101001002.bwa_gatk.bam', '101001003.bwa_gatk.bam', '101001301.bwa_gatk.bam', '101001302.bwa_gatk.bam', '101001303.bwa_gatk.bam', '101001901.bwa_gatk.bam', '101001902.bwa_gatk.bam', '101001903.bwa_gatk.bam', '101002001.bwa_gatk.bam', '101002002.bwa_gatk.bam', '101002003.bwa_gatk.bam', '101002101.bwa_gatk.bam', '101002102.bwa_gatk.bam', '101002103.bwa_gatk.bam', '101002201.bwa_gatk.bam', '101002202.bwa_gatk.bam', '101002203.bwa_gatk.bam', '101002301.bwa_gatk.bam', '101002302.bwa_gatk.bam', '101002303.bwa_gatk.bam', '101002401.bwa_gatk.bam', '101002402.bwa_gatk.bam', '101002403.bwa_gatk.bam', '101002501.bwa_gatk.bam', '101002502.bwa_gatk.bam', '101002503.bwa_gatk.bam', '101002601.bwa_gatk.bam', '101002602.bwa_gatk.bam', '101002603.bwa_gatk.bam', '101002801.bwa_gatk.bam', '101002802.bwa_gatk.bam', '101002803.bwa_gatk.bam']
b3_bams = ['103000201.bwa_gatk.bam', '103000202.bwa_gatk.bam', '103000203.bwa_gatk.bam', '103000601.bwa_gatk.bam', '103000602.bwa_gatk.bam', '103000603.bwa_gatk.bam', '103000701.bwa_gatk.bam', '103000702.bwa_gatk.bam', '106000401.bwa_gatk.bam', '106000402.bwa_gatk.bam', '106000403.bwa_gatk.bam', '106001501.bwa_gatk.bam', '106001502.bwa_gatk.bam', '106001503.bwa_gatk.bam', '107000101.bwa_gatk.bam', '107000102.bwa_gatk.bam', '107000103.bwa_gatk.bam', '107000201.bwa_gatk.bam', '107000202.bwa_gatk.bam', '107000203.bwa_gatk.bam', '107000301.bwa_gatk.bam', '107000303.bwa_gatk.bam', '107001101.bwa_gatk.bam', '107001102.bwa_gatk.bam', '107001103.bwa_gatk.bam']
b4_bams = ['101003101.bwa_gatk.bam', '101003102.bwa_gatk.bam', '101003103.bwa_gatk.bam', '101003201.bwa_gatk.bam', '101003202.bwa_gatk.bam', '101003203.bwa_gatk.bam', '103001301.bwa_gatk.bam', '103001302.bwa_gatk.bam', '103001303.bwa_gatk.bam', '103001401.bwa_gatk.bam', '103001402.bwa_gatk.bam', '103001403.bwa_gatk.bam', '103001501.bwa_gatk.bam', '103001502.bwa_gatk.bam', '103001503.bwa_gatk.bam', '106001801.bwa_gatk.bam', '106001806a.bwa_gatk.bam', '106001806b.bwa_gatk.bam', '106001809.bwa_gatk.bam', '106002001.bwa_gatk.bam', '106002002.bwa_gatk.bam', '106002003.bwa_gatk.bam', '107000302.bwa_gatk.bam', '107001201.bwa_gatk.bam', '107001202.bwa_gatk.bam', '107001203.bwa_gatk.bam', '107001301.bwa_gatk.bam', '107001302.bwa_gatk.bam', '107001303.bwa_gatk.bam', '107001501.bwa_gatk.bam', '107001502.bwa_gatk.bam', '107001503.bwa_gatk.bam']
all_bams = b1_bams + b2_bams + b3_bams + b4_bams
b1_prefix = 'batch1'
b1_hgdp_all_for_r = b1_prefix + '.hgdp_all.txt'
b1_hgdp_select_for_r = b1_prefix + '.hgdp_select.txt'
b1_hgdp_america_for_r = b1_prefix + '.hgdp_america.txt'
b14_prefix = 'batch1_4'
b14_hgdp_all_for_r = b14_prefix + '.hgdp_all.txt'
b14_hgdp_select_for_r = b14_prefix + '.hgdp_select.txt'
b14_hgdp_america_for_r = b14_prefix + '.hgdp_america.txt'

hgdp_continent_dict = {'Adygei' : 'Europe', 'Balochi' : 'Central and South Asia', 'BantuKenya' : 'Africa', 
		'BantuSouthAfrica' : 'Africa', 'Basque' : 'Europe', 'Bedouin' : 'Middle East', 
		'BiakaPygmy' : 'Africa', 'Brahui' : 'Central and South Asia', 
		'Burusho' : 'Central and South Asia', 'Cambodian' : 'East Asia', 'Colombian' : 'America', 
		'Dai' : 'East Asia', 'Daur' : 'East Asia', 'Druze' : 'Middle East', 'French' : 'Europe', 
		'Han' : 'East Asia', 'Han-NChina' : 'East Asia', 'Hazara' : 'Central and South Asia', 
		'Hezhen' : 'East Asia', 'Italian' : 'Europe', 'Japanese' : 'East Asia', 
		'Kalash' : 'Central and South Asia', 'Karitiana' : 'America', 'Lahu' : 'East Asia', 
		'Makrani' : 'Central and South Asia', 'Mandenka' : 'Africa', 'Maya' : 'America', 
		'MbutiPygmy' : 'Africa', 'Melanesian' : 'Oceania', 'Miao' : 'East Asia', 
		'Mongola' : 'East Asia', 'Mozabite' : 'Middle East', 'Naxi' : 'East Asia', 
		'Orcadian' : 'Europe', 'Oroqen' : 'East Asia', 'Palestinian' : 'Middle East', 
		'Papuan' : 'Oceania', 'Pathan' : 'Central and South Asia', 'Pima' : 'America', 
		'Russian' : 'Europe', 'San' : 'Africa', 'Sardinian' : 'Europe', 'She' : 'East Asia', 
		'Sindhi' : 'Central and South Asia', 'Surui' : 'America', 'Tu' : 'East Asia', 
		'Tujia' : 'East Asia', 'Tuscan' : 'Europe', 'Uygur' : 'Central and South Asia', 
		'Xibo' : 'East Asia', 'Yakut' : 'East Asia', 'Yi' : 'East Asia', 'Yoruba' : 'Africa'}
study_site_dict = {'101':'Site_Bogota', '103':'Site_Cali', '106':'Site_Cali', '107':'Site_Pereira'}
select_continents = ['Europe', 'America', 'Africa']



##methods
def make_bed_from_site_file(site_file, bed_file):
	with open(site_file, 'r') as site_fh, open(bed_file, 'w') as bed_fh:
		line_count = 0
		for line in site_fh:
			line_count +=1
			if line_count > 1:
				line = line.rstrip().split(delim)
				chrom = line[0]
				start = str(int(line[1]) - 1)
				end = line[1]
				# print chrom, start, end
				bed_fh.write(delim.join([chrom, start, end, '\n']))

def samtools_pileup(bam_files):
	'samtools mpileup -q 30 -Q 20 -f ../../LASER-resource/reference/hs37d5.fa -l HGDP_938.bed exampleBAM/NA12878.chrom22.recal.bam > NA12878.chrom22.pileup'
	##iterate over bamfiles
	for bam in bam_files:
		sample = bam.split('.')[0]
		pileup_file = sample + '.pileup'
		with open(pileup_file, 'w') as pu_fh:
			bcf_index = subprocess.Popen(['samtools', 'mpileup', '-q', '30', '-Q', '20', '-f', fasta, '-l', hgdp_bed, bam], stdout=pu_fh)
			bcf_index.wait()

def pileup_names_from_bams(bam_files):
	pileup_names = []
	for bam in bam_files:
		pileup_file = bam.split('.')[0] + '.pileup'
		pileup_names.append(pileup_file)
	return pileup_names

def combine_mpileup_to_seq_file(mpileup_files, out_prefix):
	'python pileup2seq.py -f hs37d5.fa -m ref.site -b target.bed -i example.id -o output A.pileup B.pileup C.pileup'
	run_pileup2seq = subprocess.Popen(['python', pileup2seq, '-f', fasta, '-m', hgdp_site, '-o', out_prefix] + mpileup_files)
	run_pileup2seq.wait()

def run_laser(geno_file, coord_file, file_prefix, pcs_req):
	run_laser_cmd = subprocess.Popen([laser, '-g', geno_file, '-c', coord_file, '-k', pcs_req, '-s', file_prefix + '.seq', '-o', file_prefix])
	run_laser_cmd.wait()

def combine_format_coord_files_for_r(sample_file_prefix, coord_file, output_file, populations_req):
	sample_file = sample_file_prefix + '.SeqPC.coord'
	with open(coord_file, "r") as ctl_fh, open(sample_file, "r") as samp_fh, open(output_file, 'w') as out_fh:
		##write header
		out_fh.write(delim.join(['continent', 'population', 'sample', 'pc1', 'pc2', 'pc3', 'pc4', 'pc5', 'pc6', 'pc7', 'pc8', '\n']))
		line_count = 0
		for line in ctl_fh:
			line = line.strip('\n').split(delim)
			line_count += 1
			if line_count > 1:
				popID = line[0]
				contID = hgdp_continent_dict[popID]
				# print line[:4]
				if populations_req == 'na':
					out_fh.write(delim.join([contID] + line[:10] + ['\n']))
				else:
					if contID in populations_req:
						out_fh.write(delim.join([contID] + line[:10] + ['\n']))

		line_count = 0
		for line in samp_fh:
			line = line.strip('\n').split(delim)
			line_count += 1
			if line_count > 1:
				sample_id = line[0]
				pcs = line[6:14]
				study_site = study_site_dict[sample_id[:3]]
				# print sample_id, pcs, study_site
				out_fh.write(delim.join(['Sample', study_site, sample_id] + pcs + ['\n']))
##run methods

##get data from bams
##make bed file from site file (only need once)
# make_bed_from_site_file(hgdp_site, hgdp_bed)
##make mpileup files
# samtools_pileup(b1_bams)
# samtools_pileup(b2_bams)
# samtools_pileup(b3_bams)
# samtools_pileup(b4_bams)

##covert pileup and combine in seq file

##batch1
##get list of pileup from bams and then combine into seq file
# pileup_files = pileup_names_from_bams(b1_bams)
# print pileup_files
# combine_mpileup_to_seq_file(pileup_files, b1_prefix)
##run laser
# run_laser(hgdp_geno, hgdp_coord, b1_prefix, '10')
# ##combine sample and control coord files for graphing
# combine_format_coord_files_for_r(b1_prefix, hgdp_coord, b1_hgdp_all_for_r, 'na')
# combine_format_coord_files_for_r(b1_prefix, hgdp_coord, b1_hgdp_select_for_r, select_continents)
# combine_format_coord_files_for_r(b1_prefix, hgdp_coord, b1_hgdp_america_for_r, ['America'])

##combined bams
##get list of pileup from bams and then combine into seq file
pileup_files = pileup_names_from_bams(all_bams)
print pileup_files
combine_mpileup_to_seq_file(pileup_files, b14_prefix)
##run laser
run_laser(hgdp_geno, hgdp_coord, b14_prefix, '10')
##combine sample and control coord files for graphing
combine_format_coord_files_for_r(b14_prefix, hgdp_coord, b14_hgdp_all_for_r, 'na')
combine_format_coord_files_for_r(b14_prefix, hgdp_coord, b14_hgdp_select_for_r, select_continents)
combine_format_coord_files_for_r(b14_prefix, hgdp_coord, b14_hgdp_america_for_r, ['America'])
