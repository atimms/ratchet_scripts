#!/tools/BioBuilds-2015.04/bin/python
import sys
import subprocess
import os
import glob

##parameters
delim = '\t'
##working directory
# working_dir = '/data/atimms/microtia_exomes/cnv_analysis'
working_dir = '/data/atimms/microtia_exomes/cnv_analysis_combined_bed'
os.chdir(working_dir)


##file names etc
ref_dir = '/data/atimms/references/'
fasta = ref_dir + 'human_g1k_v37.fasta'
##which bed file?
# exome_capture_bed = ref_dir + 'dobyns_exome.in_any_padded_target.1015.bed'
# exome_bed_for_fishing = ref_dir + 'dobyns_exome.in_all_targets.added_col.0416.bed'
exome_bed_for_fishing = ref_dir + 'CCDS.current.no_chr.bed'
##file suffixes
final_bam_suffix = '.bwa_gatk.bam'
cov_suffix = '.coverage'
rpkm_suffix = '.rpkm'

##so cp *.bwa_gatk.bam and *.bwa_gatk.bai files to this dir then run script
# bam_dirs = ['/data/atimms/microtia_exomes/batch1', '/data/atimms/microtia_exomes/batch2', '/data/atimms/microtia_exomes/batch3']
bam_dirs = ['/data/atimms/microtia_exomes/temp_bams']
##batch1-3
'''
control_rpkm_files = ['101000202.rpkm', '101000203.rpkm', '101000302.rpkm', '101000402.rpkm', '101000403.rpkm', '101000502.rpkm', '101000503.rpkm', 
		'101000703.rpkm', '103000102.rpkm', '103000103.rpkm', '101001002.rpkm', '101001003.rpkm', '101001302.rpkm', '101001303.rpkm', '101001902.rpkm', 
		'101001903.rpkm', '101002003.rpkm', '101002102.rpkm', '101002103.rpkm', '101002202.rpkm', '101002302.rpkm', '101002402.rpkm', '101002502.rpkm', 
		'101002503.rpkm', '101002602.rpkm', '101002603.rpkm', '101002802.rpkm', '101002803.rpkm', '103000202.rpkm', '103000203.rpkm', '103000602.rpkm', 
		'103000603.rpkm', '103000702.rpkm', '106000402.rpkm', '106000403.rpkm', '106001502.rpkm', '107000103.rpkm', '107000202.rpkm', 
		'107000203.rpkm', '107000302.rpkm', '107000303.rpkm', '107001102.rpkm', '107001103.rpkm']
all_rpkm_files = ['101000201.rpkm', '101000202.rpkm', '101000203.rpkm', '101000301.rpkm', '101000302.rpkm', '101000401.rpkm', '101000402.rpkm', 
		'101000403.rpkm', '101000501.rpkm', '101000502.rpkm', '101000503.rpkm', '101000601.rpkm', '101000602.rpkm', '101000701.rpkm', '101000702.rpkm', 
		'101000703.rpkm', '101000801.rpkm', '101000802.rpkm', '101000803.rpkm', '101001001.rpkm', '101001002.rpkm', '101001003.rpkm', '101001301.rpkm', 
		'101001302.rpkm', '101001303.rpkm', '101001901.rpkm', '101001902.rpkm', '101001903.rpkm', '101002001.rpkm', '101002002.rpkm', '101002003.rpkm', 
		'101002101.rpkm', '101002102.rpkm', '101002103.rpkm', '101002201.rpkm', '101002202.rpkm', '101002203.rpkm', '101002301.rpkm', '101002302.rpkm', 
		'101002303.rpkm', '101002401.rpkm', '101002402.rpkm', '101002403.rpkm', '101002501.rpkm', '101002502.rpkm', '101002503.rpkm', '101002601.rpkm', 
		'101002602.rpkm', '101002603.rpkm', '101002801.rpkm', '101002802.rpkm', '101002803.rpkm', '103000101.rpkm', '103000102.rpkm', '103000103.rpkm', 
		'103000201.rpkm', '103000202.rpkm', '103000203.rpkm', '103000601.rpkm', '103000602.rpkm', '103000603.rpkm', '103000701.rpkm', '103000702.rpkm', 
		'106000401.rpkm', '106000402.rpkm', '106000403.rpkm', '106001501.rpkm', '106001502.rpkm', '106001503.rpkm', '107000101.rpkm', '107000102.rpkm', 
		'107000103.rpkm', '107000201.rpkm', '107000202.rpkm', '107000203.rpkm', '107000301.rpkm', '107000302.rpkm', '107000303.rpkm', '107001101.rpkm', 
		'107001102.rpkm', '107001103.rpkm']
'''
##batch1-4
# '''
control_rpkm_files = ['101000202.rpkm', '101000203.rpkm', '101000302.rpkm', '101000402.rpkm', '101000403.rpkm', '101000502.rpkm', '101000503.rpkm', 
		'101000703.rpkm', '103000102.rpkm', '103000103.rpkm', '101001002.rpkm', '101001003.rpkm', '101001302.rpkm', '101001303.rpkm', '101001902.rpkm', 
		'101001903.rpkm', '101002003.rpkm', '101002102.rpkm', '101002103.rpkm', '101002202.rpkm', '101002302.rpkm', '101002402.rpkm', '101002502.rpkm', 
		'101002503.rpkm', '101002602.rpkm', '101002603.rpkm', '101002802.rpkm', '101002803.rpkm', '103000202.rpkm', '103000203.rpkm', '103000602.rpkm', 
		'103000603.rpkm', '103000702.rpkm', '106000402.rpkm', '106000403.rpkm', '106001502.rpkm', '107000103.rpkm', '107000202.rpkm', 
		'107000203.rpkm', '107000302.rpkm', '107000303.rpkm', '107001102.rpkm', '107001103.rpkm', '101003103.rpkm', '101003202.rpkm', '101003203.rpkm', 
		'103001302.rpkm', '103001303.rpkm', '103001402.rpkm', '103001403.rpkm', '103001502.rpkm', '103001503.rpkm', '106001806b.rpkm', '106002002.rpkm', 
		'106002003.rpkm', '107001202.rpkm', '107001203.rpkm', '107001302.rpkm', '107001303.rpkm', '107001502.rpkm', '107001503.rpkm']
all_rpkm_files = ['101000201.rpkm', '101000202.rpkm', '101000203.rpkm', '101000301.rpkm', '101000302.rpkm', '101000401.rpkm', '101000402.rpkm', '101000403.rpkm', 
		'101000501.rpkm', '101000502.rpkm', '101000503.rpkm', '101000601.rpkm', '101000602.rpkm', '101000701.rpkm', '101000702.rpkm', '101000703.rpkm', '101000801.rpkm', 
		'101000802.rpkm', '101000803.rpkm', '101001001.rpkm', '101001002.rpkm', '101001003.rpkm', '101001301.rpkm', '101001302.rpkm', '101001303.rpkm', '101001901.rpkm', 
		'101001902.rpkm', '101001903.rpkm', '101002001.rpkm', '101002002.rpkm', '101002003.rpkm', '101002101.rpkm', '101002102.rpkm', '101002103.rpkm', '101002201.rpkm', 
		'101002202.rpkm', '101002203.rpkm', '101002301.rpkm', '101002302.rpkm', '101002303.rpkm', '101002401.rpkm', '101002402.rpkm', '101002403.rpkm', '101002501.rpkm', 
		'101002502.rpkm', '101002503.rpkm', '101002601.rpkm', '101002602.rpkm', '101002603.rpkm', '101002801.rpkm', '101002802.rpkm', '101002803.rpkm', '101003101.rpkm', 
		'101003102.rpkm', '101003103.rpkm', '101003201.rpkm', '101003202.rpkm', '101003203.rpkm', '103000101.rpkm', '103000102.rpkm', '103000103.rpkm', '103000201.rpkm', 
		'103000202.rpkm', '103000203.rpkm', '103000601.rpkm', '103000602.rpkm', '103000603.rpkm', '103000701.rpkm', '103000702.rpkm', '103001301.rpkm', '103001302.rpkm', 
		'103001303.rpkm', '103001401.rpkm', '103001402.rpkm', '103001403.rpkm', '103001501.rpkm', '103001502.rpkm', '103001503.rpkm', '106000401.rpkm', '106000402.rpkm', 
		'106000403.rpkm', '106001501.rpkm', '106001502.rpkm', '106001503.rpkm', '106001801.rpkm', '106001806a.rpkm', '106001806b.rpkm', '106001809.rpkm', '106002001.rpkm', 
		'106002002.rpkm', '106002003.rpkm', '107000101.rpkm', '107000102.rpkm', '107000103.rpkm', '107000201.rpkm', '107000202.rpkm', '107000203.rpkm', '107000301.rpkm', 
		'107000302.rpkm', '107000303.rpkm', '107001101.rpkm', '107001102.rpkm', '107001103.rpkm', '107001201.rpkm', '107001202.rpkm', '107001203.rpkm', '107001301.rpkm', 
		'107001302.rpkm', '107001303.rpkm', '107001501.rpkm', '107001502.rpkm', '107001503.rpkm']
# '''
control_rpkm = 'microtia_exomes.controls.rpkm'



##programs
gatk = '/home/atimms/programs/FishingCNV_2.1_pipeline/GenomeAnalysisTK-2.3-9/GenomeAnalysisTK.jar'
fishingCNV = '/home/atimms/programs/FishingCNV_2.1_pipeline/FishingCNV_2.1_pipeline.jar'
fishingRscript = '/home/atimms/programs/FishingCNV_2.1_pipeline/FishingCNV.R'


##gatk coverage
def gatk_depth_of_coverage_and_convert_to_rpkm(bam_dir, bam_file_suffix, doc_file_suffix, rpkm_file_suffix, final_dir):
	# java -jar FishingCNV.jar -cc -c ctr#.coverage.sample_interval_summary -b ccds.bed -o ctr#.rpkm
	os.chdir(bam_dir)
	bam_files = glob.glob('*' + bam_file_suffix)
	cov_sum_suffix = doc_file_suffix + '.sample_interval_summary'
	for bam in bam_files:
		sample = bam.replace(bam_file_suffix, '')
		# print bam, sample
		##gatk depth of coverage
		gatk_doc = subprocess.Popen(['java', '-Xmx20g', '-Xms6g', '-jar', gatk, '-T', 'DepthOfCoverage', '-R', fasta, '-I', bam,'-L', exome_bed_for_fishing, '-o', sample + doc_file_suffix, '--minMappingQuality', '15', '--minBaseQuality', '10', '--omitDepthOutputAtEachBase', '--logging_level', 'INFO', '--summaryCoverageThreshold', '5', '--summaryCoverageThreshold', '7', '--summaryCoverageThreshold', '10', '--summaryCoverageThreshold', '15', '--summaryCoverageThreshold', '20', '--summaryCoverageThreshold', '30', '--summaryCoverageThreshold', '50' ])
		gatk_doc.wait()
		##convert to rpkm
		fish_rpkm = subprocess.Popen(['java', '-jar', fishingCNV, '-cc', '-c', sample + cov_sum_suffix, '-b', exome_bed_for_fishing, '-o', sample + rpkm_file_suffix])
		fish_rpkm.wait()
		##move to final dir
		os.rename(sample + rpkm_file_suffix, final_dir + '/' + sample + rpkm_file_suffix)
	##move to working dir
	os.chdir(final_dir)


def combine_ctl_rpkm_files(control_files, ctl_rpkm_file):
	fish_comb = subprocess.Popen(['java', '-jar', fishingCNV, '-p', '-rpkm'] + control_files +  ['-o', ctl_rpkm_file])
	fish_comb.wait()
# java -jar FishingCNV.jar -p -rpkm path/to/rpkm1.rpkm path/to/rpkm2.rpkm ... -o
# test/sample_output_small_control.ctr

def analyze_fishCNV(rpkm_files, ctl_rpkm_file, output_dir, working_dir):
	# Rscript FishingCNV.R -c controls.ctr.complete -v -s sample1.rpkm sample2.rpkm ... -o output/folder -pca
	os.chdir('/home/atimms/programs/FishingCNV_2.1_pipeline')
	dir_rpkm_files = []
	for rpkm_file in rpkm_files:
		new_rpkm_file = working_dir + '/' + rpkm_file
		dir_rpkm_files.append(new_rpkm_file)
	fish_anal = subprocess.Popen(['Rscript', fishingRscript, '-c', working_dir + '/' + ctl_rpkm_file + '.complete', '-v', '-s'] + dir_rpkm_files + ['-o', output_dir, '-pca'] )
	# fish_anal = subprocess.Popen(['Rscript', fishingRscript, '-c', working_dir + '/' + ctl_rpkm_file + '.complete', '-v', '-s'] + dir_rpkm_files + ['-o', output_dir] )
	fish_anal.wait()

def make_bed_sig_cnvs(results_dir, p_req):
	os.chdir(results_dir)
	csv_files = glob.glob('*csv')
	for csv_file in csv_files:
		holm_bed = csv_file.split('.')[0] + '_holm.bed'
		bh_bed = csv_file.split('.')[0] + '_bh.bed'
		print 'using cnv file %s to get significant files %s and %s'%(csv_file, holm_bed, bh_bed)
		with open(csv_file, 'r') as csv_fh, open(holm_bed, 'w') as holm_fh, open(bh_bed, 'w') as bh_fh:
			line_count = 0
			for line in csv_fh:
				line_count +=1
				if line_count != 1:
					line = line.rstrip().split(',')
					chrom = line[2].strip('"')
					start = line[3].strip('"')
					end = line[4].strip('"')
					genes = line[11].strip('"')
					p_holm = float(line[12].strip('"'))
					p_bh = float(line[13].strip('"'))
					# print chrom, start, end, genes, p_holm, p_bh
					if p_holm <= p_req:
						holm_fh.write(delim.join([chrom, start, end, genes, str(p_holm), '\n']))
					if p_bh <= p_req:
						bh_fh.write(delim.join([chrom, start, end, genes, str(p_bh), '\n']))
					# print line
##fishingCNV methods
##step1 - make rpkm files
# for bam_directory in bam_dirs:
# 	gatk_depth_of_coverage_and_convert_to_rpkm(bam_directory, final_bam_suffix, cov_suffix, rpkm_suffix, working_dir)

##step2 - combine control files and then analyze
# results_dir = working_dir
##move to cnv dir
# os.chdir(results_dir)
##combine controls (can do with sudo)
# combine_ctl_rpkm_files(control_rpkm_files, control_rpkm)
##analyze all samples -need to use sudo 
# analyze_fishCNV(all_rpkm_files, control_rpkm, results_dir + '/fishingcnv_results', results_dir)


##step3 - convert results to bed with sig results
make_bed_sig_cnvs(working_dir + '/fishingcnv_results', 0.05)

##get de novo variants




