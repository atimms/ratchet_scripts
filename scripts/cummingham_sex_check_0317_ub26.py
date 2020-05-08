#!/usr/bin/env python
import sys
import subprocess
import os
import glob

##parameters
delim = '\t'



bcftools = 'bcftools'
bgzip = '/home/atimms/programs/htslib-1.3/bgzip'
plink = '/home/atimms/programs/plink'



def plink_relatadness_check(vcf, file_prefix):
	##bzgip vcf file
	# if os.path.isfile(vcf):
	# 	bgzip_cmd = subprocess.Popen([bgzip, vcf])
	# 	bgzip_cmd.wait()
	##correct filtering??
	bcftools_filter = subprocess.Popen([bcftools, 'view', '-f', 'PASS', '-i', "QUAL>50 & MIN(DP)>30", '-V', 'indels', '-o', 'temp_plink.vcf.gz', '-O', 'z', vcf])
	bcftools_filter.wait()
	##generate plink file from vcf
	make_plink = subprocess.Popen([plink, '--vcf', 'temp_plink.vcf.gz', '--vcf-filter', '--const-fid', '--out', 'temp.pass_q50_dp50'])
	make_plink.wait()
	##plink prune by ld
	# plink_prune = subprocess.Popen([plink, '--bfile', 'temp.pass_q50_dp50', '--indep-pairphase', '20000', '20000', '0.5'])
	plink_prune = subprocess.Popen([plink, '--bfile', 'temp.pass_q50_dp50', '--indep', '50', '5', '2'])
	plink_prune.wait()
	##then trim 
	plink_prune_out = subprocess.Popen([plink, '--bfile', 'temp.pass_q50_dp50', '--extract', 'plink.prune.in', '--make-bed', '--out', 'temp.pruned.pass_q50_dp50'])
	plink_prune_out.wait()
	##split x
	plink_split_x = subprocess.Popen([plink, '--bfile', 'temp.pruned.pass_q50_dp50', '--split-x', 'hg19', '--make-bed', '--out', 'temp.pruned.splitx.pass_q50_dp50'])
	plink_split_x.wait()
	##check sex -results in .sexcheck file
	plink_sex = subprocess.Popen([plink, '--bfile', 'temp.pruned.splitx.pass_q50_dp50', '--check-sex', '--out', file_prefix + '.pass_q50_dp50'])
	plink_sex.wait()
	##ibd check
	# plink_ibd = subprocess.Popen([plink,  '--bfile', 'temp.pass_q50_dp50', '--genome', '--out', file_prefix + '.pass_q50_dp50'])
	# plink_ibd.wait()

project = 'cohort1.cases'
vcf_file = 'cohort1.cases.passed.vcf.gz'


plink_relatadness_check(vcf_file, project)