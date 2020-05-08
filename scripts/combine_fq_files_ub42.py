#!/usr/bin/env python
import sys
import subprocess
import os
import glob

##make sure fastq names are correct
def combine_fq_file(fq_dir, samples, r1_identifier, r2_identifier):
	os.chdir(fq_dir)
	for sample in samples:
		print sample
		r1_fq =  sample + '.r1.fq.gz'
		r2_fq = sample + '.r2.fq.gz'
		r1_to_combine = glob.glob('*/*/*/' + sample + '*' + r1_identifier)
		r2_to_combine = glob.glob('*/*/*/' + sample + '*' + r2_identifier)
		print r1_fq, r1_to_combine
		print r2_fq, r2_to_combine
		with open(r1_fq, 'w') as r1_fh:
			cat_files = subprocess.Popen(['cat'] + r1_to_combine, stdout=r1_fh)
			cat_files.wait()
		with open(r2_fq, 'w') as r2_fh:
			cat_files = subprocess.Popen(['cat'] + r2_to_combine, stdout=r2_fh)
			cat_files.wait()


##batch1
# combine_fq_file('/data/atimms/sherr_data_17/ec2_ec3_es1', ['1090-0', '1090-1', '1090-2', '1173-0-B2', '1173-1-E1', '1173-2', '1175-0', '1175-1', '1175-2', '1175-4', '1339-0', '1339-1', '1339-2', '1588-0', '1588-1', '1588-2', '1629-0', '1629-1', '1629-2'])
##batch2
# combine_fq_file('/data/atimms/sherr_data_17/sherr_b2', ['1289-0', '1289-1', '1289-2', '1318-0', '1318-1', '1318-2', '1324-0', '1324-1', '1324-2', '1416-0', '1416-1', '1416-2', '1426-0', '1426-1', '1426-2', '1512-0', '1512-1', '1512-2'])
##batch3
# combine_fq_file('/data/atimms/sherr_data_17/sherr_b3', ['1011-0', '1011-1', '1011-2', '1020-0', '1020-1', '1020-2', '1200-0', '1200-1', '1200-2', '1536-0', '1536-1', '1536-2', '1590-0', '1590-1', '1590-2'])
##batch4 (on ub26)
# combine_fq_file('/data/atimms/sherr_exomes/sherr_b4', ['1201-0', '1201-1', '1201-2', '1328-0', '1328-1', '1328-2', '1795-0', '1795-1', '1795-2', '2352-0', '2352-1', '2352-2'])

##novartis exomes on 26

combine_fq_file('/data/atimms/novartis_fcd_exomes_0317', ['LIB-01551', 'LIB-01552', 'LIB-01553', 'LIB-01554', 'LIB-01555', 'LIB-01556', 'LIB-01557', 'LIB-01558', 
		'LIB-01559', 'LIB-01560', 'LIB-01561', 'LIB-01562', 'LIB-01563', 'LIB-01564', 'LIB-01565', 'LIB-01566', 'LIB-01567', 'LIB-01568', 'LIB-01569', 'LIB-01570', 
		'LIB-01571', 'LIB-01572', 'LIB-01573', 'LIB-01574', 'LIB-01575', 'LIB-01576', 'LIB-01577', 'LIB-01578', 'LIB-01579', 'LIB-01580', 'LIB-01581', 'LIB-01582', 
		'LIB-01850', 'LIB-03057'], 'R1_001.fastq.gz', 'R2_001.fastq.gz')


##manual cat
with open('LIB-01572_LIB-03057.r1.fq.gz', 'w') as r1_fh:
	cat1 = subprocess.Popen(['cat', 'LIB-01572.r1.fq.gz', 'LIB-03057.r1.fq.gz'], stdout=r1_fh)
	cat1.wait()

with open('LIB-01572_LIB-03057.r2.fq.gz', 'w') as r2_fh:
	cat1 = subprocess.Popen(['cat', 'LIB-01572.r2.fq.gz', 'LIB-03057.r2.fq.gz'], stdout=r2_fh)
	cat1.wait()