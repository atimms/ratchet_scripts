#!/tools/BioBuilds-2015.04/bin/python
import dobyns_exome_pipeline_ub26_v4

##call all methods, supply:
##working dir, file prefix (usually proband nams), if starting from bam or fastq, sample dict, pedigree type

# dobyns_exome_pipeline_v4.call_all_exome_methods_inc_gemini(working_dir, '', 'fastq', {'':['', ''], '':['', ''], '':['', '']}, 'trio')




working_dir = '/data/atimms/zarate_bams_1016'
# dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, 'zarate1', 'bam', {'ACH51154':'Zarate-1-ACH51154_S1.bam', 'ACH51155':'Zarate-1-ACH51155_S1.bam', 'ACH51156':'Zarate-1-ACH51156_S1.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v4.call_just_gemini(working_dir, 'zarate1', 'bam', {'ACH51154':'Zarate-1-ACH51154_S1.bam', 'ACH51155':'Zarate-1-ACH51155_S1.bam', 'ACH51156':'Zarate-1-ACH51156_S1.bam'}, 'trio')

dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, 'zarate2', 'bam', {'ACH57538':'Zarate-1-ACH57538_S1.bam', 'ACH57536':'Zarate-1-ACH57536_S1.bam', 'ACH57537':'Zarate-1-ACH57537_S1.bam'}, 'trio')
dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, 'zarate3', 'bam', {'ACH62818':'Zarate-1-ACH62818_S1.bam', 'ACH62820':'Zarate-1-ACH62820_S1.bam', 'ACH62819':'Zarate-1-ACH62819_S1.bam'}, 'trio')
dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, 'zarate4', 'bam', {'ACH91376':'Zarate-1-ACH91376_S1.bam', 'ACH91365':'Zarate-1-ACH91365_S1.bam', 'ACH91354':'Zarate-1-ACH91354_S1.bam'}, 'trio')