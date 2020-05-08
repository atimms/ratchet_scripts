#!/tools/BioBuilds-2015.04/bin/python
import dobyns_exome_pipeline_ub26_v5

##call all methods, supply:
##working dir, file prefix (usually proband nams), if starting from bam or fastq, sample dict, pedigree type

##templates
##from fq and bam on trio
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '', 'fastq', {'':['', ''], '':['', ''], '':['', '']}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '', 'bam', {'':'', '':'', '':''}, 'trio')


##genwiz batch recieved 0916
# working_dir = '/data/atimms/gw_0916'
# dobyns_exome_pipeline_v4.call_all_exome_methods_inc_gemini(working_dir, '', 'fastq', {'':['', ''], '':['', ''], '':['', '']}, 'trio')
# dobyns_exome_pipeline_v4.call_all_exome_methods_inc_gemini(working_dir, 'LR16-065', 'fastq', {'LR16-065_saliva':['SC25-LR16-065-1603027-PostHy_R1_001.fastq.gz', 'SC25-LR16-065-1603027-PostHy_R2_001.fastq.gz'], 'LR16-065_eye':['SC26-LR16-065-1606076-PostHy_R1_001.fastq.gz', 'SC26-LR16-065-1606076-PostHy_R2_001.fastq.gz'], 'LR16-065f':['SC27-LR16-065-1603027-PostHy_R1_001.fastq.gz', 'SC27-LR16-065-1603027-PostHy_R2_001.fastq.gz'], 'LR16-065m':['SC28-LR16-065-1603027-PostHy_R1_001.fastq.gz', 'SC28-LR16-065-1603027-PostHy_R2_001.fastq.gz']}, 'trio')


# working_dir = '/data/atimms/jimmy_trios_0916'

# dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, '160342', 'bam', {'1591501':'1591501.bam', '1591524':'1591524.bam', '1591518':'1591518.bam'}, 'trio')
##gemini didn't work so just repeat that and issue with freebayes and bcftools
# dobyns_exome_pipeline_ub26_v4.temp_methods(working_dir, '160342', 'bam', {'1591501':'1591501.bam', '1591524':'1591524.bam', '1591518':'1591518.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, '162107', 'bam', {'1610661':'1610661.bam', '1610668':'1610668.bam', '1610671':'1610671.bam'}, 'trio')

# working_dir = '/data/atimms/test_1016'
# dobyns_exome_pipeline_ub26_v4.call_just_gemini(working_dir, '160342', 'bam', {'1591501':'1591501.bam', '1591524':'1591524.bam', '1591518':'1591518.bam'}, 'trio')

##3c exomes 1216
# working_dir = '/data/atimms/3c_exomes_1216'
# dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, '3C-10', 'bam', {'3C-10P':'3C-10P.rmdup.bam'}, 'singleton')
# dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, '3C-11', 'bam', {'3C-11P':'3C-11P.rmdup.bam'}, 'singleton')
# dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, '3C-12', 'bam', {'3C-12P':'3C-12P.rmdup.bam'}, 'singleton')
# dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, '3C-14', 'bam', {'3C-14P':'3C-14P.rmdup.bam'}, 'singleton')
# dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, '3C-7', 'bam', {'3C-7P':'3C-7P.rmdup.bam'}, 'singleton')
# dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, '3C-2', 'bam', {'3C-2F':'3C-2F.rmdup.bam', '3C-2M':'3C-2M.rmdup.bam', '3C-2P':'3C-2P.rmdup.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, '3C-3', 'bam', {'3C-3M':'3C-3M.rmdup.bam', '3C-3P':'3C-3P.rmdup.bam'}, 'duo')
# dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, '3C-4', 'bam', {'3C-4F':'3C-4F.rmdup.bam', '3C-4M':'3C-4M.rmdup.bam', '3C-4P':'3C-4P.rmdup.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, '3C-6', 'bam', {'3C-6F':'3C-6F.rmdup.bam', '3C-6M':'3C-6M.rmdup.bam', '3C-6P':'3C-6P.rmdup.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, '3C-8', 'bam', {'3C-8M':'3C-8M.rmdup.bam', '3C-8P':'3C-8P.rmdup.bam', '3C-8S':'3C-8S.rmdup.bam', '3C-8AB':'3C-8AB.rmdup.bam'}, 'trio')

##genedx exomes 0117
# working_dir = '/data/atimms/genedx_0117'
# dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, 'LR12-476', 'bam', {'LR12-476':'1304741.bam', 'LR12-476m':'1304742.bam', 'LR12-476f':'1304743.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, 'LR13-321', 'bam', {'LR13-321':'1323937.bam', 'LR13-321m':'1323938.bam', 'LR13-321f':'1323939.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, 'LR14-092', 'bam', {'LR14-092':'1313072.bam', 'LR14-092m':'1313073.bam', 'LR14-092f':'1313074.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, 'LR15-003', 'bam', {'LR15-003':'1319462.bam', 'LR15-003m':'1319463.bam', 'LR15-003f':'1319464.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, 'LR15-081', 'bam', {'LR15-081':'1332728.bam', 'LR15-081m':'1332729.bam', 'LR15-081f':'1332730.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, 'LR15-097', 'bam', {'LR15-097':'1409397.bam', 'LR15-097m':'1409398.bam', 'LR15-097f':'1409399.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, 'LR12-230', 'bam', {'LR12-230':'1464222.bam', 'LR12-230m':'1464223.bam', 'LR12-230f':'1464224.bam'}, 'trio')
##issue with one of the bams, ready now
# dobyns_exome_pipeline_ub26_v4.call_all_exome_methods_inc_gemini(working_dir, 'LR15-004', 'bam', {'LR15-004':'1418734.bam', 'LR15-004m':'1418735.bam', 'LR15-004f':'1418736.bam'}, 'trio')

##sherr data
##batch3
# working_dir = '/data/atimms/sherr_exomes/sherr_b3'
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '1011', 'fastq', {'1011-0':['1011-0.r1.fq.gz', '1011-0.r2.fq.gz'], '1011-1':['1011-1.r1.fq.gz', '1011-1.r2.fq.gz'], '1011-2':['1011-2.r1.fq.gz', '1011-2.r2.fq.gz']}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '1020', 'fastq', {'1020-0':['1020-0.r1.fq.gz', '1020-0.r2.fq.gz'], '1020-1':['1020-1.r1.fq.gz', '1020-1.r2.fq.gz'], '1020-2':['1020-2.r1.fq.gz', '1020-2.r2.fq.gz']}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '1200', 'fastq', {'1200-0':['1200-0.r1.fq.gz', '1200-0.r2.fq.gz'], '1200-1':['1200-1.r1.fq.gz', '1200-1.r2.fq.gz'], '1200-2':['1200-2.r1.fq.gz', '1200-2.r2.fq.gz']}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '1536', 'fastq', {'1536-0':['1536-0.r1.fq.gz', '1536-0.r2.fq.gz'], '1536-1':['1536-1.r1.fq.gz', '1536-1.r2.fq.gz'], '1536-2':['1536-2.r1.fq.gz', '1536-2.r2.fq.gz']}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '1590', 'fastq', {'1590-0':['1590-0.r1.fq.gz', '1590-0.r2.fq.gz'], '1590-1':['1590-1.r1.fq.gz', '1590-1.r2.fq.gz'], '1590-2':['1590-2.r1.fq.gz', '1590-2.r2.fq.gz']}, 'trio')
##batch4
# working_dir = '/data/atimms/sherr_exomes/sherr_b4'
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '1201', 'fastq', {'1201-0':['1201-0.r1.fq.gz', '1201-0.r2.fq.gz'], '1201-1':['1201-1.r1.fq.gz', '1201-1.r2.fq.gz'], '1201-2':['1201-2.r1.fq.gz', '1201-2.r2.fq.gz']}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '1328', 'fastq', {'1328-0':['1328-0.r1.fq.gz', '1328-0.r2.fq.gz'], '1328-1':['1328-1.r1.fq.gz', '1328-1.r2.fq.gz'], '1328-2':['1328-2.r1.fq.gz', '1328-2.r2.fq.gz']}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '1795', 'fastq', {'1795-0':['1795-0.r1.fq.gz', '1795-0.r2.fq.gz'], '1795-1':['1795-1.r1.fq.gz', '1795-1.r2.fq.gz'], '1795-2':['1795-2.r1.fq.gz', '1795-2.r2.fq.gz']}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '2352', 'fastq', {'2352-0':['2352-0.r1.fq.gz', '2352-0.r2.fq.gz'], '2352-1':['2352-1.r1.fq.gz', '2352-1.r2.fq.gz'], '2352-2':['2352-2.r1.fq.gz', '2352-2.r2.fq.gz']}, 'trio')
##batch5 i.e. uw
# working_dir = '/data/atimms/sherr_exomes/sherr_b5'
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '1036', 'bam', {'1036-0':'1036-0.74706.bam', '1036-1':'1036-1.74707.bam', '1036-2':'1036-2.74708.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '1270', 'bam', {'1270-0':'1270-0.74709.bam', '1270-1':'1270-1.74710.bam', '1270-2':'1270-2.74711.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '1351', 'bam', {'1351-0':'1351-0.74712.bam', '1351-1':'1351-1.74713.bam', '1351-2':'1351-2.74714.bam'}, 'trio')
##batch5 i.e. clinical part1
# working_dir = '/data/atimms/sherr_exomes/sherr_b6'
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '1250', 'bam', {'1250-0':'1250-0_1436154.bam', '1250-1':'1250-1_1436155.bam', '1250-2':'1250-2_1436156.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '1257', 'bam', {'1257-0':'1257-0_1431400.bam', '1257-1':'1257-1_1431401.bam', '1257-2':'1257-2_1431402.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '1327', 'bam', {'1327-0':'1327-0_1416303.bam', '1327-1':'1327-1_1416304.bam', '1327-2':'1327-2_1416305.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '1510', 'bam', {'1510-0':'1510-0_1414811.bam', '1510-1':'1510-1_1414812.bam', '1510-2':'1510-2_1414813.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '1533', 'bam', {'1533-0':'1533-0_1331266.bam', '1533-1':'1533-1_1331267.bam', '1533-2':'1533-2_1331268.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '1630', 'bam', {'1630-0':'1630-0_1312634.bam', '1630-1':'1630-1_1312635.bam', '1630-2':'1630-2_1312636.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '1810', 'bam', {'1810-0':'1810-0_1463264.bam', '1810-1':'1810-1_1463265.bam', '1810-2':'1810-2_1463266.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '1814', 'bam', {'1814-0':'1814-0_1314161.bam', '1814-1':'1814-1_1314162.bam', '1814-2':'1814-2_1314163.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '1841', 'bam', {'1841-0':'1841-0_1337979.bam', '1841-1':'1841-1_1337980.bam', '1841-2':'1841-2_1337981.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '1940', 'bam', {'1940-0':'1940-0_1320702.bam', '1940-1':'1940-1_1320703.bam', '1940-2':'1940-2_1320704.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '1978', 'bam', {'1978-0':'1978-0_1314652.bam', '1978-1':'1978-1_1314653.bam', '1978-2':'1978-2_1314654.bam'}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '2177', 'bam', {'2177-0':'2177-0_1464193.bam', '2177-1':'2177-1_1464194.bam', '2177-2':'2177-2_1464195.bam'}, 'trio')
##mssing bam files
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, '2210', 'bam', {'2210-0':'2210-0_1470668.bam', '2210-1':'', '2210-2':''}, 'trio')

##novartis fcd exomes
# working_dir = '/data/atimms/novartis_fcd_exomes_0317'
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, 'LR12-243', 'fastq', {'LR12-243_brain':['LIB-01552.r1.fq.gz', 'LIB-01552.r2.fq.gz'], 'LR12-243_gl':['LIB-01553.r1.fq.gz', 'LIB-01553.r2.fq.gz'], 'LR12-243f':['LIB-01554.r1.fq.gz', 'LIB-01554.r2.fq.gz'], 'LR12-243m':['LIB-01555.r1.fq.gz', 'LIB-01555.r2.fq.gz']}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, 'LR12-245', 'fastq', {'LR12-245_brain':['LIB-01556.r1.fq.gz', 'LIB-01556.r2.fq.gz'], 'LR12-245_gl':['LIB-01557.r1.fq.gz', 'LIB-01557.r2.fq.gz'], 'LR12-245f':['LIB-01558.r1.fq.gz', 'LIB-01558.r2.fq.gz'], 'LR12-245m':['LIB-01559.r1.fq.gz', 'LIB-01559.r2.fq.gz']}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, 'LR12-250', 'fastq', {'LR12-250_brain':['LIB-01560.r1.fq.gz', 'LIB-01560.r2.fq.gz'], 'LR12-250_gl':['LIB-01561.r1.fq.gz', 'LIB-01561.r2.fq.gz'], 'LR12-250f':['LIB-01562.r1.fq.gz', 'LIB-01562.r2.fq.gz'], 'LR12-250m':['LIB-01563.r1.fq.gz', 'LIB-01563.r2.fq.gz']}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, 'LR12-255', 'fastq', {'LR12-255_brain':['LIB-01564.r1.fq.gz', 'LIB-01564.r2.fq.gz'], 'LR12-255_gl':['LIB-01565.r1.fq.gz', 'LIB-01565.r2.fq.gz'], 'LR12-255f':['LIB-01566.r1.fq.gz', 'LIB-01566.r2.fq.gz'], 'LR12-255m':['LIB-01567.r1.fq.gz', 'LIB-01567.r2.fq.gz']}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, 'LR12-259', 'fastq', {'LR12-259_brain':['LIB-01568.r1.fq.gz', 'LIB-01568.r2.fq.gz'], 'LR12-259_gl':['LIB-01569.r1.fq.gz', 'LIB-01569.r2.fq.gz'], 'LR12-259f':['LIB-01570.r1.fq.gz', 'LIB-01570.r2.fq.gz'], 'LR12-259m':['LIB-01571.r1.fq.gz', 'LIB-01571.r2.fq.gz']}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, 'LR12-260', 'fastq', {'LR12-260_brain':['LIB-01572_LIB-03057.r1.fq.gz', 'LIB-01572_LIB-03057.r2.fq.gz'], 'LR12-260_gl':['LIB-01573.r1.fq.gz', 'LIB-01573.r2.fq.gz'], 'LR12-260f':['LIB-01574.r1.fq.gz', 'LIB-01574.r2.fq.gz'], 'LR12-260m':['LIB-01575.r1.fq.gz', 'LIB-01575.r2.fq.gz']}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, 'LR12-269', 'fastq', {'LR12-269_brain':['LIB-01576.r1.fq.gz', 'LIB-01576.r2.fq.gz'], 'LR12-269_gl':['LIB-01577.r1.fq.gz', 'LIB-01577.r2.fq.gz'], 'LR12-269f':['LIB-01578.r1.fq.gz', 'LIB-01578.r2.fq.gz'], 'LR12-269m':['LIB-01579.r1.fq.gz', 'LIB-01579.r2.fq.gz']}, 'trio')
# dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, 'LR13-354', 'fastq', {'LR13-354_brain':['LIB-01580.r1.fq.gz', 'LIB-01580.r2.fq.gz'], 'LR13-354_gl':['LIB-01581.r1.fq.gz', 'LIB-01581.r2.fq.gz'], 'LR13-354f':['LIB-01582.r1.fq.gz', 'LIB-01582.r2.fq.gz'], 'LR13-354m':['LIB-01551.r1.fq.gz', 'LIB-01551.r2.fq.gz']}, 'trio')

##LR17-170
working_dir = '/data/atimms/dobyns_exomes_after_0317/LR17-170'
dobyns_exome_pipeline_ub26_v5.call_all_exome_methods_inc_gemini(working_dir, 'LR17-170', 'bam', {'LR17-170':'1657569.bam', 'LR17-170f':'1657585.bam', 'LR17-170m':'1657578.bam'}, 'trio')




