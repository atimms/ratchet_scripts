#!/usr/bin/python
import dobyns_exome_pipeline_cybertron_v6
import dobyns_exome_mosaic_pipeline_cybertron_v3

##presume biobuilds.2016.11 is loaded by modules

##call all methods, supply:
##working dir, file prefix (usually proband nams), if starting from bam or fastq, sample dict, pedigree type



##templates -- 
##must have ped file in place (pedname + .ped)
##fq, trio and single
##dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, '', 'fastq', {'':['', ''], '':['', ''], '':['', '']}, 'trio')
##dobyns_exome_pipeline_v1.call_all_exome_methods_inc_gemini(working_dir, '', 'fastq', {'':['', '']}, 'na')
##bam, trio and single
##dobyns_exome_pipeline_v1.call_all_exome_methods_inc_gemini(working_dir, '', 'bam', {'':'', '':'', '':''}, 'trio')
##dobyns_exome_pipeline_v1.call_all_exome_methods_inc_gemini(working_dir, '', 'bam', {'':''}, 'na')
# dobyns_exome_pipeline_v4.call_all_exome_methods_inc_gemini(working_dir, '', 'bam', {'':'', '':'', '':''}, 'trio')
##mosaic
##paired i.e. mutect2
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'ped_paired', ['aff_bam', 'unaff_bam'], 'sex', '', 'paired')
##trio i.e. mosaic hunter
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'ped', ['pro_bam', 'dad_bam', 'mom_bam'], 'sex', 'pro_coverage', 'trio')



##exomes after 0317



##test ped
# working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR06-339', 'fastq', {'LR06-339': ['LR06-339-1208013-GRA3_R1_001.fastq.gz', 'LR06-339-1208013-GRA3_R2_001.fastq.gz'] , 'LR06-339f':['LR06-339f-1208014-GRB3_R1_001.fastq.gz','LR06-339f-1208014-GRB3_R2_001.fastq.gz'], 'LR06-339m':['LR06-339m-1208015-GRC3_R1_001.fastq.gz','LR06-339m-1208015-GRC3_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_after_mapping(working_dir, 'LR06-339', 'fastq', {'LR06-339': ['LR06-339-1208013-GRA3_R1_001.fastq.gz', 'LR06-339-1208013-GRA3_R2_001.fastq.gz'] , 'LR06-339f':['LR06-339f-1208014-GRB3_R1_001.fastq.gz','LR06-339f-1208014-GRB3_R2_001.fastq.gz'], 'LR06-339m':['LR06-339m-1208015-GRC3_R1_001.fastq.gz','LR06-339m-1208015-GRC3_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR06-339', 'fastq', {'LR06-339': ['LR06-339-1208013-GRA3_R1_001.fastq.gz', 'LR06-339-1208013-GRA3_R2_001.fastq.gz'] , 'LR06-339f':['LR06-339f-1208014-GRB3_R1_001.fastq.gz','LR06-339f-1208014-GRB3_R2_001.fastq.gz'], 'LR06-339m':['LR06-339m-1208015-GRC3_R1_001.fastq.gz','LR06-339m-1208015-GRC3_R2_001.fastq.gz']}, 'trio')

##trying to run all at once
'''
print 'running this bit'
##kims trios - 11654
# working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/kims_trios'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR16-451', 'fastq', {'LR16-451':['LR16-451-1703029-GRB10_R1_001.fastq.gz', 'LR16-451-1703029-GRB10_R2_001.fastq.gz'], 'LR16-451f':['LR16-451f-1703030-GRC10_R1_001.fastq.gz', 'LR16-451f-1703030-GRC10_R2_001.fastq.gz'], 'LR16-451m':['LR16-451m-1703031-GRD10_R1_001.fastq.gz', 'LR16-451m-1703031-GRD10_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR16-079', 'fastq', {'LR16-079':['LR16-070-1603022-GRE10_R1_001.fastq.gz', 'LR16-070-1603022-GRE10_R2_001.fastq.gz'], 'LR16-079f':['LR16-070f-1603023-GRF10_R1_001.fastq.gz', 'LR16-070f-1603023-GRF10_R2_001.fastq.gz'], 'LR16-079m':['LR16-070m-1603024-GRG10_R1_001.fastq.gz', 'LR16-070m-1603024-GRG10_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR14-221', 'fastq', {'LR14-221':['LR14-221-1702004-GRH10_R1_001.fastq.gz', 'LR14-221-1702004-GRH10_R2_001.fastq.gz'], 'LR14-221f':['LR14-221f-1702005-GRA11_R1_001.fastq.gz', 'LR14-221f-1702005-GRA11_R2_001.fastq.gz'], 'LR14-221m':['LR14-221m-1702006-GRB11_R1_001.fastq.gz', 'LR14-221m-1702006-GRB11_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR14-071', 'fastq', {'LR14-071':['LR14-071-1403129-GRC11_R1_001.fastq.gz', 'LR14-071-1403129-GRC11_R2_001.fastq.gz'], 'LR14-071f':['LR14-071f-1403130-GRD11_R1_001.fastq.gz', 'LR14-071f-1403130-GRD11_R2_001.fastq.gz'], 'LR14-071m':['LR14-071m-1403131-GRE11_R1_001.fastq.gz', 'LR14-071m-1403131-GRE11_R2_001.fastq.gz']}, 'trio')
##ghayda's peds - b1 - 11603
# working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b1'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR11-233', 'fastq', {'LR11-233a1':['LR11-233a1-1306082-GRB6_R1_001.fastq.gz', 'LR11-233a1-1306082-GRB6_R2_001.fastq.gz'], 'LR11-233a2':['LR11-233a2-1306083-GRC6_R1_001.fastq.gz', 'LR11-233a2-1306083-GRC6_R2_001.fastq.gz'], 'LR11-233f':['LR11-233f-1306084-GRD6_R1_001.fastq.gz', 'LR11-233f-1306084-GRD6_R2_001.fastq.gz'], 'LR11-233m':['LR11-233m-1306085-GRE6_R1_001.fastq.gz', 'LR11-233m-1306085-GRE6_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR13-013', 'fastq', {'LR13-013':['LR13-013.r1.fq.gz', 'LR13-013.r2.fq.gz'], 'LR13-013f':['LR13-013f-1302031-GRA7_R1_001.fastq.gz', 'LR13-013f-1302031-GRA7_R2_001.fastq.gz'], 'LR13-013m':['LR13-013m.r1.fq.gz', 'LR13-013m.r2.fq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR11-003', 'fastq', {'LR11-003':['LR11-003.r1.fq.gz', 'LR11-003.r2.fq.gz'], 'LR11-003f':['LR11-003f-1311097-GRH5_R1_001.fastq.gz', 'LR11-003f-1311097-GRH5_R2_001.fastq.gz'], 'LR11-003m':['LR11-003m-1201089-GRA6_R1_001.fastq.gz', 'LR11-003m-1201089-GRA6_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR11-347', 'fastq', {'LR11-347_1112100':['LR11-347-1112100-GRB2_R1_001.fastq.gz', 'LR11-347-1112100-GRB2_R2_001.fastq.gz'], 'LR11-347_1209066':['LR11-347-1209066-GRC2_R1_001.fastq.gz', 'LR11-347-1209066-GRC2_R2_001.fastq.gz'], 'LR11-347f':['LR11-347f-1207087-GRD2_R1_001.fastq.gz', 'LR11-347f-1207087-GRD2_R2_001.fastq.gz'], 'LR11-347m':['LR11-347m-1112101-GRE2_R1_001.fastq.gz', 'LR11-347m-1112101-GRE2_R2_001.fastq.gz']}, 'trio')
##ghayda's peds - b2 - 11655
working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b2'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR15-061', 'fastq', {'LR15-061_1504023':['LR15-061-1504023-GRB4_R1_001.fastq.gz', 'LR15-061-1504023-GRB4_R2_001.fastq.gz'], 'LR15-061_1701022':['LR15-061-1701022-GRC4_R1_001.fastq.gz', 'LR15-061-1701022-GRC4_R2_001.fastq.gz'], 'LR15-061f':['LR15-061f-1701023-GRD4_R1_001.fastq.gz', 'LR15-061f-1701023-GRD4_R2_001.fastq.gz'], 'LR15-061m':['LR15-061m-1701024-GRE4_R1_001.fastq.gz', 'LR15-061m-1701024-GRE4_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR14-166', 'fastq', {'LR14-166_1406023':['LR14-166-1406023-GRA5_R1_001.fastq.gz', 'LR14-166-1406023-GRA5_R2_001.fastq.gz'], 'LR14-166_1406095':['LR14-166-1406095-GRB5_R1_001.fastq.gz', 'LR14-166-1406095-GRB5_R2_001.fastq.gz'], 'LR14-166f':['LR14-166f-1701035-GRC5_R1_001.fastq.gz', 'LR14-166f-1701035-GRC5_R2_001.fastq.gz'], 'LR14-166m':['LR14-166m-1406024-GRD5_R1_001.fastq.gz', 'LR14-166m-1406024-GRD5_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR14-211', 'fastq', {'LR14-211':['LR14-211-1408047-GRA1_R1_001.fastq.gz', 'LR14-211-1408047-GRA1_R2_001.fastq.gz'], 'LR14-211f':['LR14-211f-1612060-GRB1_R1_001.fastq.gz', 'LR14-211f-1612060-GRB1_R2_001.fastq.gz'], 'LR14-211m':['LR14-211m-1612061-GRC1_R1_001.fastq.gz', 'LR14-211m-1612061-GRC1_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR12-346', 'fastq', {'LR12-346':['LR12-346-1702019-GRD1_R1_001.fastq.gz', 'LR12-346-1702019-GRD1_R2_001.fastq.gz'], 'LR12-346f':['LR12-346f-1212044-GRE1_R1_001.fastq.gz', 'LR12-346f-1212044-GRE1_R2_001.fastq.gz'], 'LR12-346m':['LR12-346m-1212045-GRF1_R1_001.fastq.gz', 'LR12-346m-1212045-GRF1_R2_001.fastq.gz']}, 'trio')
##ghayda's peds - b3 - 11656
# working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b3'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR14-291', 'fastq', {'LR14-291':['LR14-291-1412042-GRG1_R1_001.fastq.gz', 'LR14-291-1412042-GRG1_R2_001.fastq.gz'], 'LR14-291f':['LR14-291f-1412043-GRH1_R1_001.fastq.gz', 'LR14-291f-1412043-GRH1_R2_001.fastq.gz'], 'LR14-291m':['LR14-291m-1412044-GRA2_R1_001.fastq.gz', 'LR14-291m-1412044-GRA2_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR13-160', 'fastq', {'LR13-160':['LR13-160-1312009-GRF2_R1_001.fastq.gz', 'LR13-160-1312009-GRF2_R2_001.fastq.gz'], 'LR13-160f':['LR13-160f-1311114-GRG2_R1_001.fastq.gz', 'LR13-160f-1311114-GRG2_R2_001.fastq.gz'], 'LR13-160m':['LR13-160m-1311115-GRH2_R1_001.fastq.gz', 'LR13-160m-1311115-GRH2_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR13-150', 'fastq', {'LR13-150':['LR13-150-1404116-GRD3_R1_001.fastq.gz', 'LR13-150-1404116-GRD3_R2_001.fastq.gz'], 'LR13-150f':['LR13-150f-1404117-GRE3_R1_001.fastq.gz', 'LR13-150f-1404117-GRE3_R2_001.fastq.gz'], 'LR13-150m':['LR13-150m-1404118-GRF3_R1_001.fastq.gz', 'LR13-150m-1404118-GRF3_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR15-011', 'fastq', {'LR15-011':['LR15-011-1502109-GRG3_R1_001.fastq.gz', 'LR15-011-1502109-GRG3_R2_001.fastq.gz'], 'LR15-011f':['LR15-011f-1502110-GRH3_R1_001.fastq.gz', 'LR15-011f-1502110-GRH3_R2_001.fastq.gz'], 'LR15-011m':['LR15-011m-1506039-A4_R1_001.fastq.gz', 'LR15-011m-1506039-A4_R2_001.fastq.gz']}, 'trio')
##ghayda's peds - b4 - 11657
# working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b4'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR12-354', 'fastq', {'LR12-354':['LR12-354-1209041-GRF4_R1_001.fastq.gz', 'LR12-354-1209041-GRF4_R2_001.fastq.gz'], 'LR12-354f':['LR12-354f-1209042-GRG4_R1_001.fastq.gz', 'LR12-354f-1209042-GRG4_R2_001.fastq.gz'], 'LR12-354m':['LR12-354m-1209043-GRH4_R1_001.fastq.gz', 'LR12-354m-1209043-GRH4_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR15-236', 'fastq', {'LR15-236':['LR15-236-1508056-GRE5_R1_001.fastq.gz', 'LR15-236-1508056-GRE5_R2_001.fastq.gz'], 'LR15-236f':['LR15-236f-1508057-GRF5_R1_001.fastq.gz', 'LR15-236f-1508057-GRF5_R2_001.fastq.gz'], 'LR15-236m':['LR15-236m-1508058-GRG5_R1_001.fastq.gz', 'LR15-236m-1508058-GRG5_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR12-225', 'fastq', {'LR12-225':['LR12-225-1206067-GRF6_R1_001.fastq.gz', 'LR12-225-1206067-GRF6_R2_001.fastq.gz'], 'LR12-225f':['LR12-225f-1206068-GRG6_R1_001.fastq.gz', 'LR12-225f-1206068-GRG6_R2_001.fastq.gz'], 'LR12-225m':['LR12-225m-1206069-GRH6_R1_001.fastq.gz', 'LR12-225m-1206069-GRH6_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR15-155', 'fastq', {'LR15-155':['LR15-155-1506068-GRB7_R1_001.fastq.gz', 'LR15-155-1506068-GRB7_R2_001.fastq.gz'], 'LR15-155f':['LR15-155f-1507092-GRC7_R1_001.fastq.gz', 'LR15-155f-1507092-GRC7_R2_001.fastq.gz'], 'LR15-155m':['LR15-155m-1507093-GRD7_R1_001.fastq.gz', 'LR15-155m-1507093-GRD7_R2_001.fastq.gz']}, 'trio')
##ghayda's peds - b5 - 11658
# working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b5'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR13-300', 'fastq', {'LR13-300':['LR13-300-1605011-GRE7_R1_001.fastq.gz', 'LR13-300-1605011-GRE7_R2_001.fastq.gz'], 'LR13-300f':['LR13-300f-1605012-GRF7_R1_001.fastq.gz', 'LR13-300f-1605012-GRF7_R2_001.fastq.gz'], 'LR13-300m':['LR13-300m-1605013-GRG7_R1_001.fastq.gz', 'LR13-300m-1605013-GRG7_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR11-124', 'fastq', {'LR11-124':['LR11-124-1105054-GRH7_R1_001.fastq.gz', 'LR11-124-1105054-GRH7_R2_001.fastq.gz'], 'LR11-124f':['LR11-124f--1105055-GRA8_R1_001.fastq.gz', 'LR11-124f--1105055-GRA8_R2_001.fastq.gz'], 'LR11-124m':['LR11-124m--1105056-GRB8_R1_001.fastq.gz', 'LR11-124m--1105056-GRB8_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR12-304', 'fastq', {'LR12-304':['LR12-304-1301017-GRC8_R1_001.fastq.gz', 'LR12-304-1301017-GRC8_R2_001.fastq.gz'], 'LR12-304f':['LR12-304f-1301018-GRD8_R1_001.fastq.gz', 'LR12-304f-1301018-GRD8_R2_001.fastq.gz'], 'LR12-304m':['LR12-304m-1301019-GRE8_R1_001.fastq.gz', 'LR12-304m-1301019-GRE8_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR13-061', 'fastq', {'LR13-061':['LR13-061-1303027-GRF8_R1_001.fastq.gz', 'LR13-061-1303027-GRF8_R2_001.fastq.gz'], 'LR13-061f':['LR13-061f-1303028-GRG8_R1_001.fastq.gz', 'LR13-061f-1303028-GRG8_R2_001.fastq.gz'], 'LR13-061m':['LR13-061m-1303029-GRH8_R1_001.fastq.gz', 'LR13-061m-1303029-GRH8_R2_001.fastq.gz']}, 'trio')
##ghayda's peds - b6 - 11659
# working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b6'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR10-246', 'fastq', {'LR10-246':['LR10-246-1308074-GRA9_R1_001.fastq.gz', 'LR10-246-1308074-GRA9_R2_001.fastq.gz'], 'LR10-246f':['LR10-246f-1504075-GRB9_R1_001.fastq.gz', 'LR10-246f-1504075-GRB9_R2_001.fastq.gz'], 'LR10-246m':['LR10-246m-1308075-GRC9_R1_001.fastq.gz', 'LR10-246m-1308075-GRC9_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR12-492', 'fastq', {'LR12-492':['LR12-492-1405059-GRD9_R1_001.fastq.gz', 'LR12-492-1405059-GRD9_R2_001.fastq.gz'], 'LR12-492f':['LR12-492f-1405060-GRE9_R1_001.fastq.gz', 'LR12-492f-1405060-GRE9_R2_001.fastq.gz'], 'LR12-492m':['LR12-492m-1405061-GRF9_R1_001.fastq.gz', 'LR12-492m-1405061-GRF9_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR13-270', 'fastq', {'LR13-270':['LR13-270-1507110-GRG9_R1_001.fastq.gz', 'LR13-270-1507110-GRG9_R2_001.fastq.gz'], 'LR13-270f':['LR13-270f-1507111-GRH9_R1_001.fastq.gz', 'LR13-270f-1507111-GRH9_R2_001.fastq.gz'], 'LR13-270m':['LR13-270m-1507112-GRA10_R1_001.fastq.gz', 'LR13-270m-1507112-GRA10_R2_001.fastq.gz']}, 'trio')
'''
'''
##without gemini
##kims trios - 11891
# working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/kims_trios'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR16-451', 'fastq', {'LR16-451':['LR16-451-1703029-GRB10_R1_001.fastq.gz', 'LR16-451-1703029-GRB10_R2_001.fastq.gz'], 'LR16-451f':['LR16-451f-1703030-GRC10_R1_001.fastq.gz', 'LR16-451f-1703030-GRC10_R2_001.fastq.gz'], 'LR16-451m':['LR16-451m-1703031-GRD10_R1_001.fastq.gz', 'LR16-451m-1703031-GRD10_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR16-079', 'fastq', {'LR16-079':['LR16-070-1603022-GRE10_R1_001.fastq.gz', 'LR16-070-1603022-GRE10_R2_001.fastq.gz'], 'LR16-079f':['LR16-070f-1603023-GRF10_R1_001.fastq.gz', 'LR16-070f-1603023-GRF10_R2_001.fastq.gz'], 'LR16-079m':['LR16-070m-1603024-GRG10_R1_001.fastq.gz', 'LR16-070m-1603024-GRG10_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR14-221', 'fastq', {'LR14-221':['LR14-221-1702004-GRH10_R1_001.fastq.gz', 'LR14-221-1702004-GRH10_R2_001.fastq.gz'], 'LR14-221f':['LR14-221f-1702005-GRA11_R1_001.fastq.gz', 'LR14-221f-1702005-GRA11_R2_001.fastq.gz'], 'LR14-221m':['LR14-221m-1702006-GRB11_R1_001.fastq.gz', 'LR14-221m-1702006-GRB11_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR14-071', 'fastq', {'LR14-071':['LR14-071-1403129-GRC11_R1_001.fastq.gz', 'LR14-071-1403129-GRC11_R2_001.fastq.gz'], 'LR14-071f':['LR14-071f-1403130-GRD11_R1_001.fastq.gz', 'LR14-071f-1403130-GRD11_R2_001.fastq.gz'], 'LR14-071m':['LR14-071m-1403131-GRE11_R1_001.fastq.gz', 'LR14-071m-1403131-GRE11_R2_001.fastq.gz']}, 'trio')
##ghayda's peds - b1 - 11892
# working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b1'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR11-233', 'fastq', {'LR11-233a1':['LR11-233a1-1306082-GRB6_R1_001.fastq.gz', 'LR11-233a1-1306082-GRB6_R2_001.fastq.gz'], 'LR11-233a2':['LR11-233a2-1306083-GRC6_R1_001.fastq.gz', 'LR11-233a2-1306083-GRC6_R2_001.fastq.gz'], 'LR11-233f':['LR11-233f-1306084-GRD6_R1_001.fastq.gz', 'LR11-233f-1306084-GRD6_R2_001.fastq.gz'], 'LR11-233m':['LR11-233m-1306085-GRE6_R1_001.fastq.gz', 'LR11-233m-1306085-GRE6_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR13-013', 'fastq', {'LR13-013':['LR13-013.r1.fq.gz', 'LR13-013.r2.fq.gz'], 'LR13-013f':['LR13-013f-1302031-GRA7_R1_001.fastq.gz', 'LR13-013f-1302031-GRA7_R2_001.fastq.gz'], 'LR13-013m':['LR13-013m.r1.fq.gz', 'LR13-013m.r2.fq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR11-003', 'fastq', {'LR11-003':['LR11-003.r1.fq.gz', 'LR11-003.r2.fq.gz'], 'LR11-003f':['LR11-003f-1311097-GRH5_R1_001.fastq.gz', 'LR11-003f-1311097-GRH5_R2_001.fastq.gz'], 'LR11-003m':['LR11-003m-1201089-GRA6_R1_001.fastq.gz', 'LR11-003m-1201089-GRA6_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR11-347', 'fastq', {'LR11-347_1112100':['LR11-347-1112100-GRB2_R1_001.fastq.gz', 'LR11-347-1112100-GRB2_R2_001.fastq.gz'], 'LR11-347_1209066':['LR11-347-1209066-GRC2_R1_001.fastq.gz', 'LR11-347-1209066-GRC2_R2_001.fastq.gz'], 'LR11-347f':['LR11-347f-1207087-GRD2_R1_001.fastq.gz', 'LR11-347f-1207087-GRD2_R2_001.fastq.gz'], 'LR11-347m':['LR11-347m-1112101-GRE2_R1_001.fastq.gz', 'LR11-347m-1112101-GRE2_R2_001.fastq.gz']}, 'trio')
##ghayda's peds - b2 - 11893
# working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b2'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR15-061', 'fastq', {'LR15-061_1504023':['LR15-061-1504023-GRB4_R1_001.fastq.gz', 'LR15-061-1504023-GRB4_R2_001.fastq.gz'], 'LR15-061_1701022':['LR15-061-1701022-GRC4_R1_001.fastq.gz', 'LR15-061-1701022-GRC4_R2_001.fastq.gz'], 'LR15-061f':['LR15-061f-1701023-GRD4_R1_001.fastq.gz', 'LR15-061f-1701023-GRD4_R2_001.fastq.gz'], 'LR15-061m':['LR15-061m-1701024-GRE4_R1_001.fastq.gz', 'LR15-061m-1701024-GRE4_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR14-166', 'fastq', {'LR14-166_1406023':['LR14-166-1406023-GRA5_R1_001.fastq.gz', 'LR14-166-1406023-GRA5_R2_001.fastq.gz'], 'LR14-166_1406095':['LR14-166-1406095-GRB5_R1_001.fastq.gz', 'LR14-166-1406095-GRB5_R2_001.fastq.gz'], 'LR14-166f':['LR14-166f-1701035-GRC5_R1_001.fastq.gz', 'LR14-166f-1701035-GRC5_R2_001.fastq.gz'], 'LR14-166m':['LR14-166m-1406024-GRD5_R1_001.fastq.gz', 'LR14-166m-1406024-GRD5_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR14-211', 'fastq', {'LR14-211':['LR14-211-1408047-GRA1_R1_001.fastq.gz', 'LR14-211-1408047-GRA1_R2_001.fastq.gz'], 'LR14-211f':['LR14-211f-1612060-GRB1_R1_001.fastq.gz', 'LR14-211f-1612060-GRB1_R2_001.fastq.gz'], 'LR14-211m':['LR14-211m-1612061-GRC1_R1_001.fastq.gz', 'LR14-211m-1612061-GRC1_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR12-346', 'fastq', {'LR12-346':['LR12-346-1702019-GRD1_R1_001.fastq.gz', 'LR12-346-1702019-GRD1_R2_001.fastq.gz'], 'LR12-346f':['LR12-346f-1212044-GRE1_R1_001.fastq.gz', 'LR12-346f-1212044-GRE1_R2_001.fastq.gz'], 'LR12-346m':['LR12-346m-1212045-GRF1_R1_001.fastq.gz', 'LR12-346m-1212045-GRF1_R2_001.fastq.gz']}, 'trio')
##ghayda's peds - b3 - 11894
# working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b3'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR14-291', 'fastq', {'LR14-291':['LR14-291-1412042-GRG1_R1_001.fastq.gz', 'LR14-291-1412042-GRG1_R2_001.fastq.gz'], 'LR14-291f':['LR14-291f-1412043-GRH1_R1_001.fastq.gz', 'LR14-291f-1412043-GRH1_R2_001.fastq.gz'], 'LR14-291m':['LR14-291m-1412044-GRA2_R1_001.fastq.gz', 'LR14-291m-1412044-GRA2_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR13-160', 'fastq', {'LR13-160':['LR13-160-1312009-GRF2_R1_001.fastq.gz', 'LR13-160-1312009-GRF2_R2_001.fastq.gz'], 'LR13-160f':['LR13-160f-1311114-GRG2_R1_001.fastq.gz', 'LR13-160f-1311114-GRG2_R2_001.fastq.gz'], 'LR13-160m':['LR13-160m-1311115-GRH2_R1_001.fastq.gz', 'LR13-160m-1311115-GRH2_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR13-150', 'fastq', {'LR13-150':['LR13-150-1404116-GRD3_R1_001.fastq.gz', 'LR13-150-1404116-GRD3_R2_001.fastq.gz'], 'LR13-150f':['LR13-150f-1404117-GRE3_R1_001.fastq.gz', 'LR13-150f-1404117-GRE3_R2_001.fastq.gz'], 'LR13-150m':['LR13-150m-1404118-GRF3_R1_001.fastq.gz', 'LR13-150m-1404118-GRF3_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR15-011', 'fastq', {'LR15-011':['LR15-011-1502109-GRG3_R1_001.fastq.gz', 'LR15-011-1502109-GRG3_R2_001.fastq.gz'], 'LR15-011f':['LR15-011f-1502110-GRH3_R1_001.fastq.gz', 'LR15-011f-1502110-GRH3_R2_001.fastq.gz'], 'LR15-011m':['LR15-011m-1506039-A4_R1_001.fastq.gz', 'LR15-011m-1506039-A4_R2_001.fastq.gz']}, 'trio')
##ghayda's peds - b4 - 11895
# working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b4'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR12-354', 'fastq', {'LR12-354':['LR12-354-1209041-GRF4_R1_001.fastq.gz', 'LR12-354-1209041-GRF4_R2_001.fastq.gz'], 'LR12-354f':['LR12-354f-1209042-GRG4_R1_001.fastq.gz', 'LR12-354f-1209042-GRG4_R2_001.fastq.gz'], 'LR12-354m':['LR12-354m-1209043-GRH4_R1_001.fastq.gz', 'LR12-354m-1209043-GRH4_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR15-236', 'fastq', {'LR15-236':['LR15-236-1508056-GRE5_R1_001.fastq.gz', 'LR15-236-1508056-GRE5_R2_001.fastq.gz'], 'LR15-236f':['LR15-236f-1508057-GRF5_R1_001.fastq.gz', 'LR15-236f-1508057-GRF5_R2_001.fastq.gz'], 'LR15-236m':['LR15-236m-1508058-GRG5_R1_001.fastq.gz', 'LR15-236m-1508058-GRG5_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR12-225', 'fastq', {'LR12-225':['LR12-225-1206067-GRF6_R1_001.fastq.gz', 'LR12-225-1206067-GRF6_R2_001.fastq.gz'], 'LR12-225f':['LR12-225f-1206068-GRG6_R1_001.fastq.gz', 'LR12-225f-1206068-GRG6_R2_001.fastq.gz'], 'LR12-225m':['LR12-225m-1206069-GRH6_R1_001.fastq.gz', 'LR12-225m-1206069-GRH6_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR15-155', 'fastq', {'LR15-155':['LR15-155-1506068-GRB7_R1_001.fastq.gz', 'LR15-155-1506068-GRB7_R2_001.fastq.gz'], 'LR15-155f':['LR15-155f-1507092-GRC7_R1_001.fastq.gz', 'LR15-155f-1507092-GRC7_R2_001.fastq.gz'], 'LR15-155m':['LR15-155m-1507093-GRD7_R1_001.fastq.gz', 'LR15-155m-1507093-GRD7_R2_001.fastq.gz']}, 'trio')
##ghayda's peds - b5 - 11896
# working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b5'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR13-300', 'fastq', {'LR13-300':['LR13-300-1605011-GRE7_R1_001.fastq.gz', 'LR13-300-1605011-GRE7_R2_001.fastq.gz'], 'LR13-300f':['LR13-300f-1605012-GRF7_R1_001.fastq.gz', 'LR13-300f-1605012-GRF7_R2_001.fastq.gz'], 'LR13-300m':['LR13-300m-1605013-GRG7_R1_001.fastq.gz', 'LR13-300m-1605013-GRG7_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR11-124', 'fastq', {'LR11-124':['LR11-124-1105054-GRH7_R1_001.fastq.gz', 'LR11-124-1105054-GRH7_R2_001.fastq.gz'], 'LR11-124f':['LR11-124f--1105055-GRA8_R1_001.fastq.gz', 'LR11-124f--1105055-GRA8_R2_001.fastq.gz'], 'LR11-124m':['LR11-124m--1105056-GRB8_R1_001.fastq.gz', 'LR11-124m--1105056-GRB8_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR12-304', 'fastq', {'LR12-304':['LR12-304-1301017-GRC8_R1_001.fastq.gz', 'LR12-304-1301017-GRC8_R2_001.fastq.gz'], 'LR12-304f':['LR12-304f-1301018-GRD8_R1_001.fastq.gz', 'LR12-304f-1301018-GRD8_R2_001.fastq.gz'], 'LR12-304m':['LR12-304m-1301019-GRE8_R1_001.fastq.gz', 'LR12-304m-1301019-GRE8_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR13-061', 'fastq', {'LR13-061':['LR13-061-1303027-GRF8_R1_001.fastq.gz', 'LR13-061-1303027-GRF8_R2_001.fastq.gz'], 'LR13-061f':['LR13-061f-1303028-GRG8_R1_001.fastq.gz', 'LR13-061f-1303028-GRG8_R2_001.fastq.gz'], 'LR13-061m':['LR13-061m-1303029-GRH8_R1_001.fastq.gz', 'LR13-061m-1303029-GRH8_R2_001.fastq.gz']}, 'trio')
##ghayda's peds - b6 - 12577
working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b6'
dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR10-246', 'fastq', {'LR10-246':['LR10-246-1308074-GRA9_R1_001.fastq.gz', 'LR10-246-1308074-GRA9_R2_001.fastq.gz'], 'LR10-246f':['LR10-246f-1504075-GRB9_R1_001.fastq.gz', 'LR10-246f-1504075-GRB9_R2_001.fastq.gz'], 'LR10-246m':['LR10-246m-1308075-GRC9_R1_001.fastq.gz', 'LR10-246m-1308075-GRC9_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR12-492', 'fastq', {'LR12-492':['LR12-492-1405059-GRD9_R1_001.fastq.gz', 'LR12-492-1405059-GRD9_R2_001.fastq.gz'], 'LR12-492f':['LR12-492f-1405060-GRE9_R1_001.fastq.gz', 'LR12-492f-1405060-GRE9_R2_001.fastq.gz'], 'LR12-492m':['LR12-492m-1405061-GRF9_R1_001.fastq.gz', 'LR12-492m-1405061-GRF9_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_without_gemini(working_dir, 'LR13-270', 'fastq', {'LR13-270':['LR13-270-1507110-GRG9_R1_001.fastq.gz', 'LR13-270-1507110-GRG9_R2_001.fastq.gz'], 'LR13-270f':['LR13-270f-1507111-GRH9_R1_001.fastq.gz', 'LR13-270f-1507111-GRH9_R2_001.fastq.gz'], 'LR13-270m':['LR13-270m-1507112-GRA10_R1_001.fastq.gz', 'LR13-270m-1507112-GRA10_R2_001.fastq.gz']}, 'trio')
'''

##just gemini
'''
##kims trios - 11891
working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/kims_trios'
##issue
# dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR16-451', 'fastq', {'LR16-451':['LR16-451-1703029-GRB10_R1_001.fastq.gz', 'LR16-451-1703029-GRB10_R2_001.fastq.gz'], 'LR16-451f':['LR16-451f-1703030-GRC10_R1_001.fastq.gz', 'LR16-451f-1703030-GRC10_R2_001.fastq.gz'], 'LR16-451m':['LR16-451m-1703031-GRD10_R1_001.fastq.gz', 'LR16-451m-1703031-GRD10_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR16-079', 'fastq', {'LR16-079':['LR16-070-1603022-GRE10_R1_001.fastq.gz', 'LR16-070-1603022-GRE10_R2_001.fastq.gz'], 'LR16-079f':['LR16-070f-1603023-GRF10_R1_001.fastq.gz', 'LR16-070f-1603023-GRF10_R2_001.fastq.gz'], 'LR16-079m':['LR16-070m-1603024-GRG10_R1_001.fastq.gz', 'LR16-070m-1603024-GRG10_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR14-221', 'fastq', {'LR14-221':['LR14-221-1702004-GRH10_R1_001.fastq.gz', 'LR14-221-1702004-GRH10_R2_001.fastq.gz'], 'LR14-221f':['LR14-221f-1702005-GRA11_R1_001.fastq.gz', 'LR14-221f-1702005-GRA11_R2_001.fastq.gz'], 'LR14-221m':['LR14-221m-1702006-GRB11_R1_001.fastq.gz', 'LR14-221m-1702006-GRB11_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR14-071', 'fastq', {'LR14-071':['LR14-071-1403129-GRC11_R1_001.fastq.gz', 'LR14-071-1403129-GRC11_R2_001.fastq.gz'], 'LR14-071f':['LR14-071f-1403130-GRD11_R1_001.fastq.gz', 'LR14-071f-1403130-GRD11_R2_001.fastq.gz'], 'LR14-071m':['LR14-071m-1403131-GRE11_R1_001.fastq.gz', 'LR14-071m-1403131-GRE11_R2_001.fastq.gz']}, 'trio')
##ghayda's peds - b1 - 11892
working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b1'
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR11-233', 'fastq', {'LR11-233a1':['LR11-233a1-1306082-GRB6_R1_001.fastq.gz', 'LR11-233a1-1306082-GRB6_R2_001.fastq.gz'], 'LR11-233a2':['LR11-233a2-1306083-GRC6_R1_001.fastq.gz', 'LR11-233a2-1306083-GRC6_R2_001.fastq.gz'], 'LR11-233f':['LR11-233f-1306084-GRD6_R1_001.fastq.gz', 'LR11-233f-1306084-GRD6_R2_001.fastq.gz'], 'LR11-233m':['LR11-233m-1306085-GRE6_R1_001.fastq.gz', 'LR11-233m-1306085-GRE6_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR13-013', 'fastq', {'LR13-013':['LR13-013.r1.fq.gz', 'LR13-013.r2.fq.gz'], 'LR13-013f':['LR13-013f-1302031-GRA7_R1_001.fastq.gz', 'LR13-013f-1302031-GRA7_R2_001.fastq.gz'], 'LR13-013m':['LR13-013m.r1.fq.gz', 'LR13-013m.r2.fq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR11-003', 'fastq', {'LR11-003':['LR11-003.r1.fq.gz', 'LR11-003.r2.fq.gz'], 'LR11-003f':['LR11-003f-1311097-GRH5_R1_001.fastq.gz', 'LR11-003f-1311097-GRH5_R2_001.fastq.gz'], 'LR11-003m':['LR11-003m-1201089-GRA6_R1_001.fastq.gz', 'LR11-003m-1201089-GRA6_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR11-347', 'fastq', {'LR11-347_1112100':['LR11-347-1112100-GRB2_R1_001.fastq.gz', 'LR11-347-1112100-GRB2_R2_001.fastq.gz'], 'LR11-347_1209066':['LR11-347-1209066-GRC2_R1_001.fastq.gz', 'LR11-347-1209066-GRC2_R2_001.fastq.gz'], 'LR11-347f':['LR11-347f-1207087-GRD2_R1_001.fastq.gz', 'LR11-347f-1207087-GRD2_R2_001.fastq.gz'], 'LR11-347m':['LR11-347m-1112101-GRE2_R1_001.fastq.gz', 'LR11-347m-1112101-GRE2_R2_001.fastq.gz']}, 'trio')

##ghayda's peds - b2 - 11893
working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b2'
##issue
# dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR15-061', 'fastq', {'LR15-061_1504023':['LR15-061-1504023-GRB4_R1_001.fastq.gz', 'LR15-061-1504023-GRB4_R2_001.fastq.gz'], 'LR15-061_1701022':['LR15-061-1701022-GRC4_R1_001.fastq.gz', 'LR15-061-1701022-GRC4_R2_001.fastq.gz'], 'LR15-061f':['LR15-061f-1701023-GRD4_R1_001.fastq.gz', 'LR15-061f-1701023-GRD4_R2_001.fastq.gz'], 'LR15-061m':['LR15-061m-1701024-GRE4_R1_001.fastq.gz', 'LR15-061m-1701024-GRE4_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR14-166', 'fastq', {'LR14-166_1406023':['LR14-166-1406023-GRA5_R1_001.fastq.gz', 'LR14-166-1406023-GRA5_R2_001.fastq.gz'], 'LR14-166_1406095':['LR14-166-1406095-GRB5_R1_001.fastq.gz', 'LR14-166-1406095-GRB5_R2_001.fastq.gz'], 'LR14-166f':['LR14-166f-1701035-GRC5_R1_001.fastq.gz', 'LR14-166f-1701035-GRC5_R2_001.fastq.gz'], 'LR14-166m':['LR14-166m-1406024-GRD5_R1_001.fastq.gz', 'LR14-166m-1406024-GRD5_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR14-211', 'fastq', {'LR14-211':['LR14-211-1408047-GRA1_R1_001.fastq.gz', 'LR14-211-1408047-GRA1_R2_001.fastq.gz'], 'LR14-211f':['LR14-211f-1612060-GRB1_R1_001.fastq.gz', 'LR14-211f-1612060-GRB1_R2_001.fastq.gz'], 'LR14-211m':['LR14-211m-1612061-GRC1_R1_001.fastq.gz', 'LR14-211m-1612061-GRC1_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR12-346', 'fastq', {'LR12-346':['LR12-346-1702019-GRD1_R1_001.fastq.gz', 'LR12-346-1702019-GRD1_R2_001.fastq.gz'], 'LR12-346f':['LR12-346f-1212044-GRE1_R1_001.fastq.gz', 'LR12-346f-1212044-GRE1_R2_001.fastq.gz'], 'LR12-346m':['LR12-346m-1212045-GRF1_R1_001.fastq.gz', 'LR12-346m-1212045-GRF1_R2_001.fastq.gz']}, 'trio')
##ghayda's peds - b3 - 
working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b3'
##issue
# dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR14-291', 'fastq', {'LR14-291':['LR14-291-1412042-GRG1_R1_001.fastq.gz', 'LR14-291-1412042-GRG1_R2_001.fastq.gz'], 'LR14-291f':['LR14-291f-1412043-GRH1_R1_001.fastq.gz', 'LR14-291f-1412043-GRH1_R2_001.fastq.gz'], 'LR14-291m':['LR14-291m-1412044-GRA2_R1_001.fastq.gz', 'LR14-291m-1412044-GRA2_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR13-160', 'fastq', {'LR13-160':['LR13-160-1312009-GRF2_R1_001.fastq.gz', 'LR13-160-1312009-GRF2_R2_001.fastq.gz'], 'LR13-160f':['LR13-160f-1311114-GRG2_R1_001.fastq.gz', 'LR13-160f-1311114-GRG2_R2_001.fastq.gz'], 'LR13-160m':['LR13-160m-1311115-GRH2_R1_001.fastq.gz', 'LR13-160m-1311115-GRH2_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR13-150', 'fastq', {'LR13-150':['LR13-150-1404116-GRD3_R1_001.fastq.gz', 'LR13-150-1404116-GRD3_R2_001.fastq.gz'], 'LR13-150f':['LR13-150f-1404117-GRE3_R1_001.fastq.gz', 'LR13-150f-1404117-GRE3_R2_001.fastq.gz'], 'LR13-150m':['LR13-150m-1404118-GRF3_R1_001.fastq.gz', 'LR13-150m-1404118-GRF3_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR15-011', 'fastq', {'LR15-011':['LR15-011-1502109-GRG3_R1_001.fastq.gz', 'LR15-011-1502109-GRG3_R2_001.fastq.gz'], 'LR15-011f':['LR15-011f-1502110-GRH3_R1_001.fastq.gz', 'LR15-011f-1502110-GRH3_R2_001.fastq.gz'], 'LR15-011m':['LR15-011m-1506039-A4_R1_001.fastq.gz', 'LR15-011m-1506039-A4_R2_001.fastq.gz']}, 'trio')

##ghayda's peds - b4 - 
working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b4'
##issue
# dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR12-354', 'fastq', {'LR12-354':['LR12-354-1209041-GRF4_R1_001.fastq.gz', 'LR12-354-1209041-GRF4_R2_001.fastq.gz'], 'LR12-354f':['LR12-354f-1209042-GRG4_R1_001.fastq.gz', 'LR12-354f-1209042-GRG4_R2_001.fastq.gz'], 'LR12-354m':['LR12-354m-1209043-GRH4_R1_001.fastq.gz', 'LR12-354m-1209043-GRH4_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR15-236', 'fastq', {'LR15-236':['LR15-236-1508056-GRE5_R1_001.fastq.gz', 'LR15-236-1508056-GRE5_R2_001.fastq.gz'], 'LR15-236f':['LR15-236f-1508057-GRF5_R1_001.fastq.gz', 'LR15-236f-1508057-GRF5_R2_001.fastq.gz'], 'LR15-236m':['LR15-236m-1508058-GRG5_R1_001.fastq.gz', 'LR15-236m-1508058-GRG5_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR12-225', 'fastq', {'LR12-225':['LR12-225-1206067-GRF6_R1_001.fastq.gz', 'LR12-225-1206067-GRF6_R2_001.fastq.gz'], 'LR12-225f':['LR12-225f-1206068-GRG6_R1_001.fastq.gz', 'LR12-225f-1206068-GRG6_R2_001.fastq.gz'], 'LR12-225m':['LR12-225m-1206069-GRH6_R1_001.fastq.gz', 'LR12-225m-1206069-GRH6_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR15-155', 'fastq', {'LR15-155':['LR15-155-1506068-GRB7_R1_001.fastq.gz', 'LR15-155-1506068-GRB7_R2_001.fastq.gz'], 'LR15-155f':['LR15-155f-1507092-GRC7_R1_001.fastq.gz', 'LR15-155f-1507092-GRC7_R2_001.fastq.gz'], 'LR15-155m':['LR15-155m-1507093-GRD7_R1_001.fastq.gz', 'LR15-155m-1507093-GRD7_R2_001.fastq.gz']}, 'trio')

##ghayda's peds - b5 - 
working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b5'
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR13-300', 'fastq', {'LR13-300':['LR13-300-1605011-GRE7_R1_001.fastq.gz', 'LR13-300-1605011-GRE7_R2_001.fastq.gz'], 'LR13-300f':['LR13-300f-1605012-GRF7_R1_001.fastq.gz', 'LR13-300f-1605012-GRF7_R2_001.fastq.gz'], 'LR13-300m':['LR13-300m-1605013-GRG7_R1_001.fastq.gz', 'LR13-300m-1605013-GRG7_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR11-124', 'fastq', {'LR11-124':['LR11-124-1105054-GRH7_R1_001.fastq.gz', 'LR11-124-1105054-GRH7_R2_001.fastq.gz'], 'LR11-124f':['LR11-124f--1105055-GRA8_R1_001.fastq.gz', 'LR11-124f--1105055-GRA8_R2_001.fastq.gz'], 'LR11-124m':['LR11-124m--1105056-GRB8_R1_001.fastq.gz', 'LR11-124m--1105056-GRB8_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR12-304', 'fastq', {'LR12-304':['LR12-304-1301017-GRC8_R1_001.fastq.gz', 'LR12-304-1301017-GRC8_R2_001.fastq.gz'], 'LR12-304f':['LR12-304f-1301018-GRD8_R1_001.fastq.gz', 'LR12-304f-1301018-GRD8_R2_001.fastq.gz'], 'LR12-304m':['LR12-304m-1301019-GRE8_R1_001.fastq.gz', 'LR12-304m-1301019-GRE8_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR13-061', 'fastq', {'LR13-061':['LR13-061-1303027-GRF8_R1_001.fastq.gz', 'LR13-061-1303027-GRF8_R2_001.fastq.gz'], 'LR13-061f':['LR13-061f-1303028-GRG8_R1_001.fastq.gz', 'LR13-061f-1303028-GRG8_R2_001.fastq.gz'], 'LR13-061m':['LR13-061m-1303029-GRH8_R1_001.fastq.gz', 'LR13-061m-1303029-GRH8_R2_001.fastq.gz']}, 'trio')
##ghayda's peds - b6 - 
working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b6'
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR10-246', 'fastq', {'LR10-246':['LR10-246-1308074-GRA9_R1_001.fastq.gz', 'LR10-246-1308074-GRA9_R2_001.fastq.gz'], 'LR10-246f':['LR10-246f-1504075-GRB9_R1_001.fastq.gz', 'LR10-246f-1504075-GRB9_R2_001.fastq.gz'], 'LR10-246m':['LR10-246m-1308075-GRC9_R1_001.fastq.gz', 'LR10-246m-1308075-GRC9_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR12-492', 'fastq', {'LR12-492':['LR12-492-1405059-GRD9_R1_001.fastq.gz', 'LR12-492-1405059-GRD9_R2_001.fastq.gz'], 'LR12-492f':['LR12-492f-1405060-GRE9_R1_001.fastq.gz', 'LR12-492f-1405060-GRE9_R2_001.fastq.gz'], 'LR12-492m':['LR12-492m-1405061-GRF9_R1_001.fastq.gz', 'LR12-492m-1405061-GRF9_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR13-270', 'fastq', {'LR13-270':['LR13-270-1507110-GRG9_R1_001.fastq.gz', 'LR13-270-1507110-GRG9_R2_001.fastq.gz'], 'LR13-270f':['LR13-270f-1507111-GRH9_R1_001.fastq.gz', 'LR13-270f-1507111-GRH9_R2_001.fastq.gz'], 'LR13-270m':['LR13-270m-1507112-GRA10_R1_001.fastq.gz', 'LR13-270m-1507112-GRA10_R2_001.fastq.gz']}, 'trio')


##a lot of vcfs didn't get intersected so rpt
# working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/kims_trios'
# dobyns_exome_pipeline_cybertron_v6.intersect_vcfs_and_gemini(working_dir, 'LR16-451', 'fastq', {'LR16-451':['LR16-451-1703029-GRB10_R1_001.fastq.gz', 'LR16-451-1703029-GRB10_R2_001.fastq.gz'], 'LR16-451f':['LR16-451f-1703030-GRC10_R1_001.fastq.gz', 'LR16-451f-1703030-GRC10_R2_001.fastq.gz'], 'LR16-451m':['LR16-451m-1703031-GRD10_R1_001.fastq.gz', 'LR16-451m-1703031-GRD10_R2_001.fastq.gz']}, 'trio')
# working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b2'
# dobyns_exome_pipeline_cybertron_v6.intersect_vcfs_and_gemini(working_dir, 'LR15-061', 'fastq', {'LR15-061_1504023':['LR15-061-1504023-GRB4_R1_001.fastq.gz', 'LR15-061-1504023-GRB4_R2_001.fastq.gz'], 'LR15-061_1701022':['LR15-061-1701022-GRC4_R1_001.fastq.gz', 'LR15-061-1701022-GRC4_R2_001.fastq.gz'], 'LR15-061f':['LR15-061f-1701023-GRD4_R1_001.fastq.gz', 'LR15-061f-1701023-GRD4_R2_001.fastq.gz'], 'LR15-061m':['LR15-061m-1701024-GRE4_R1_001.fastq.gz', 'LR15-061m-1701024-GRE4_R2_001.fastq.gz']}, 'trio')
working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b3'
dobyns_exome_pipeline_cybertron_v6.intersect_vcfs_and_gemini(working_dir, 'LR14-291', 'fastq', {'LR14-291':['LR14-291-1412042-GRG1_R1_001.fastq.gz', 'LR14-291-1412042-GRG1_R2_001.fastq.gz'], 'LR14-291f':['LR14-291f-1412043-GRH1_R1_001.fastq.gz', 'LR14-291f-1412043-GRH1_R2_001.fastq.gz'], 'LR14-291m':['LR14-291m-1412044-GRA2_R1_001.fastq.gz', 'LR14-291m-1412044-GRA2_R2_001.fastq.gz']}, 'trio')
working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b4'
dobyns_exome_pipeline_cybertron_v6.intersect_vcfs_and_gemini(working_dir, 'LR12-354', 'fastq', {'LR12-354':['LR12-354-1209041-GRF4_R1_001.fastq.gz', 'LR12-354-1209041-GRF4_R2_001.fastq.gz'], 'LR12-354f':['LR12-354f-1209042-GRG4_R1_001.fastq.gz', 'LR12-354f-1209042-GRG4_R2_001.fastq.gz'], 'LR12-354m':['LR12-354m-1209043-GRH4_R1_001.fastq.gz', 'LR12-354m-1209043-GRH4_R2_001.fastq.gz']}, 'trio')
'''

'''
##rpt plink on all
##kims trios - 11891
working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/kims_trios'
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR16-451', 'fastq', {'LR16-451':['LR16-451-1703029-GRB10_R1_001.fastq.gz', 'LR16-451-1703029-GRB10_R2_001.fastq.gz'], 'LR16-451f':['LR16-451f-1703030-GRC10_R1_001.fastq.gz', 'LR16-451f-1703030-GRC10_R2_001.fastq.gz'], 'LR16-451m':['LR16-451m-1703031-GRD10_R1_001.fastq.gz', 'LR16-451m-1703031-GRD10_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR16-079', 'fastq', {'LR16-079':['LR16-070-1603022-GRE10_R1_001.fastq.gz', 'LR16-070-1603022-GRE10_R2_001.fastq.gz'], 'LR16-079f':['LR16-070f-1603023-GRF10_R1_001.fastq.gz', 'LR16-070f-1603023-GRF10_R2_001.fastq.gz'], 'LR16-079m':['LR16-070m-1603024-GRG10_R1_001.fastq.gz', 'LR16-070m-1603024-GRG10_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR14-221', 'fastq', {'LR14-221':['LR14-221-1702004-GRH10_R1_001.fastq.gz', 'LR14-221-1702004-GRH10_R2_001.fastq.gz'], 'LR14-221f':['LR14-221f-1702005-GRA11_R1_001.fastq.gz', 'LR14-221f-1702005-GRA11_R2_001.fastq.gz'], 'LR14-221m':['LR14-221m-1702006-GRB11_R1_001.fastq.gz', 'LR14-221m-1702006-GRB11_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR14-071', 'fastq', {'LR14-071':['LR14-071-1403129-GRC11_R1_001.fastq.gz', 'LR14-071-1403129-GRC11_R2_001.fastq.gz'], 'LR14-071f':['LR14-071f-1403130-GRD11_R1_001.fastq.gz', 'LR14-071f-1403130-GRD11_R2_001.fastq.gz'], 'LR14-071m':['LR14-071m-1403131-GRE11_R1_001.fastq.gz', 'LR14-071m-1403131-GRE11_R2_001.fastq.gz']}, 'trio')
##ghayda's peds - b1 - 11892
working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b1'
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR11-233', 'fastq', {'LR11-233a1':['LR11-233a1-1306082-GRB6_R1_001.fastq.gz', 'LR11-233a1-1306082-GRB6_R2_001.fastq.gz'], 'LR11-233a2':['LR11-233a2-1306083-GRC6_R1_001.fastq.gz', 'LR11-233a2-1306083-GRC6_R2_001.fastq.gz'], 'LR11-233f':['LR11-233f-1306084-GRD6_R1_001.fastq.gz', 'LR11-233f-1306084-GRD6_R2_001.fastq.gz'], 'LR11-233m':['LR11-233m-1306085-GRE6_R1_001.fastq.gz', 'LR11-233m-1306085-GRE6_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR13-013', 'fastq', {'LR13-013':['LR13-013.r1.fq.gz', 'LR13-013.r2.fq.gz'], 'LR13-013f':['LR13-013f-1302031-GRA7_R1_001.fastq.gz', 'LR13-013f-1302031-GRA7_R2_001.fastq.gz'], 'LR13-013m':['LR13-013m.r1.fq.gz', 'LR13-013m.r2.fq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR11-003', 'fastq', {'LR11-003':['LR11-003.r1.fq.gz', 'LR11-003.r2.fq.gz'], 'LR11-003f':['LR11-003f-1311097-GRH5_R1_001.fastq.gz', 'LR11-003f-1311097-GRH5_R2_001.fastq.gz'], 'LR11-003m':['LR11-003m-1201089-GRA6_R1_001.fastq.gz', 'LR11-003m-1201089-GRA6_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR11-347', 'fastq', {'LR11-347_1112100':['LR11-347-1112100-GRB2_R1_001.fastq.gz', 'LR11-347-1112100-GRB2_R2_001.fastq.gz'], 'LR11-347_1209066':['LR11-347-1209066-GRC2_R1_001.fastq.gz', 'LR11-347-1209066-GRC2_R2_001.fastq.gz'], 'LR11-347f':['LR11-347f-1207087-GRD2_R1_001.fastq.gz', 'LR11-347f-1207087-GRD2_R2_001.fastq.gz'], 'LR11-347m':['LR11-347m-1112101-GRE2_R1_001.fastq.gz', 'LR11-347m-1112101-GRE2_R2_001.fastq.gz']}, 'trio')
##ghayda's peds - b2 - 11893
working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b2'
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR15-061', 'fastq', {'LR15-061_1504023':['LR15-061-1504023-GRB4_R1_001.fastq.gz', 'LR15-061-1504023-GRB4_R2_001.fastq.gz'], 'LR15-061_1701022':['LR15-061-1701022-GRC4_R1_001.fastq.gz', 'LR15-061-1701022-GRC4_R2_001.fastq.gz'], 'LR15-061f':['LR15-061f-1701023-GRD4_R1_001.fastq.gz', 'LR15-061f-1701023-GRD4_R2_001.fastq.gz'], 'LR15-061m':['LR15-061m-1701024-GRE4_R1_001.fastq.gz', 'LR15-061m-1701024-GRE4_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR14-166', 'fastq', {'LR14-166_1406023':['LR14-166-1406023-GRA5_R1_001.fastq.gz', 'LR14-166-1406023-GRA5_R2_001.fastq.gz'], 'LR14-166_1406095':['LR14-166-1406095-GRB5_R1_001.fastq.gz', 'LR14-166-1406095-GRB5_R2_001.fastq.gz'], 'LR14-166f':['LR14-166f-1701035-GRC5_R1_001.fastq.gz', 'LR14-166f-1701035-GRC5_R2_001.fastq.gz'], 'LR14-166m':['LR14-166m-1406024-GRD5_R1_001.fastq.gz', 'LR14-166m-1406024-GRD5_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR14-211', 'fastq', {'LR14-211':['LR14-211-1408047-GRA1_R1_001.fastq.gz', 'LR14-211-1408047-GRA1_R2_001.fastq.gz'], 'LR14-211f':['LR14-211f-1612060-GRB1_R1_001.fastq.gz', 'LR14-211f-1612060-GRB1_R2_001.fastq.gz'], 'LR14-211m':['LR14-211m-1612061-GRC1_R1_001.fastq.gz', 'LR14-211m-1612061-GRC1_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR12-346', 'fastq', {'LR12-346':['LR12-346-1702019-GRD1_R1_001.fastq.gz', 'LR12-346-1702019-GRD1_R2_001.fastq.gz'], 'LR12-346f':['LR12-346f-1212044-GRE1_R1_001.fastq.gz', 'LR12-346f-1212044-GRE1_R2_001.fastq.gz'], 'LR12-346m':['LR12-346m-1212045-GRF1_R1_001.fastq.gz', 'LR12-346m-1212045-GRF1_R2_001.fastq.gz']}, 'trio')
##ghayda's peds - b3 - 11894
working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b3'
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR14-291', 'fastq', {'LR14-291':['LR14-291-1412042-GRG1_R1_001.fastq.gz', 'LR14-291-1412042-GRG1_R2_001.fastq.gz'], 'LR14-291f':['LR14-291f-1412043-GRH1_R1_001.fastq.gz', 'LR14-291f-1412043-GRH1_R2_001.fastq.gz'], 'LR14-291m':['LR14-291m-1412044-GRA2_R1_001.fastq.gz', 'LR14-291m-1412044-GRA2_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR13-160', 'fastq', {'LR13-160':['LR13-160-1312009-GRF2_R1_001.fastq.gz', 'LR13-160-1312009-GRF2_R2_001.fastq.gz'], 'LR13-160f':['LR13-160f-1311114-GRG2_R1_001.fastq.gz', 'LR13-160f-1311114-GRG2_R2_001.fastq.gz'], 'LR13-160m':['LR13-160m-1311115-GRH2_R1_001.fastq.gz', 'LR13-160m-1311115-GRH2_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR13-150', 'fastq', {'LR13-150':['LR13-150-1404116-GRD3_R1_001.fastq.gz', 'LR13-150-1404116-GRD3_R2_001.fastq.gz'], 'LR13-150f':['LR13-150f-1404117-GRE3_R1_001.fastq.gz', 'LR13-150f-1404117-GRE3_R2_001.fastq.gz'], 'LR13-150m':['LR13-150m-1404118-GRF3_R1_001.fastq.gz', 'LR13-150m-1404118-GRF3_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR15-011', 'fastq', {'LR15-011':['LR15-011-1502109-GRG3_R1_001.fastq.gz', 'LR15-011-1502109-GRG3_R2_001.fastq.gz'], 'LR15-011f':['LR15-011f-1502110-GRH3_R1_001.fastq.gz', 'LR15-011f-1502110-GRH3_R2_001.fastq.gz'], 'LR15-011m':['LR15-011m-1506039-A4_R1_001.fastq.gz', 'LR15-011m-1506039-A4_R2_001.fastq.gz']}, 'trio')
##ghayda's peds - b4 - 11895
working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b4'
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR12-354', 'fastq', {'LR12-354':['LR12-354-1209041-GRF4_R1_001.fastq.gz', 'LR12-354-1209041-GRF4_R2_001.fastq.gz'], 'LR12-354f':['LR12-354f-1209042-GRG4_R1_001.fastq.gz', 'LR12-354f-1209042-GRG4_R2_001.fastq.gz'], 'LR12-354m':['LR12-354m-1209043-GRH4_R1_001.fastq.gz', 'LR12-354m-1209043-GRH4_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR15-236', 'fastq', {'LR15-236':['LR15-236-1508056-GRE5_R1_001.fastq.gz', 'LR15-236-1508056-GRE5_R2_001.fastq.gz'], 'LR15-236f':['LR15-236f-1508057-GRF5_R1_001.fastq.gz', 'LR15-236f-1508057-GRF5_R2_001.fastq.gz'], 'LR15-236m':['LR15-236m-1508058-GRG5_R1_001.fastq.gz', 'LR15-236m-1508058-GRG5_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR12-225', 'fastq', {'LR12-225':['LR12-225-1206067-GRF6_R1_001.fastq.gz', 'LR12-225-1206067-GRF6_R2_001.fastq.gz'], 'LR12-225f':['LR12-225f-1206068-GRG6_R1_001.fastq.gz', 'LR12-225f-1206068-GRG6_R2_001.fastq.gz'], 'LR12-225m':['LR12-225m-1206069-GRH6_R1_001.fastq.gz', 'LR12-225m-1206069-GRH6_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR15-155', 'fastq', {'LR15-155':['LR15-155-1506068-GRB7_R1_001.fastq.gz', 'LR15-155-1506068-GRB7_R2_001.fastq.gz'], 'LR15-155f':['LR15-155f-1507092-GRC7_R1_001.fastq.gz', 'LR15-155f-1507092-GRC7_R2_001.fastq.gz'], 'LR15-155m':['LR15-155m-1507093-GRD7_R1_001.fastq.gz', 'LR15-155m-1507093-GRD7_R2_001.fastq.gz']}, 'trio')
##ghayda's peds - b5 - 11896
working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b5'
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR13-300', 'fastq', {'LR13-300':['LR13-300-1605011-GRE7_R1_001.fastq.gz', 'LR13-300-1605011-GRE7_R2_001.fastq.gz'], 'LR13-300f':['LR13-300f-1605012-GRF7_R1_001.fastq.gz', 'LR13-300f-1605012-GRF7_R2_001.fastq.gz'], 'LR13-300m':['LR13-300m-1605013-GRG7_R1_001.fastq.gz', 'LR13-300m-1605013-GRG7_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR11-124', 'fastq', {'LR11-124':['LR11-124-1105054-GRH7_R1_001.fastq.gz', 'LR11-124-1105054-GRH7_R2_001.fastq.gz'], 'LR11-124f':['LR11-124f--1105055-GRA8_R1_001.fastq.gz', 'LR11-124f--1105055-GRA8_R2_001.fastq.gz'], 'LR11-124m':['LR11-124m--1105056-GRB8_R1_001.fastq.gz', 'LR11-124m--1105056-GRB8_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR12-304', 'fastq', {'LR12-304':['LR12-304-1301017-GRC8_R1_001.fastq.gz', 'LR12-304-1301017-GRC8_R2_001.fastq.gz'], 'LR12-304f':['LR12-304f-1301018-GRD8_R1_001.fastq.gz', 'LR12-304f-1301018-GRD8_R2_001.fastq.gz'], 'LR12-304m':['LR12-304m-1301019-GRE8_R1_001.fastq.gz', 'LR12-304m-1301019-GRE8_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR13-061', 'fastq', {'LR13-061':['LR13-061-1303027-GRF8_R1_001.fastq.gz', 'LR13-061-1303027-GRF8_R2_001.fastq.gz'], 'LR13-061f':['LR13-061f-1303028-GRG8_R1_001.fastq.gz', 'LR13-061f-1303028-GRG8_R2_001.fastq.gz'], 'LR13-061m':['LR13-061m-1303029-GRH8_R1_001.fastq.gz', 'LR13-061m-1303029-GRH8_R2_001.fastq.gz']}, 'trio')
##ghayda's peds - b6 - 12577
working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/b6'
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR10-246', 'fastq', {'LR10-246':['LR10-246-1308074-GRA9_R1_001.fastq.gz', 'LR10-246-1308074-GRA9_R2_001.fastq.gz'], 'LR10-246f':['LR10-246f-1504075-GRB9_R1_001.fastq.gz', 'LR10-246f-1504075-GRB9_R2_001.fastq.gz'], 'LR10-246m':['LR10-246m-1308075-GRC9_R1_001.fastq.gz', 'LR10-246m-1308075-GRC9_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR12-492', 'fastq', {'LR12-492':['LR12-492-1405059-GRD9_R1_001.fastq.gz', 'LR12-492-1405059-GRD9_R2_001.fastq.gz'], 'LR12-492f':['LR12-492f-1405060-GRE9_R1_001.fastq.gz', 'LR12-492f-1405060-GRE9_R2_001.fastq.gz'], 'LR12-492m':['LR12-492m-1405061-GRF9_R1_001.fastq.gz', 'LR12-492m-1405061-GRF9_R2_001.fastq.gz']}, 'trio')
dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR13-270', 'fastq', {'LR13-270':['LR13-270-1507110-GRG9_R1_001.fastq.gz', 'LR13-270-1507110-GRG9_R2_001.fastq.gz'], 'LR13-270f':['LR13-270f-1507111-GRH9_R1_001.fastq.gz', 'LR13-270f-1507111-GRH9_R2_001.fastq.gz'], 'LR13-270m':['LR13-270m-1507112-GRA10_R1_001.fastq.gz', 'LR13-270m-1507112-GRA10_R2_001.fastq.gz']}, 'trio')
'''

##sample switches
# working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/gw_0517_ghaydas_trios'
##plink wasn't run
# dobyns_exome_pipeline_cybertron_v6.just_plink(working_dir, 'LR06-339', 'fastq', {'LR06-339': ['LR06-339-1208013-GRA3_R1_001.fastq.gz', 'LR06-339-1208013-GRA3_R2_001.fastq.gz'] , 'LR06-339f':['LR06-339f-1208014-GRB3_R1_001.fastq.gz','LR06-339f-1208014-GRB3_R2_001.fastq.gz'], 'LR06-339m':['LR06-339m-1208015-GRC3_R1_001.fastq.gz','LR06-339m-1208015-GRC3_R2_001.fastq.gz']}, 'trio')
##father and proband switched
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR14-291', 'fastq', {'LR14-291f':['LR14-291-1412042-GRG1_R1_001.fastq.gz', 'LR14-291-1412042-GRG1_R2_001.fastq.gz'], 'LR14-291':['LR14-291f-1412043-GRH1_R1_001.fastq.gz', 'LR14-291f-1412043-GRH1_R2_001.fastq.gz'], 'LR14-291m':['LR14-291m-1412044-GRA2_R1_001.fastq.gz', 'LR14-291m-1412044-GRA2_R2_001.fastq.gz']}, 'trio')
# working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0517/gw_0517_ghaydas_trios/duos'
##father not related, run as duo
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR14-166', 'fastq', {'LR14-166_1406023':['LR14-166-1406023-GRA5_R1_001.fastq.gz', 'LR14-166-1406023-GRA5_R2_001.fastq.gz'], 'LR14-166_1406095':['LR14-166-1406095-GRB5_R1_001.fastq.gz', 'LR14-166-1406095-GRB5_R2_001.fastq.gz'], 'LR14-166m':['LR14-166m-1406024-GRD5_R1_001.fastq.gz', 'LR14-166m-1406024-GRD5_R2_001.fastq.gz']}, 'duo')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR12-346', 'fastq', {'LR12-346':['LR12-346-1702019-GRD1_R1_001.fastq.gz', 'LR12-346-1702019-GRD1_R2_001.fastq.gz'], 'LR12-346m':['LR12-346m-1212045-GRF1_R1_001.fastq.gz', 'LR12-346m-1212045-GRF1_R2_001.fastq.gz']}, 'duo')

##eric 0617 exomes
'''
working_dir = '/home/atimms/ngs_data/exomes/eric_0617/162107'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, '162107', 'bam', {'1610661':'1610661.bam', '1610668':'1610668.bam', '1610671':'1610671.bam'}, 'trio')
working_dir = '/home/atimms/ngs_data/exomes/eric_0617/CIIT0002'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'CIIT0002', 'bam', {'SFA1A':'SFA1A.138733.bam', 'SFA1R':'SFA1R.138736.bam', 'SFA1F':'SFA1F.138734.bam', 'SFA1M':'SFA1M.138735.bam'}, 'trio')
working_dir = '/home/atimms/ngs_data/exomes/eric_0617/CIIT0005'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'CIIT0005', 'bam', {'SCID1A':'SCID1A.138729.bam', 'SCID1R':'SCID1R.138732.bam', 'SCID1F':'SCID1F.138730.bam', 'SCID1M':'SCID1M.138731.bam'}, 'trio')
working_dir = '/home/atimms/ngs_data/exomes/eric_0617/CIIT0057'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'CIIT0057', 'bam', {'CIIT0057A1':'CIIT0057A1.144108.bam', 'CIIT0057A2':'CIIT0057A2.144109.bam', 'CIIT0057F':'CIIT0057F.144111.bam', 'CIIT0057M':'CIIT0057M.144110.bam'}, 'trio')
working_dir = '/home/atimms/ngs_data/exomes/eric_0617/SCH006'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'SCH006', 'bam', {'SCH006_3':'SCH006_3.116878.bam', 'SCH006_2':'SCH006_2.116877.bam', 'SCH006_1':'SCH006_1.116876.bam'}, 'trio')
working_dir = '/home/atimms/ngs_data/exomes/eric_0617/alps_1'
dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'alps_1', 'bam', {'CIIT0008A0159SA':'CIIT0008A0159SA.162317.bam', 'CIIT0008R10142BL':'CIIT0008R10142BL.162320.bam', 'CIIT0008R30197SA':'CIIT0008R30197SA.162322.bam', 'CIIT0008R20509SA':'CIIT0008R20509SA.165493.bam', 'CIIT0008F0508SA':'CIIT0008F0508SA.162318.bam', 'CIIT0008M0196SA':'CIIT0008M0196SA.162319.bam'}, 'trio')

##forgot the ped files
working_dir = '/home/atimms/ngs_data/exomes/eric_0617/CIIT0002'
dobyns_exome_pipeline_cybertron_v6.intersect_vcfs_and_gemini(working_dir, 'CIIT0002', 'bam', {'SFA1A':'SFA1A.138733.bam', 'SFA1R':'SFA1R.138736.bam', 'SFA1F':'SFA1F.138734.bam', 'SFA1M':'SFA1M.138735.bam'}, 'trio')
working_dir = '/home/atimms/ngs_data/exomes/eric_0617/CIIT0005'
dobyns_exome_pipeline_cybertron_v6.intersect_vcfs_and_gemini(working_dir, 'CIIT0005', 'bam', {'SCID1A':'SCID1A.138729.bam', 'SCID1R':'SCID1R.138732.bam', 'SCID1F':'SCID1F.138730.bam', 'SCID1M':'SCID1M.138731.bam'}, 'trio')

working_dir = '/home/atimms/ngs_data/exomes/eric_0617/CIIT0057'
dobyns_exome_pipeline_cybertron_v6.intersect_vcfs_and_gemini(working_dir, 'CIIT0057', 'bam', {'CIIT0057A1':'CIIT0057A1.144108.bam', 'CIIT0057A2':'CIIT0057A2.144109.bam', 'CIIT0057F':'CIIT0057F.144111.bam', 'CIIT0057M':'CIIT0057M.144110.bam'}, 'trio')
working_dir = '/home/atimms/ngs_data/exomes/eric_0617/SCH006'
dobyns_exome_pipeline_cybertron_v6.intersect_vcfs_and_gemini(working_dir, 'SCH006', 'bam', {'SCH006_3':'SCH006_3.116878.bam', 'SCH006_2':'SCH006_2.116877.bam', 'SCH006_1':'SCH006_1.116876.bam'}, 'trio')
working_dir = '/home/atimms/ngs_data/exomes/eric_0617/alps_1'
dobyns_exome_pipeline_cybertron_v6.intersect_vcfs_and_gemini(working_dir, 'alps_1', 'bam', {'CIIT0008A0159SA':'CIIT0008A0159SA.162317.bam', 'CIIT0008R10142BL':'CIIT0008R10142BL.162320.bam', 'CIIT0008R30197SA':'CIIT0008R30197SA.162322.bam', 'CIIT0008R20509SA':'CIIT0008R20509SA.165493.bam', 'CIIT0008F0508SA':'CIIT0008F0508SA.162318.bam', 'CIIT0008M0196SA':'CIIT0008M0196SA.162319.bam'}, 'trio')
working_dir = '/home/atimms/ngs_data/exomes/eric_0617/CIIT0057'
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'CIIT0057', 'bam', {'CIIT0057A1':'CIIT0057A1.144108.bam', 'CIIT0057A2':'CIIT0057A2.144109.bam', 'CIIT0057F':'CIIT0057F.144111.bam', 'CIIT0057M':'CIIT0057M.144110.bam'}, 'trio')

##with mom as dominant
working_dir = '/home/atimms/ngs_data/exomes/eric_0617/alps_1_dom'
dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'alps_1', 'bam', {'CIIT0008A0159SA':'CIIT0008A0159SA.162317.bam', 'CIIT0008R10142BL':'CIIT0008R10142BL.162320.bam', 'CIIT0008R30197SA':'CIIT0008R30197SA.162322.bam', 'CIIT0008R20509SA':'CIIT0008R20509SA.165493.bam', 'CIIT0008F0508SA':'CIIT0008F0508SA.162318.bam', 'CIIT0008M0196SA':'CIIT0008M0196SA.162319.bam'}, 'trio')
'''


##ghayda peds 0617 - mosaic analysis
working_dir = '/home/atimms/ngs_data/exomes/dobyns_mosaic_peds/ghayda_0617'
##mutect -- 6/26
##14476
##these 2 still need to run
# dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR11-347_paired', ['LR11-347_1209066.bwa_gatk.bam', 'LR11-347_1112100.bwa_gatk.bam'], 'F', '', 'paired')
# dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR15-061_paired', ['LR15-061_1504023.bwa_gatk.bam', 'LR15-061_1701022.bwa_gatk.bam'], 'M', '', 'paired')
##extra one runnig
# dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR14-226_skin', ['LR14-226_skin.bwa_gatk.bam', 'LR14-226_saliva.bwa_gatk.bam'], 'F', '174', 'paired')

######
# dobyns_exome_mosaic_pipeline_cybertron_v2.run_mosaic_variant_calling(working_dir, 'LR14-166_paired', ['LR14-166_1406095.bwa_gatk.bam', 'LR14-166_1406023.bwa_gatk.bam'], 'F', '', 'paired')
# dobyns_exome_mosaic_pipeline_cybertron_v2.run_mosaic_variant_calling(working_dir, 'LR14-285_skin-affected', ['LR14-285_skin-affected.bwa_gatk.bam', 'LR14-285_skin-normal.bwa_gatk.bam'], 'F', '', 'paired')
# dobyns_exome_mosaic_pipeline_cybertron_v2.run_mosaic_variant_calling(working_dir, 'LR14-285_vasculargrowth', ['LR14-285_vasculargrowth.bwa_gatk.bam', 'LR14-285_skin-normal.bwa_gatk.bam'], 'F', '', 'paired')
# dobyns_exome_mosaic_pipeline_cybertron_v2.run_mosaic_variant_calling(working_dir, 'LR16-053_tumor_fetal', ['LR16-053_Tumor.bwa_gatk.bam', 'LR16-053_FetalTissue.bwa_gatk.bam'], 'sex', '', 'paired')
##done
# dobyns_exome_mosaic_pipeline_cybertron_v2.run_mosaic_variant_calling(working_dir, 'LR16-053_tumor_shoulder', ['LR16-053_Tumor.bwa_gatk.bam', 'LR16-053_leftshoulder.bwa_gatk.bam'], 'sex', '', 'paired')
# dobyns_exome_mosaic_pipeline_cybertron_v2.run_mosaic_variant_calling(working_dir, 'LR16-053_tumor_testis', ['LR16-053_Tumor.bwa_gatk.bam', 'LR16-053_Testis.bwa_gatk.bam'], 'sex', '', 'paired')

##mosaic hunter -- 6/26
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR03-223', ['LR03-223.bwa_gatk.bam', 'LR03-223f.bwa_gatk.bam', 'LR03-223m.bwa_gatk.bam'], 'F', '89', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR03-130', ['LR03-130.bwa_gatk.bam', 'LR03-130f.bwa_gatk.bam', 'LR03-130m.bwa_gatk.bam'], 'M', '100', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR14-226_skin', ['LR14-226_skin.bwa_gatk.bam', 'LR14-226f.bwa_gatk.bam', 'LR14-226m.bwa_gatk.bam'], 'F', '174', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR14-226_saliva', ['LR14-226_saliva.bwa_gatk.bam', 'LR14-226f.bwa_gatk.bam', 'LR14-226m.bwa_gatk.bam'], 'F', '474', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR14-285_skin-affected', ['LR14-285_skin-affected.bwa_gatk.bam', 'LR14-285f.bwa_gatk.bam', 'LR14-285m.bwa_gatk.bam'], 'F', '220', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR14-285_skin-normal', ['LR14-285_skin-normal.bwa_gatk.bam', 'LR14-285f.bwa_gatk.bam', 'LR14-285m.bwa_gatk.bam'], 'F', '225', 'trio')
##13922
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR14-285_vasculargrowth', ['LR14-285_vasculargrowth.bwa_gatk.bam', 'LR14-285f.bwa_gatk.bam', 'LR14-285m.bwa_gatk.bam'], 'F', '212', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR11-347_1112100', ['LR11-347_1112100.bwa_gatk.bam', 'LR11-347f.bwa_gatk.bam', 'LR11-347m.bwa_gatk.bam'], 'F', '408', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR11-347_1209066', ['LR11-347_1209066.bwa_gatk.bam', 'LR11-347f.bwa_gatk.bam', 'LR11-347m.bwa_gatk.bam'], 'F', '384', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR15-061_1504023', ['LR15-061_1504023.bwa_gatk.bam', 'LR15-061f.bwa_gatk.bam', 'LR15-061m.bwa_gatk.bam'], 'M', '441', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR15-061_1701022', ['LR15-061_1701022.bwa_gatk.bam', 'LR15-061f.bwa_gatk.bam', 'LR15-061m.bwa_gatk.bam'], 'M', '301', 'trio')

##vardict -- to rpt with new filters
# dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR03-223', ['LR03-223.bwa_gatk.bam', 'LR03-223f.bwa_gatk.bam', 'LR03-223m.bwa_gatk.bam'], 'F', '89', 'vardict_trio')
# dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR03-130', ['LR03-130.bwa_gatk.bam', 'LR03-130f.bwa_gatk.bam', 'LR03-130m.bwa_gatk.bam'], 'M', '100', 'vardict_trio')
# dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR14-226_skin', ['LR14-226_skin.bwa_gatk.bam', 'LR14-226f.bwa_gatk.bam', 'LR14-226m.bwa_gatk.bam'], 'F', '174', 'vardict_trio')
# dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR14-226_saliva', ['LR14-226_saliva.bwa_gatk.bam', 'LR14-226f.bwa_gatk.bam', 'LR14-226m.bwa_gatk.bam'], 'F', '474', 'vardict_trio')
# dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR14-285_skin-affected', ['LR14-285_skin-affected.bwa_gatk.bam', 'LR14-285f.bwa_gatk.bam', 'LR14-285m.bwa_gatk.bam'], 'F', '220', 'vardict_trio')
# dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR14-285_skin-normal', ['LR14-285_skin-normal.bwa_gatk.bam', 'LR14-285f.bwa_gatk.bam', 'LR14-285m.bwa_gatk.bam'], 'F', '225', 'vardict_trio')
# dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR14-285_vasculargrowth', ['LR14-285_vasculargrowth.bwa_gatk.bam', 'LR14-285f.bwa_gatk.bam', 'LR14-285m.bwa_gatk.bam'], 'F', '212', 'vardict_trio')
# dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR11-347_1112100', ['LR11-347_1112100.bwa_gatk.bam', 'LR11-347f.bwa_gatk.bam', 'LR11-347m.bwa_gatk.bam'], 'F', '408', 'vardict_trio')
# dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR11-347_1209066', ['LR11-347_1209066.bwa_gatk.bam', 'LR11-347f.bwa_gatk.bam', 'LR11-347m.bwa_gatk.bam'], 'F', '384', 'vardict_trio')
# dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR15-061_1504023', ['LR15-061_1504023.bwa_gatk.bam', 'LR15-061f.bwa_gatk.bam', 'LR15-061m.bwa_gatk.bam'], 'M', '441', 'vardict_trio')
# dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR15-061_1701022', ['LR15-061_1701022.bwa_gatk.bam', 'LR15-061f.bwa_gatk.bam', 'LR15-061m.bwa_gatk.bam'], 'M', '301', 'vardict_trio')

##pisces -- to rpt with new filters (14465)
'''
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR03-223', ['LR03-223.bwa_gatk.bam', 'LR03-223f.bwa_gatk.bam', 'LR03-223m.bwa_gatk.bam'], 'F', '89', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR03-130', ['LR03-130.bwa_gatk.bam', 'LR03-130f.bwa_gatk.bam', 'LR03-130m.bwa_gatk.bam'], 'M', '100', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR14-226_skin', ['LR14-226_skin.bwa_gatk.bam', 'LR14-226f.bwa_gatk.bam', 'LR14-226m.bwa_gatk.bam'], 'F', '174', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR14-226_saliva', ['LR14-226_saliva.bwa_gatk.bam', 'LR14-226f.bwa_gatk.bam', 'LR14-226m.bwa_gatk.bam'], 'F', '474', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR14-285_skin-affected', ['LR14-285_skin-affected.bwa_gatk.bam', 'LR14-285f.bwa_gatk.bam', 'LR14-285m.bwa_gatk.bam'], 'F', '220', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR14-285_skin-normal', ['LR14-285_skin-normal.bwa_gatk.bam', 'LR14-285f.bwa_gatk.bam', 'LR14-285m.bwa_gatk.bam'], 'F', '225', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR14-285_vasculargrowth', ['LR14-285_vasculargrowth.bwa_gatk.bam', 'LR14-285f.bwa_gatk.bam', 'LR14-285m.bwa_gatk.bam'], 'F', '212', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR11-347_1112100', ['LR11-347_1112100.bwa_gatk.bam', 'LR11-347f.bwa_gatk.bam', 'LR11-347m.bwa_gatk.bam'], 'F', '408', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR11-347_1209066', ['LR11-347_1209066.bwa_gatk.bam', 'LR11-347f.bwa_gatk.bam', 'LR11-347m.bwa_gatk.bam'], 'F', '384', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR15-061_1504023', ['LR15-061_1504023.bwa_gatk.bam', 'LR15-061f.bwa_gatk.bam', 'LR15-061m.bwa_gatk.bam'], 'M', '441', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR15-061_1701022', ['LR15-061_1701022.bwa_gatk.bam', 'LR15-061f.bwa_gatk.bam', 'LR15-061m.bwa_gatk.bam'], 'M', '301', 'pisces_trio')
'''


##kims peds 0617 - mosaic analysis
working_dir = '/home/atimms/ngs_data/exomes/dobyns_mosaic_peds/kim_0617'
##mosaic hunter
##13866
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, '3C-2P', ['3C-2P.bwa_gatk.bam', '3C-2F.bwa_gatk.bam', '3C-2M.bwa_gatk.bam'], 'M', '155', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, '3C-6P', ['3C-6P.bwa_gatk.bam', '3C-6F.bwa_gatk.bam', '3C-6M.bwa_gatk.bam'], 'F', '147', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR03-305', ['LR03-305.bwa_gatk.bam', 'LR03-305f.bwa_gatk.bam', 'LR03-305m.bwa_gatk.bam'], 'M', '79', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR04-186', ['LR04-186.bwa_gatk.bam', 'LR04-186f.bwa_gatk.bam', 'LR04-186m.bwa_gatk.bam'], 'M', '50', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR05-203a2', ['LR05-203a2_sample2.bwa_gatk.bam ', 'LR05-203f2.bwa_gatk.bam', 'LR05-203m.bwa_gatk.bam'], 'F', '54', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR05-354', ['LR05-354.bwa_gatk.bam', 'LR05-354f.bwa_gatk.bam', 'LR05-354m.bwa_gatk.bam'], 'F', '88', 'trio')
##13867
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR08-396', ['LR08-396.bwa_gatk.bam', 'LR08-396f.bwa_gatk.bam', 'LR08-396m.bwa_gatk.bam'], 'F', '50', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR11-169', ['LR11-169.bwa_gatk.bam', 'LR11-169f.bwa_gatk.bam', 'LR11-169m.bwa_gatk.bam'], 'F', '50', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR12-032', ['LR12-032.bwa_gatk.bam', 'LR12-032f.bwa_gatk.bam', 'LR12-032m.bwa_gatk.bam'], 'F', '55', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR12-115', ['LR12-115.bwa_gatk.bam', 'LR12-115f.bwa_gatk.bam', 'LR12-115m.bwa_gatk.bam'], 'M', '62', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR12-443', ['LR12-443.bwa_gatk.bam', 'LR12-443f.bwa_gatk.bam', 'LR12-443m.bwa_gatk.bam'], 'F', '74', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR12-463', ['LR12-463.bwa_gatk.bam', 'LR12-463f.bwa_gatk.bam', 'LR12-463m.bwa_gatk.bam'], 'M', '55', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR12-464', ['LR12-464.bwa_gatk.bam', 'LR12-464f.bwa_gatk.bam', 'LR12-464m.bwa_gatk.bam'], 'M', '49', 'trio')
##13868
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR13-002', ['LR13-002.bwa_gatk.bam', 'LR13-002f.bwa_gatk.bam', 'LR13-002m.bwa_gatk.bam'], 'M', '50', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR13-085', ['LR13-085.bwa_gatk.bam', 'LR13-085f.bwa_gatk.bam', 'LR13-085m.bwa_gatk.bam'], 'F', '58', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR13-153', ['LR13-153.bwa_gatk.bam', 'LR13-153f.bwa_gatk.bam', 'LR13-153m.bwa_gatk.bam'], 'M', '41', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR13-200', ['LR13-200.bwa_gatk.bam', 'LR13-200f.bwa_gatk.bam', 'LR13-200m.bwa_gatk.bam'], 'M', '45', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR13-315', ['LR13-315.bwa_gatk.bam', 'LR13-315f.bwa_gatk.bam', 'LR13-315m.bwa_gatk.bam'], 'M', '75', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR03-055', ['LR03-055.bwa_gatk.bam', 'LR03-055f.bwa_gatk.bam', 'LR03-055m.bwa_gatk.bam'], 'F', '108', 'trio')
##13869
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR03-274', ['LR03-274.bwa_gatk.bam', 'LR03-274f.bwa_gatk.bam', 'LR03-274m.bwa_gatk.bam'], 'F', '68', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR04-020', ['LR04-020.bwa_gatk.bam', 'LR04-020f.bwa_gatk.bam', 'LR04-020m.bwa_gatk.bam'], 'M', '98', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR04-371', ['LR04-371.bwa_gatk.bam', 'LR04-371f.bwa_gatk.bam', 'LR04-371m.bwa_gatk.bam'], 'M', '95', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR09-280', ['LR09-280.bwa_gatk.bam', 'LR09-280f.bwa_gatk.bam', 'LR09-280m.bwa_gatk.bam'], 'M', '73', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR03-077', ['LR03-077.bwa_gatk.bam', 'LR03-077f.bwa_gatk.bam', 'LR03-077m.bwa_gatk.bam'], 'F', '53', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR03-169', ['LR03-169.bwa_gatk.bam', 'LR03-169f.bwa_gatk.bam', 'LR03-169m.bwa_gatk.bam'], 'M', '66', 'trio')
##13871
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR03-223', ['LR03-223.bwa_gatk.bam', 'LR03-223f.bwa_gatk.bam', 'LR03-223m.bwa_gatk.bam'], 'F', '89', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR03-278', ['LR03-278.bwa_gatk.bam', 'LR03-278f.bwa_gatk.bam', 'LR03-278m.bwa_gatk.bam'], 'M', '52', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR03-298', ['LR03-298.bwa_gatk.bam', 'LR03-298f.bwa_gatk.bam', 'LR03-298m.bwa_gatk.bam'], 'M', '117', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR04-017', ['LR04-017.bwa_gatk.bam', 'LR04-017f.bwa_gatk.bam', 'LR04-017m.bwa_gatk.bam'], 'M', '83', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR04-084', ['LR04-084.bwa_gatk.bam', 'LR04-084f.bwa_gatk.bam', 'LR04-084m.bwa_gatk.bam'], 'F', '117', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR05-118', ['LR05-118.bwa_gatk.bam', 'LR05-118f.bwa_gatk.bam', 'LR05-118m.bwa_gatk.bam'], 'F', '46', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR05-396', ['LR05-396.bwa_gatk.bam', 'LR05-396f.bwa_gatk.bam', 'LR05-396m.bwa_gatk.bam'], 'F', '123', 'trio')
##13873
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR01-079', ['LR01-079.bwa_gatk.bam', 'LR01-079f.bwa_gatk.bam', 'LR01-079m.bwa_gatk.bam'], 'M', '47', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR04-399', ['LR04-399.bwa_gatk.bam', 'LR04-399f.bwa_gatk.bam', 'LR04-399m.bwa_gatk.bam'], 'F', '46', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR05-398', ['LR05-398.bwa_gatk.bam', 'LR05-398f.bwa_gatk.bam', 'LR05-398m.bwa_gatk.bam'], 'M', '41', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR12-313a2', ['LR12-313a2.bwa_gatk.bam', 'LR12-313f.bwa_gatk.bam', 'LR12-313m.bwa_gatk.bam'], 'M', '47', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR12-439', ['LR12-439.bwa_gatk.bam', 'LR12-439f.bwa_gatk.bam', 'LR12-439m.bwa_gatk.bam'], 'M', '46', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR06-157', ['LR06-157.bwa_gatk.bam', 'LR06-157f.bwa_gatk.bam', 'LR06-157m.bwa_gatk.bam'], 'M', '272', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR10-222', ['LR10-222.bwa_gatk.bam', 'LR10-222f.bwa_gatk.bam', 'LR10-222m.bwa_gatk.bam'], 'M', '258', 'trio')
##13874
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR11-152', ['LR11-152.bwa_gatk.bam', 'LR11-152f.bwa_gatk.bam', 'LR11-152m.bwa_gatk.bam'], 'M', '284', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR03-120', ['LR03-120.bwa_gatk.bam', 'LR03-120f.bwa_gatk.bam', 'LR03-120m.bwa_gatk.bam'], 'M', '56', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR04-376', ['LR04-376.bwa_gatk.bam', 'LR04-376f.bwa_gatk.bam', 'LR04-376m.bwa_gatk.bam'], 'M', '53', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR06-278', ['LR06-278.bwa_gatk.bam', 'LR06-278f.bwa_gatk.bam', 'LR06-278m.bwa_gatk.bam'], 'M', '53', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR08-002', ['LR08-002.bwa_gatk.bam', 'LR08-002f.bwa_gatk.bam', 'LR08-002m.bwa_gatk.bam'], 'M', '52', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR12-316', ['LR12-316.bwa_gatk.bam', 'LR12-316f.bwa_gatk.bam', 'LR12-316m.bwa_gatk.bam'], 'F', '95', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR10-026', ['LR10-026.bwa_gatk.bam', 'LR10-026f.bwa_gatk.bam', 'LR10-026m.bwa_gatk.bam'], 'M', '188', 'trio')
# dobyns_exome_mosaic_pipeline_cybertron_v1.run_mosaic_variant_calling(working_dir, 'LR10-228a1', ['LR10-228a1.bwa_gatk.bam', 'LR10-228f.bwa_gatk.bam', 'LR10-228m.bwa_gatk.bam'], 'F', '169', 'trio')

##vardict
##14461
'''
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, '3C-2P', ['3C-2P.bwa_gatk.bam', '3C-2F.bwa_gatk.bam', '3C-2M.bwa_gatk.bam'], 'M', '155', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, '3C-6P', ['3C-6P.bwa_gatk.bam', '3C-6F.bwa_gatk.bam', '3C-6M.bwa_gatk.bam'], 'F', '147', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR03-305', ['LR03-305.bwa_gatk.bam', 'LR03-305f.bwa_gatk.bam', 'LR03-305m.bwa_gatk.bam'], 'M', '79', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR04-186', ['LR04-186.bwa_gatk.bam', 'LR04-186f.bwa_gatk.bam', 'LR04-186m.bwa_gatk.bam'], 'M', '50', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR05-203a2', ['LR05-203a2_sample2.bwa_gatk.bam', 'LR05-203f2.bwa_gatk.bam', 'LR05-203m.bwa_gatk.bam'], 'F', '54', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR05-354', ['LR05-354.bwa_gatk.bam', 'LR05-354f.bwa_gatk.bam', 'LR05-354m.bwa_gatk.bam'], 'F', '88', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR08-396', ['LR08-396.bwa_gatk.bam', 'LR08-396f.bwa_gatk.bam', 'LR08-396m.bwa_gatk.bam'], 'F', '50', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR11-169', ['LR11-169.bwa_gatk.bam', 'LR11-169f.bwa_gatk.bam', 'LR11-169m.bwa_gatk.bam'], 'F', '50', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR12-032', ['LR12-032.bwa_gatk.bam', 'LR12-032f.bwa_gatk.bam', 'LR12-032m.bwa_gatk.bam'], 'F', '55', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR12-115', ['LR12-115.bwa_gatk.bam', 'LR12-115f.bwa_gatk.bam', 'LR12-115m.bwa_gatk.bam'], 'M', '62', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR12-443', ['LR12-443.bwa_gatk.bam', 'LR12-443f.bwa_gatk.bam', 'LR12-443m.bwa_gatk.bam'], 'F', '74', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR12-463', ['LR12-463.bwa_gatk.bam', 'LR12-463f.bwa_gatk.bam', 'LR12-463m.bwa_gatk.bam'], 'M', '55', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR12-464', ['LR12-464.bwa_gatk.bam', 'LR12-464f.bwa_gatk.bam', 'LR12-464m.bwa_gatk.bam'], 'M', '49', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR13-002', ['LR13-002.bwa_gatk.bam', 'LR13-002f.bwa_gatk.bam', 'LR13-002m.bwa_gatk.bam'], 'M', '50', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR13-085', ['LR13-085.bwa_gatk.bam', 'LR13-085f.bwa_gatk.bam', 'LR13-085m.bwa_gatk.bam'], 'F', '58', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR13-153', ['LR13-153.bwa_gatk.bam', 'LR13-153f.bwa_gatk.bam', 'LR13-153m.bwa_gatk.bam'], 'M', '41', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR13-200', ['LR13-200.bwa_gatk.bam', 'LR13-200f.bwa_gatk.bam', 'LR13-200m.bwa_gatk.bam'], 'M', '45', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR13-315', ['LR13-315.bwa_gatk.bam', 'LR13-315f.bwa_gatk.bam', 'LR13-315m.bwa_gatk.bam'], 'M', '75', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR03-055', ['LR03-055.bwa_gatk.bam', 'LR03-055f.bwa_gatk.bam', 'LR03-055m.bwa_gatk.bam'], 'F', '108', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR03-274', ['LR03-274.bwa_gatk.bam', 'LR03-274f.bwa_gatk.bam', 'LR03-274m.bwa_gatk.bam'], 'F', '68', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR04-020', ['LR04-020.bwa_gatk.bam', 'LR04-020f.bwa_gatk.bam', 'LR04-020m.bwa_gatk.bam'], 'M', '98', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR04-371', ['LR04-371.bwa_gatk.bam', 'LR04-371f.bwa_gatk.bam', 'LR04-371m.bwa_gatk.bam'], 'M', '95', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR09-280', ['LR09-280.bwa_gatk.bam', 'LR09-280f.bwa_gatk.bam', 'LR09-280m.bwa_gatk.bam'], 'M', '73', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR03-077', ['LR03-077.bwa_gatk.bam', 'LR03-077f.bwa_gatk.bam', 'LR03-077m.bwa_gatk.bam'], 'F', '53', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR03-169', ['LR03-169.bwa_gatk.bam', 'LR03-169f.bwa_gatk.bam', 'LR03-169m.bwa_gatk.bam'], 'M', '66', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR03-223', ['LR03-223.bwa_gatk.bam', 'LR03-223f.bwa_gatk.bam', 'LR03-223m.bwa_gatk.bam'], 'F', '89', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR03-278', ['LR03-278.bwa_gatk.bam', 'LR03-278f.bwa_gatk.bam', 'LR03-278m.bwa_gatk.bam'], 'M', '52', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR03-298', ['LR03-298.bwa_gatk.bam', 'LR03-298f.bwa_gatk.bam', 'LR03-298m.bwa_gatk.bam'], 'M', '117', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR04-017', ['LR04-017.bwa_gatk.bam', 'LR04-017f.bwa_gatk.bam', 'LR04-017m.bwa_gatk.bam'], 'M', '83', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR04-084', ['LR04-084.bwa_gatk.bam', 'LR04-084f.bwa_gatk.bam', 'LR04-084m.bwa_gatk.bam'], 'F', '117', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR05-118', ['LR05-118.bwa_gatk.bam', 'LR05-118f.bwa_gatk.bam', 'LR05-118m.bwa_gatk.bam'], 'F', '46', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR05-396', ['LR05-396.bwa_gatk.bam', 'LR05-396f.bwa_gatk.bam', 'LR05-396m.bwa_gatk.bam'], 'F', '123', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR01-079', ['LR01-079.bwa_gatk.bam', 'LR01-079f.bwa_gatk.bam', 'LR01-079m.bwa_gatk.bam'], 'M', '47', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR04-399', ['LR04-399.bwa_gatk.bam', 'LR04-399f.bwa_gatk.bam', 'LR04-399m.bwa_gatk.bam'], 'F', '46', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR05-398', ['LR05-398.bwa_gatk.bam', 'LR05-398f.bwa_gatk.bam', 'LR05-398m.bwa_gatk.bam'], 'M', '41', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR12-313a2', ['LR12-313a2.bwa_gatk.bam', 'LR12-313f.bwa_gatk.bam', 'LR12-313m.bwa_gatk.bam'], 'M', '47', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR12-439', ['LR12-439.bwa_gatk.bam', 'LR12-439f.bwa_gatk.bam', 'LR12-439m.bwa_gatk.bam'], 'M', '46', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR06-157', ['LR06-157.bwa_gatk.bam', 'LR06-157f.bwa_gatk.bam', 'LR06-157m.bwa_gatk.bam'], 'M', '272', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR10-222', ['LR10-222.bwa_gatk.bam', 'LR10-222f.bwa_gatk.bam', 'LR10-222m.bwa_gatk.bam'], 'M', '258', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR11-152', ['LR11-152.bwa_gatk.bam', 'LR11-152f.bwa_gatk.bam', 'LR11-152m.bwa_gatk.bam'], 'M', '284', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR03-120', ['LR03-120.bwa_gatk.bam', 'LR03-120f.bwa_gatk.bam', 'LR03-120m.bwa_gatk.bam'], 'M', '56', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR04-376', ['LR04-376.bwa_gatk.bam', 'LR04-376f.bwa_gatk.bam', 'LR04-376m.bwa_gatk.bam'], 'M', '53', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR06-278', ['LR06-278.bwa_gatk.bam', 'LR06-278f.bwa_gatk.bam', 'LR06-278m.bwa_gatk.bam'], 'M', '53', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR08-002', ['LR08-002.bwa_gatk.bam', 'LR08-002f.bwa_gatk.bam', 'LR08-002m.bwa_gatk.bam'], 'M', '52', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR12-316', ['LR12-316.bwa_gatk.bam', 'LR12-316f.bwa_gatk.bam', 'LR12-316m.bwa_gatk.bam'], 'F', '95', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR10-026', ['LR10-026.bwa_gatk.bam', 'LR10-026f.bwa_gatk.bam', 'LR10-026m.bwa_gatk.bam'], 'M', '188', 'vardict_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR10-228a1', ['LR10-228a1.bwa_gatk.bam', 'LR10-228f.bwa_gatk.bam', 'LR10-228m.bwa_gatk.bam'], 'F', '169', 'vardict_trio')
'''
##pisces
##14468
'''
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, '3C-2P', ['3C-2P.bwa_gatk.bam', '3C-2F.bwa_gatk.bam', '3C-2M.bwa_gatk.bam'], 'M', '155', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, '3C-6P', ['3C-6P.bwa_gatk.bam', '3C-6F.bwa_gatk.bam', '3C-6M.bwa_gatk.bam'], 'F', '147', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR03-305', ['LR03-305.bwa_gatk.bam', 'LR03-305f.bwa_gatk.bam', 'LR03-305m.bwa_gatk.bam'], 'M', '79', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR04-186', ['LR04-186.bwa_gatk.bam', 'LR04-186f.bwa_gatk.bam', 'LR04-186m.bwa_gatk.bam'], 'M', '50', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR05-203a2', ['LR05-203a2_sample2.bwa_gatk.bam', 'LR05-203f2.bwa_gatk.bam', 'LR05-203m.bwa_gatk.bam'], 'F', '54', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR05-354', ['LR05-354.bwa_gatk.bam', 'LR05-354f.bwa_gatk.bam', 'LR05-354m.bwa_gatk.bam'], 'F', '88', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR08-396', ['LR08-396.bwa_gatk.bam', 'LR08-396f.bwa_gatk.bam', 'LR08-396m.bwa_gatk.bam'], 'F', '50', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR11-169', ['LR11-169.bwa_gatk.bam', 'LR11-169f.bwa_gatk.bam', 'LR11-169m.bwa_gatk.bam'], 'F', '50', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR12-032', ['LR12-032.bwa_gatk.bam', 'LR12-032f.bwa_gatk.bam', 'LR12-032m.bwa_gatk.bam'], 'F', '55', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR12-115', ['LR12-115.bwa_gatk.bam', 'LR12-115f.bwa_gatk.bam', 'LR12-115m.bwa_gatk.bam'], 'M', '62', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR12-443', ['LR12-443.bwa_gatk.bam', 'LR12-443f.bwa_gatk.bam', 'LR12-443m.bwa_gatk.bam'], 'F', '74', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR12-463', ['LR12-463.bwa_gatk.bam', 'LR12-463f.bwa_gatk.bam', 'LR12-463m.bwa_gatk.bam'], 'M', '55', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR12-464', ['LR12-464.bwa_gatk.bam', 'LR12-464f.bwa_gatk.bam', 'LR12-464m.bwa_gatk.bam'], 'M', '49', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR13-002', ['LR13-002.bwa_gatk.bam', 'LR13-002f.bwa_gatk.bam', 'LR13-002m.bwa_gatk.bam'], 'M', '50', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR13-085', ['LR13-085.bwa_gatk.bam', 'LR13-085f.bwa_gatk.bam', 'LR13-085m.bwa_gatk.bam'], 'F', '58', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR13-153', ['LR13-153.bwa_gatk.bam', 'LR13-153f.bwa_gatk.bam', 'LR13-153m.bwa_gatk.bam'], 'M', '41', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR13-200', ['LR13-200.bwa_gatk.bam', 'LR13-200f.bwa_gatk.bam', 'LR13-200m.bwa_gatk.bam'], 'M', '45', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR13-315', ['LR13-315.bwa_gatk.bam', 'LR13-315f.bwa_gatk.bam', 'LR13-315m.bwa_gatk.bam'], 'M', '75', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR03-055', ['LR03-055.bwa_gatk.bam', 'LR03-055f.bwa_gatk.bam', 'LR03-055m.bwa_gatk.bam'], 'F', '108', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR03-274', ['LR03-274.bwa_gatk.bam', 'LR03-274f.bwa_gatk.bam', 'LR03-274m.bwa_gatk.bam'], 'F', '68', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR04-020', ['LR04-020.bwa_gatk.bam', 'LR04-020f.bwa_gatk.bam', 'LR04-020m.bwa_gatk.bam'], 'M', '98', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR04-371', ['LR04-371.bwa_gatk.bam', 'LR04-371f.bwa_gatk.bam', 'LR04-371m.bwa_gatk.bam'], 'M', '95', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR09-280', ['LR09-280.bwa_gatk.bam', 'LR09-280f.bwa_gatk.bam', 'LR09-280m.bwa_gatk.bam'], 'M', '73', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR03-077', ['LR03-077.bwa_gatk.bam', 'LR03-077f.bwa_gatk.bam', 'LR03-077m.bwa_gatk.bam'], 'F', '53', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR03-169', ['LR03-169.bwa_gatk.bam', 'LR03-169f.bwa_gatk.bam', 'LR03-169m.bwa_gatk.bam'], 'M', '66', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR03-223', ['LR03-223.bwa_gatk.bam', 'LR03-223f.bwa_gatk.bam', 'LR03-223m.bwa_gatk.bam'], 'F', '89', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR03-278', ['LR03-278.bwa_gatk.bam', 'LR03-278f.bwa_gatk.bam', 'LR03-278m.bwa_gatk.bam'], 'M', '52', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR03-298', ['LR03-298.bwa_gatk.bam', 'LR03-298f.bwa_gatk.bam', 'LR03-298m.bwa_gatk.bam'], 'M', '117', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR04-017', ['LR04-017.bwa_gatk.bam', 'LR04-017f.bwa_gatk.bam', 'LR04-017m.bwa_gatk.bam'], 'M', '83', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR04-084', ['LR04-084.bwa_gatk.bam', 'LR04-084f.bwa_gatk.bam', 'LR04-084m.bwa_gatk.bam'], 'F', '117', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR05-118', ['LR05-118.bwa_gatk.bam', 'LR05-118f.bwa_gatk.bam', 'LR05-118m.bwa_gatk.bam'], 'F', '46', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR05-396', ['LR05-396.bwa_gatk.bam', 'LR05-396f.bwa_gatk.bam', 'LR05-396m.bwa_gatk.bam'], 'F', '123', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR01-079', ['LR01-079.bwa_gatk.bam', 'LR01-079f.bwa_gatk.bam', 'LR01-079m.bwa_gatk.bam'], 'M', '47', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR04-399', ['LR04-399.bwa_gatk.bam', 'LR04-399f.bwa_gatk.bam', 'LR04-399m.bwa_gatk.bam'], 'F', '46', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR05-398', ['LR05-398.bwa_gatk.bam', 'LR05-398f.bwa_gatk.bam', 'LR05-398m.bwa_gatk.bam'], 'M', '41', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR12-313a2', ['LR12-313a2.bwa_gatk.bam', 'LR12-313f.bwa_gatk.bam', 'LR12-313m.bwa_gatk.bam'], 'M', '47', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR12-439', ['LR12-439.bwa_gatk.bam', 'LR12-439f.bwa_gatk.bam', 'LR12-439m.bwa_gatk.bam'], 'M', '46', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR06-157', ['LR06-157.bwa_gatk.bam', 'LR06-157f.bwa_gatk.bam', 'LR06-157m.bwa_gatk.bam'], 'M', '272', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR10-222', ['LR10-222.bwa_gatk.bam', 'LR10-222f.bwa_gatk.bam', 'LR10-222m.bwa_gatk.bam'], 'M', '258', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR11-152', ['LR11-152.bwa_gatk.bam', 'LR11-152f.bwa_gatk.bam', 'LR11-152m.bwa_gatk.bam'], 'M', '284', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR03-120', ['LR03-120.bwa_gatk.bam', 'LR03-120f.bwa_gatk.bam', 'LR03-120m.bwa_gatk.bam'], 'M', '56', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR04-376', ['LR04-376.bwa_gatk.bam', 'LR04-376f.bwa_gatk.bam', 'LR04-376m.bwa_gatk.bam'], 'M', '53', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR06-278', ['LR06-278.bwa_gatk.bam', 'LR06-278f.bwa_gatk.bam', 'LR06-278m.bwa_gatk.bam'], 'M', '53', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR08-002', ['LR08-002.bwa_gatk.bam', 'LR08-002f.bwa_gatk.bam', 'LR08-002m.bwa_gatk.bam'], 'M', '52', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR12-316', ['LR12-316.bwa_gatk.bam', 'LR12-316f.bwa_gatk.bam', 'LR12-316m.bwa_gatk.bam'], 'F', '95', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR10-026', ['LR10-026.bwa_gatk.bam', 'LR10-026f.bwa_gatk.bam', 'LR10-026m.bwa_gatk.bam'], 'M', '188', 'pisces_trio')
dobyns_exome_mosaic_pipeline_cybertron_v3.run_mosaic_variant_calling(working_dir, 'LR10-228a1', ['LR10-228a1.bwa_gatk.bam', 'LR10-228f.bwa_gatk.bam', 'LR10-228m.bwa_gatk.bam'], 'F', '169', 'pisces_trio')
'''

##kim/zach trios 
working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0717'
## 14430, 14432 -- 14433
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR01-194', 'fastq', {'LR01-194':['LR01-194_S1_R1_001.fastq.gz', 'LR01-194_S1_R2_001.fastq.gz'], 'LR01-194f':['LR01-194f_S2_R1_001.fastq.gz', 'LR01-194f_S2_R2_001.fastq.gz'], 'LR01-194m':['LR01-194m_S3_R1_001.fastq.gz', 'LR01-194m_S3_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR01-381', 'fastq', {'LR01-381':['LR01-381-2_S40_R1_001.fastq.gz', 'LR01-381-2_S40_R2_001.fastq.gz'], 'LR01-381f':['LR01-381f_S4_R1_001.fastq.gz', 'LR01-381f_S4_R2_001.fastq.gz'], 'LR01-381m':['LR01-381m-2_S41_R1_001.fastq.gz', 'LR01-381m-2_S41_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR02-291', 'fastq', {'LR02-291':['LR02-291_S5_R1_001.fastq.gz', 'LR02-291_S5_R2_001.fastq.gz'], 'LR02-291f':['LR02-291f_S6_R1_001.fastq.gz', 'LR02-291f_S6_R2_001.fastq.gz'], 'LR02-291m':['LR02-291m_S7_R1_001.fastq.gz', 'LR02-291m_S7_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR04-290', 'fastq', {'LR04-290':['LR04-290_S8_R1_001.fastq.gz', 'LR04-290_S8_R2_001.fastq.gz'], 'LR04-290f':['LR04-290f_S9_R1_001.fastq.gz', 'LR04-290f_S9_R2_001.fastq.gz'], 'LR04-290m':['LR04-290m-2_S42_R1_001.fastq.gz', 'LR04-290m-2_S42_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR05-055', 'fastq', {'LR05-055a2':['LR05-055a2_S10_R1_001.fastq.gz', 'LR05-055a2_S10_R2_001.fastq.gz'], 'LR05-055f':['LR05-055f_S11_R1_001.fastq.gz', 'LR05-055f_S11_R2_001.fastq.gz'], 'LR05-055m':['LR05-055m_S12_R1_001.fastq.gz', 'LR05-055m_S12_R2_001.fastq.gz']}, 'trio')
# working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0717/b2'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR05-069', 'fastq', {'LR05-069':['LR05-069_S13_R1_001.fastq.gz', 'LR05-069_S13_R2_001.fastq.gz'], 'LR05-069f':['LR05-069f_S14_R1_001.fastq.gz', 'LR05-069f_S14_R2_001.fastq.gz'], 'LR05-069m':['LR05-069m-2_S43_R1_001.fastq.gz', 'LR05-069m-2_S43_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR05-197', 'fastq', {'LR05-197':['LR05-197_S15_R1_001.fastq.gz', 'LR05-197_S15_R2_001.fastq.gz'], 'LR05-197f':['LR05-197f_S16_R1_001.fastq.gz', 'LR05-197f_S16_R2_001.fastq.gz'], 'LR05-197m':['LR05-197m_S17_R1_001.fastq.gz', 'LR05-197m_S17_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR05-201', 'fastq', {'LR05-201':['LR05-201_S18_R1_001.fastq.gz', 'LR05-201_S18_R2_001.fastq.gz'], 'LR05-201f':['LR05-201f_S19_R1_001.fastq.gz', 'LR05-201f_S19_R2_001.fastq.gz'], 'LR05-201m':['LR05-201m_S20_R1_001.fastq.gz', 'LR05-201m_S20_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR08-299', 'fastq', {'LR08-299':['LR08-299_S21_R1_001.fastq.gz', 'LR08-299_S21_R2_001.fastq.gz'], 'LR08-299f':['LR08-299f_S22_R1_001.fastq.gz', 'LR08-299f_S22_R2_001.fastq.gz'], 'LR08-299m':['LR08-299m_S23_R1_001.fastq.gz', 'LR08-299m_S23_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR09-016', 'fastq', {'LR09-016':['LR09-016_S24_R1_001.fastq.gz', 'LR09-016_S24_R2_001.fastq.gz'], 'LR09-016f':['LR09-016f_S25_R1_001.fastq.gz', 'LR09-016f_S25_R2_001.fastq.gz'], 'LR09-016m':['LR09-016m_S26_R1_001.fastq.gz', 'LR09-016m_S26_R2_001.fastq.gz']}, 'trio')
# working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0717/b3' #14448
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR17-038', 'fastq', {'LR17-038':['LR17-038_S15_R1_001.fastq.gz', 'LR17-038_S15_R2_001.fastq.gz'], 'LR17-038f':['LR17-038f_S16_R1_001.fastq.gz', 'LR17-038f_S16_R2_001.fastq.gz'], 'LR17-038m':['LR17-038m_S17_R1_001.fastq.gz', 'LR17-038m_S17_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR16-003', 'fastq', {'LR16-003':['LR16-003_S4_R1_001.fastq.gz', 'LR16-003_S4_R2_001.fastq.gz'], 'LR16-003f':['LR16-003f_S5_R1_001.fastq.gz', 'LR16-003f_S5_R2_001.fastq.gz'], 'LR16-003m':['LR16-003m_S6_R1_001.fastq.gz', 'LR16-003m_S6_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR12-188', 'fastq', {'LR12-188':['LR12-188_S37_R1_001.fastq.gz', 'LR12-188_S37_R2_001.fastq.gz'], 'LR12-188f':['LR12-188f_S38_R1_001.fastq.gz', 'LR12-188f_S38_R2_001.fastq.gz'], 'LR12-188m':['LR12-188m_S39_R1_001.fastq.gz', 'LR12-188m_S39_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR11-025', 'fastq', {'LR11-025':['LR11-025_S34_R1_001.fastq.gz', 'LR11-025_S34_R2_001.fastq.gz'], 'LR11-025f':['LR11-025f_S35_R1_001.fastq.gz', 'LR11-025f_S35_R2_001.fastq.gz'], 'LR11-025m':['LR11-025m_S36_R1_001.fastq.gz', 'LR11-025m_S36_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR10-084', 'fastq', {'LR10-084':['LR10-084_S31_R1_001.fastq.gz', 'LR10-084_S31_R2_001.fastq.gz'], 'LR10-084f':['LR10-084f_S32_R1_001.fastq.gz', 'LR10-084f_S32_R2_001.fastq.gz'], 'LR10-084m':['LR10-084m_S33_R1_001.fastq.gz', 'LR10-084m_S33_R2_001.fastq.gz']}, 'trio')
# working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0717/b4'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR03-395', 'fastq', {'LR03-395':['LR03-395_S1_R1_001.fastq.gz', 'LR03-395_S1_R2_001.fastq.gz'], 'LR03-395f':['LR03-395f_S2_R1_001.fastq.gz', 'LR03-395f_S2_R2_001.fastq.gz'], 'LR03-395m':['LR03-395m_S3_R1_001.fastq.gz', 'LR03-395m_S3_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR09-064', 'fastq', {'LR09-064a1':['LR09-064a1_S27_R1_001.fastq.gz', 'LR09-064a1_S27_R2_001.fastq.gz'], 'LR09-064a2':['LR09-064a2_S28_R1_001.fastq.gz', 'LR09-064a2_S28_R2_001.fastq.gz'], 'LR09-064f':['LR09-064f_S29_R1_001.fastq.gz', 'LR09-064f_S29_R2_001.fastq.gz'], 'LR09-064m':['LR09-064m_S30_R1_001.fastq.gz', 'LR09-064m_S30_R2_001.fastq.gz']}, 'trio')
# working_dir = '/home/atimms/ngs_data/exomes/dobyns_genewiz_0717/b5'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR04-239', 'fastq', {'LR04-239':['LR04-239_S7_R1_001.fastq.gz', 'LR04-239_S7_R2_001.fastq.gz'], 'LR04-239s1':['LR04-239s1_S10_R1_001.fastq.gz', 'LR04-239s1_S10_R2_001.fastq.gz'], 'LR04-239f':['LR04-239f_S8_R1_001.fastq.gz', 'LR04-239f_S8_R2_001.fastq.gz'], 'LR04-239m':['LR04-239m_S9_R1_001.fastq.gz', 'LR04-239m_S9_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR12-308', 'fastq', {'LR12-308a1':['LR12-308a1_S11_R1_001.fastq.gz', 'LR12-308a1_S11_R2_001.fastq.gz'], 'LR12-308s1':['LR12-308s1_S14_R1_001.fastq.gz', 'LR12-308s1_S14_R2_001.fastq.gz'], 'LR12-308f':['LR12-308f_S12_R1_001.fastq.gz', 'LR12-308f_S12_R2_001.fastq.gz'], 'LR12-308m':['LR12-308m_S13_R1_001.fastq.gz', 'LR12-308m_S13_R2_001.fastq.gz']}, 'trio')
##mom's sample was dad.. rpt as duo --14469
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR04-290', 'fastq', {'LR04-290':['LR04-290_S8_R1_001.fastq.gz', 'LR04-290_S8_R2_001.fastq.gz'], 'LR04-290f':['LR04-290f_S9_R1_001.fastq.gz', 'LR04-290f_S9_R2_001.fastq.gz']}, 'duo')


##ghayda_ped --14475
working_dir = '/home/atimms/ngs_data/exomes/ghayda_LR17-332_0717'
# dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR17-332', 'fastq', {'LR17-332f':['LG-01_S1_R1_001.fastq.gz', 'LG-01_S1_R2_001.fastq.gz'], 'LR17-332m':['LG-02_S2_R1_001.fastq.gz', 'LG-02_S2_R2_001.fastq.gz'], 'LR17-332a1':['LG-03_S3_R1_001.fastq.gz', 'LG-03_S3_R2_001.fastq.gz'], 'LR17-332a2':['LG-04_S4_R1_001.fastq.gz', 'LG-04_S4_R2_001.fastq.gz'], 'LR17-332m2':['LG-05_S5_R1_001.fastq.gz', 'LG-05_S5_R2_001.fastq.gz'], 'LR17-332a3':['LG-06_S6_R1_001.fastq.gz', 'LG-06_S6_R2_001.fastq.gz'], 'LR17-332a4':['LG-07_S7_R1_001.fastq.gz', 'LG-07_S7_R2_001.fastq.gz']}, 'trio')
# dobyns_exome_pipeline_cybertron_v6.call_just_gemini(working_dir, 'LR17-332', 'fastq', {'LR17-332f':['LG-01_S1_R1_001.fastq.gz', 'LG-01_S1_R2_001.fastq.gz'], 'LR17-332m':['LG-02_S2_R1_001.fastq.gz', 'LG-02_S2_R2_001.fastq.gz'], 'LR17-332a1':['LG-03_S3_R1_001.fastq.gz', 'LG-03_S3_R2_001.fastq.gz'], 'LR17-332a2':['LG-04_S4_R1_001.fastq.gz', 'LG-04_S4_R2_001.fastq.gz'], 'LR17-332m2':['LG-05_S5_R1_001.fastq.gz', 'LG-05_S5_R2_001.fastq.gz'], 'LR17-332a3':['LG-06_S6_R1_001.fastq.gz', 'LG-06_S6_R2_001.fastq.gz'], 'LR17-332a4':['LG-07_S7_R1_001.fastq.gz', 'LG-07_S7_R2_001.fastq.gz']}, 'trio')

##new ped for ghayda
working_dir = '/home/atimms/ngs_data/exomes/ghayda_LR17-076_0817'
dobyns_exome_pipeline_cybertron_v6.call_all_exome_methods_inc_gemini(working_dir, 'LR17-076', 'bam', {'LR17-076':'1715896.bam', 'LR17-076f':'1716714.bam', 'LR17-076m':'1716707.bam'}, 'trio')



