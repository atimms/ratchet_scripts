#!/tools/BioBuilds-2015.04/bin/python
import filtering_annotated
import os
import subprocess
import glob
import re

##set input variables and parameters
delim = '\t'
##working directory
# working_dir = '/data/atimms/cs_rnaseq_1116'
working_dir = '/data/atimms/backedup/cs_rnaseq_1116'
os.chdir(working_dir)

##files etc
project_name = 'cs_rnaseq_1116.c1_and_2'
c1_case_vcf = 'cases_newmerge.vcf.gz'
c1_case_passed_vcf = 'cohort1.cases.passed.vcf.gz'
c1_controls_vcf = 'controls_newmerge.vcf.gz'
c1_controls_passed_vcf = 'cohort1.controls.passed.vcf.gz'
c2_case_vcf = 'css_rnaseq_1016.cases.vcf.gz'
c2_case_passed_vcf = 'css_rnaseq_1016.cases.passed.vcf.gz'
c2_controls_vcf = 'css_rnaseq_1016.controls.vcf.gz'
c2_controls_passed_vcf = 'css_rnaseq_1016.controls.passed.vcf.gz'
exac_vcf = 'ExAC.r0.3.sites.vep.tidy.vcf.gz'
exac_passed_vcf = 'ExAC.r0.3.sites.vep.tidy.passed.vcf.gz'
comb_annotated_cases = project_name + '.cases.annotated.xls'
comb_annotated_ctls = project_name + '.ctls.annotated.xls'
exac_annotated = 'exac.all.annotated.txt'
exac_annotated_passed = 'exac.passed.annotated.txt'

##progrmas
bcftools = '/home/atimms/programs/bcftools-1.2/bcftools'
bgzip = '/home/atimms/programs/samtools-bcftools-htslib-1.0_x64-linux/bin/bgzip'
vt = '/data/atimms/gemini/vt/vt'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
gatk = '/tools/GATK/GenomeAnalysisTK.jar'

##sample info dict
sample_dict = {'95519':['1002', 'F', 'Metopic Synostosis', '1'], '95571':['1003', 'F', 'Metopic Synostosis', '1'], '95611':['1004', 'M', 'Coronal Synostosis R', '1'], 
		'95663':['1006', 'M', 'Metopic Synostosis', '1'], '95639':['1007', 'F', 'Metopic Synostosis', '1'], '95569':['1008', 'M', 'Metopic Synostosis', '1'], 
		'95632':['1009', 'M', 'Lambdoid', '1'], '95510':['1011', 'M', 'Sagittal Synostosis', '1'], '95646':['1012', 'M', 'Sagittal Synostosis', '1'], 
		'95457':['1013', 'F', 'Sagittal Synostosis', '1'], '95598':['1017', 'M', 'Coronal Synostosis R', '1'], '95671':['1018', 'M', 'Metopic Synostosis', '1'], 
		'95585':['1019', 'M', 'Metopic Synostosis', '1'], '95633':['1021', 'F', 'Coronal Synostosis L', '1'], '95575':['1022', 'F', 'Sagittal Synostosis', '1'], 
		'95705':['1024', 'M', 'Sagittal Synostosis', '1'], '95649':['1026', 'M', 'Metopic Synostosis', '1'], '95578':['1027', 'F', 'Sagittal Synostosis', '1'], 
		'95616':['1028', 'F', 'Coronal Synostosis L', '1'], '95512':['1031', 'M', 'Sagittal Synostosis', '1'], '95654':['1032', 'M', 'Sagittal Synostosis', '1'], 
		'95543':['1034', 'M', 'Metopic Synostosis', '1'], '95594':['1036', 'F', 'Metopic Synostosis', '1'], '95624':['1037', 'F', 'Sagittal Synostosis', '1'], 
		'95623':['1038', 'M', 'Sagittal Synostosis', '1'], '95475':['1040', 'M', 'Lambdoid', '1'], '95644':['1041', 'F', 'Coronal Synostosis L', '1'], 
		'95601':['1044', 'F', 'Coronal Synostosis R', '1'], '95562':['1047', 'M', 'Sagittal Synostosis', '1'], '95620':['1048', 'F', 'Coronal Synostosis L', '1'], 
		'95447':['1050', 'F', 'Sagittal Synostosis', '1'], '95570':['1052', 'M', 'Sagittal Synostosis', '1'], '95603':['1053', 'F', 'Coronal Synostosis L', '1'], 
		'95530':['1056', 'M', 'Sagittal Synostosis', '1'], '95481':['1060', 'M', 'Lambdoid', '1'], '95514':['1061', 'F', 'Coronal Synostosis R', '1'], 
		'95470':['1062', 'F', 'Sagittal Synostosis', '1'], '95501':['1063', 'M', 'Sagittal Synostosis', '1'], '95531':['1065', 'M', 'Metopic Synostosis', '1'], 
		'95550':['1066', 'F', 'Coronal Synostosis R', '1'], '95686':['1067', 'F', 'Coronal Synostosis L', '1'], '95693':['1068', 'M', 'Sagittal Synostosis', '1'], 
		'95628':['1070', 'M', 'Metopic Synostosis', '1'], '95683':['1072', 'F', 'Coronal Synostosis R', '1'], '95460':['1074', 'F', 'Metopic Synostosis', '1'], 
		'95494':['1076', 'M', 'Sagittal Synostosis', '1'], '95498':['1078', 'M', 'Sagittal Synostosis', '1'], '95452':['1079', 'F', 'Metopic Synostosis', '1'], 
		'95576':['1082', 'M', 'Metopic Synostosis', '1'], '95567':['1084', 'M', 'Coronal Synostosis L', '1'], '95661':['1085', 'F', 'Coronal Synostosis L', '1'], 
		'95655':['1088', 'M', 'Metopic Synostosis', '1'], '95592':['1089', 'M', 'Metopic Synostosis', '1'], '95657':['1090', 'M', 'Coronal Synostosis L', '1'], 
		'95680':['1091', 'M', 'Sagittal Synostosis', '1'], '95629':['1092', 'F', 'Sagittal Synostosis', '1'], '95526':['1095', 'M', 'Metopic Synostosis', '1'], 
		'95492':['1096', 'M', 'Metopic Synostosis', '1'], '95496':['1097', 'F', 'Coronal Synostosis L', '1'], '95593':['1098', 'M', 'Coronal Synostosis R', '1'], 
		'95506':['2001', 'F', 'Coronal Synostosis L', '1'], '95681':['2003', 'F', 'Coronal Synostosis L', '1'], '95551':['2004', 'F', 'Lambdoid', '1'], 
		'95615':['2005', 'M', 'Metopic Synostosis', '1'], '95617':['2006', 'M', 'Metopic Synostosis', '1'], '95586':['2007', 'M', 'Sagittal Synostosis', '1'], 
		'95631':['2008', 'M', 'Metopic Synostosis', '1'], '95525':['2009', 'M', 'Metopic Synostosis', '1'], '95565':['2010', 'F', 'Metopic Synostosis', '1'], 
		'95609':['2012', 'M', 'Coronal Synostosis R', '1'], '95566':['2013', 'F', 'Lambdoid', '1'], '95516':['2014', 'F', 'Coronal Synostosis L', '1'], 
		'95658':['2015', 'M', 'Sagittal Synostosis', '1'], '95684':['2023', 'M', 'Sagittal Synostosis', '1'], '95691':['2024', 'M', 'Metopic Synostosis', '1'], 
		'95484':['2026', 'M', 'Sagittal Synostosis', '1'], '95507':['2027', 'F', 'Sagittal Synostosis', '1'], '95463':['2028', 'M', 'Lambdoid', '1'], 
		'95698':['2030', 'M', 'Sagittal Synostosis', '1'], '95668':['2031', 'M', 'Lambdoid', '1'], '95518':['2032', 'F', 'Metopic Synostosis', '1'], 
		'95468':['2034', 'F', 'Coronal Synostosis R', '1'], '95606':['2035', 'F', 'Coronal Synostosis R', '1'], '95511':['2037', 'M', 'Sagittal Synostosis', '1'], 
		'95573':['2038', 'F', 'Sagittal Synostosis', '1'], '95653':['2040', 'M', 'Sagittal Synostosis', '1'], '95638':['2041', 'M', 'Metopic Synostosis', '1'], 
		'95679':['2042', 'M', 'Sagittal Synostosis', '1'], '95488':['2045', 'M', 'Metopic Synostosis', '1'], '95703':['2048', 'M', 'Sagittal Synostosis', '1'], 
		'95597':['2049', 'F', 'Coronal Synostosis L', '1'], '95509':['2050', 'F', 'Coronal Synostosis L', '1'], '95677':['2052', 'M', 'Coronal Synostosis L', '1'], 
		'95545':['2059', 'M', 'Sagittal Synostosis', '1'], '95635':['2060', 'M', 'Sagittal Synostosis', '1'], '95482':['2062', 'F', 'Coronal Synostosis R', '1'], 
		'95558':['2063', 'M', 'Coronal Synostosis R', '1'], '95479':['2064', 'M', 'Sagittal Synostosis', '1'], '95540':['2066', 'M', 'Coronal Synostosis L', '1'], 
		'95634':['2068', 'M', 'Metopic Synostosis', '1'], '95472':['2070', 'F', 'Sagittal Synostosis', '1'], '95690':['2071', 'F', 'Coronal Synostosis R', '1'], 
		'95700':['2073', 'M', 'Sagittal Synostosis', '1'], '95466':['2074', 'F', 'Sagittal Synostosis', '1'], '95596':['2079', 'F', 'Coronal Synostosis L', '1'], 
		'95552':['2083', 'M', 'Coronal Synostosis L', '1'], '95640':['2084', 'F', 'Coronal Synostosis L', '1'], '95699':['2085', 'M', 'Metopic Synostosis', '1'], 
		'95445':['2086', 'F', 'Coronal Synostosis R', '1'], '95689':['2088', 'M', 'Metopic Synostosis', '1'], '95527':['2089', 'F', 'Coronal Synostosis R', '1'], 
		'95500':['2090', 'M', 'Sagittal Synostosis', '1'], '95450':['2092', 'M', 'Coronal Synostosis R', '1'], '95485':['2094', 'F', 'Coronal Synostosis R', '1'], 
		'95453':['2095', 'M', 'Coronal Synostosis R', '1'], '95591':['2101', 'F', 'Metopic Synostosis', '1'], '95619':['3001', 'M', 'Sagittal Synostosis', '1'], 
		'95490':['3003', 'F', 'Coronal Synostosis R', '1'], '95583':['3005', 'F', 'Sagittal Synostosis', '1'], '95676':['3006', 'M', 'Lambdoid', '1'], 
		'95448':['3007', 'F', 'Coronal Synostosis R', '1'], '95604':['3008', 'M', 'Sagittal Synostosis', '1'], '95522':['3010', 'M', 'Sagittal Synostosis', '1'], 
		'95602':['3011', 'M', 'Sagittal Synostosis', '1'], '95542':['3012', 'M', 'Lambdoid', '1'], '95454':['3013', 'M', 'Sagittal Synostosis', '1'], 
		'95547':['3014', 'F', 'Coronal Synostosis L', '1'], '95694':['3015', 'F', 'Sagittal Synostosis', '1'], '95696':['3016', 'M', 'Sagittal Synostosis', '1'], 
		'95579':['3018', 'M', 'Sagittal Synostosis', '1'], '95537':['3019', 'M', 'Sagittal Synostosis', '1'], '95471':['3020', 'M', 'Sagittal Synostosis', '1'], 
		'95535':['3021', 'M', 'Sagittal Synostosis', '1'], '95627':['3022', 'F', 'Sagittal Synostosis', '1'], '95622':['3023', 'M', 'Metopic Synostosis', '1'], 
		'95524':['3025', 'M', 'Sagittal Synostosis', '1'], '95641':['3026', 'M', 'Coronal Synostosis R', '1'], '95589':['3027', 'M', 'Coronal Synostosis R', '1'], 
		'95560':['3028', 'M', 'Sagittal Synostosis', '1'], '95652':['3029', 'M', 'Sagittal Synostosis', '1'], '95491':['3030', 'M', 'Coronal Synostosis R', '1'], 
		'95546':['3031', 'F', 'Sagittal Synostosis', '1'], '95544':['3032', 'M', 'Sagittal Synostosis', '1'], '95697':['3033', 'M', 'Sagittal Synostosis', '1'], 
		'95474':['3034', 'M', 'Sagittal Synostosis', '1'], '95614':['3035', 'M', 'Metopic Synostosis', '1'], '95487':['4001', 'M', 'Sagittal Synostosis', '1'], 
		'95515':['4002', 'F', 'Sagittal Synostosis', '1'], '95656':['4003', 'F', 'Metopic Synostosis', '1'], '95580':['4005', 'F', 'Sagittal Synostosis', '1'], 
		'95505':['4006', 'M', 'Sagittal Synostosis', '1'], '95476':['4007', 'M', 'Metopic Synostosis', '1'], '95600':['4009', 'M', 'Sagittal Synostosis', '1'], 
		'95520':['4010', 'M', 'Sagittal Synostosis', '1'], '95477':['4011', 'M', 'Sagittal Synostosis', '1'], '95536':['4013', 'M', 'Sagittal Synostosis', '1'], 
		'95665':['4014', 'M', 'Sagittal Synostosis', '1'], '95446':['4015', 'F', 'Sagittal Synostosis', '1'], '95587':['4016', 'M', 'Sagittal Synostosis', '1'], 
		'95503':['4017', 'M', 'Sagittal Synostosis', '1'], '95584':['4018', 'F', 'Metopic Synostosis', '1'], '95662':['4019', 'M', 'Metopic Synostosis', '1'], 
		'95678':['4020', 'M', 'Metopic Synostosis', '1'], '95502':['4022', 'M', 'Coronal Synostosis R', '1'], '95555':['4025', 'M', 'Sagittal Synostosis', '1'], 
		'95534':['4027', 'F', 'Sagittal Synostosis', '1'], '95642':['4032', 'F', 'Metopic Synostosis', '1'], '95549':['4033', 'F', 'Sagittal Synostosis', '1'], 
		'95636':['4035', 'M', 'Metopic Synostosis', '1'], '95469':['4037', 'M', 'Metopic Synostosis', '1'], '95556':['4038', 'F', 'Metopic Synostosis', '1'], 
		'95595':['4039', 'M', 'Sagittal Synostosis', '1'], '95702':['4040', 'F', 'Sagittal Synostosis', '1'], '95673':['1059/C1177', 'M', 'Sagittal Synostosis', '1'], 
		'95669':['1071/C1290', 'M', 'Sagittal Synostosis', '1'], '95685':['1075/C1269', 'M', 'Sagittal Synostosis', '1'], '95664':['1081/C1183', 'F', 'Sagittal Synostosis', '1'], 
		'95532':['AUT01', 'M', 'Control', '1'], '95541':['AUT05', 'M', 'Control', '1'], '95590':['AUT06', 'M', 'Control', '1'], '95667':['AUT07', 'M', 'Control', '1'], 
		'95660':['AUT10', 'M', 'Control', '1'], '95704':['AUT13', 'F', 'Control', '1'], '95574':['AUT14', 'M', 'Control', '1'], '95499':['AUT19', 'M', 'Control', '1'], 
		'95521':['C1001', 'F', 'Coronal Synostosis R', '1'], '95688':['C1004', 'M', 'Sagittal Synostosis', '1'], '95493':['C1053', 'F', 'Coronal Synostosis R', '1'], 
		'95651':['C1080', 'M', 'Sagittal Synostosis', '1'], '95458':['C1116', 'M', 'Sagittal Synostosis', '1'], '95701':['C1125/SS14OC', 'M', 'Sagittal Synostosis', '1'], 
		'95557':['C1162', 'M', 'Metopic Synostosis', '1'], '95581':['C1174', 'M', 'Sagittal Synostosis', '1'], '95572':['C1180', 'M', 'Metopic Synostosis', '1'], 
		'95605':['C1186', 'M', 'Metopic Synostosis', '1'], '95486':['C1192', 'F', 'Sagittal Synostosis', '1'], '95670':['C1199', 'M', 'Sagittal Synostosis', '1'], 
		'95650':['C1202', 'M', 'Sagittal Synostosis', '1'], '95612':['C1205', 'M', 'Lambdoid', '1'], '95626':['C1214', 'M', 'Metopic Synostosis', '1'], 
		'95462':['C1221', 'M', 'Sagittal Synostosis', '1'], '95610':['C1253', 'M', 'Sagittal Synostosis', '1'], '95564':['C1258', 'M', 'Sagittal Synostosis', '1'], 
		'95675':['C1272', 'M', 'Sagittal Synostosis', '1'], '95464':['C1275', 'M', 'Sagittal Synostosis', '1'], '95588':['C1284', 'M', 'Sagittal Synostosis', '1'], 
		'95528':['C1296', 'M', 'Sagittal Synostosis', '1'], '95495':['C1299', 'M', 'Lambdoid', '1'], '95529':['C1302', 'M', 'Sagittal Synostosis', '1'], 
		'95456':['C1308', 'M', 'Sagittal Synostosis', '1'], '95630':['C1311', 'M', 'Sagittal Synostosis', '1'], '95559':['C1317', 'M', 'Metopic Synostosis', '1'], 
		'95467':['C1330', 'M', 'Sagittal Synostosis', '1'], '95682':['C1336', 'F', 'Coronal Synostosis L', '1'], '95637':['C1365', 'M', 'Lambdoid', '1'], 
		'95563':['C1374', 'M', 'Metopic Synostosis', '1'], '95599':['C1388', 'M', 'Coronal Synostosis R', '1'], '95459':['COR02', 'M', 'Coronal Synostosis R', '1'], 
		'95643':['OST01', 'M', 'Control', '1'], '95561':['OST02', 'M', 'Control', '1'], '95554':['OST03', 'M', 'Control', '1'], '95692':['OST04', 'M', 'Control', '1'], 
		'95523':['OST06', 'M', 'Control', '1'], '95517':['OST07', 'F', 'Control', '1'], '95687':['OST09', 'F', 'Control', '1'], '95674':['OST15', 'M', 'Control', '1'], 
		'95647':['OST17', 'M', 'Control', '1'], '95497':['OST18', 'F', 'Control', '1'], '95577':['OST19', 'F', 'Control', '1'], '95508':['OST20', 'M', 'Control', '1'], 
		'95607':['OST22', 'M', 'Control', '1'], '95480':['OST25', 'F', 'Control', '1'], '95648':['OST26', 'F', 'Control', '1'], '95582':['OST29', 'F', 'Control', '1'], 
		'95621':['OST31', 'M', 'Control', '1'], '95659':['OST32', 'M', 'Control', '1'], '95461':['OST33', 'M', 'Control', '1'], '95695':['OST35', 'M', 'Control', '1'], 
		'95538':['OST36', 'F', 'Control', '1'], '95451':['OST37', 'F', 'Control', '1'], '95568':['OST41', 'F', 'Control', '1'], '95483':['OST42', 'M', 'Control', '1'], 
		'95618':['OST43', 'M', 'Control', '1'], '95504':['OST47', 'F', 'Control', '1'], '95455':['OST49', 'M', 'Control', '1'], '95666':['OST50', 'F', 'Control', '1'], 
		'95608':['OST51', 'M', 'Control', '1'], '95465':['OST52', 'M', 'Control', '1'], '95613':['OST54', 'F', 'Control', '1'], '95473':['OST55', 'M', 'Control', '1'], 
		'95539':['OST66', 'M', 'Control', '1'], '95533':['OST67', 'M', 'Control', '1'], '95489':['OST68', 'M', 'Control', '1'], '95645':['SS-10-BP', 'M', 'Sagittal Synostosis', '1'], 
		'95672':['STL04', 'M', 'Control', '1'], '95513':['STL05', 'M', 'Control', '1'], '95625':['STL07', 'M', 'Control', '1'], '95449':['STL09', 'F', 'Control', '1'], 
		'95548':['STL11', 'M', 'Control', '1'], '95478':['STL12', 'M', 'Control', '1'], '95553':['STL15', 'M', 'Control', '1'], '163499':['C3010', 'M', 'Sagittal', '2'], 
		'163500':['AUT30', 'M', 'Control', '2'], '163501':['OST78', 'F', 'Control', '2'], '163502':['C1489', 'F', 'Coronal Synostosis-R', '2'], '163503':['C1600', 'M', 'Metopic', '2'], 
		'163504':['C3027', 'F', 'Coronal', '2'], '163505':['C1698', 'M', 'Coronal Synostosis-R', '2'], '163506':['C1656', 'F', 'Sagittal', '2'], '163507':['C3024', 'M', 'Sagittal', '2'], 
		'163508':['C1889', 'F', 'Coronal Synostosis-L', '2'], '163509':['C1861', 'F', 'Coronal Synostosis-R', '2'], '163510':['C1864', 'F', 'Sagittal', '2'], '163511':['C1763', 'M', 'Metopic', '2'], 
		'163512':['STL24', 'F', 'Control', '2'], '163513':['C1671', 'M', 'Metopic', '2'], '163514':['STL19', 'M', 'Control', '2'], '163515':['STL23', 'F', 'Control', '2'], 
		'163516':['C1638', 'F', 'Metopic', '2'], '163517':['C3035', 'M', 'Metopic', '2'], '163518':['C1841', 'M', 'Sagittal', '2'], '163519':['C1683', 'M', 'Sagittal', '2'], 
		'163520':['C1635', 'M', 'Sagittal', '2'], '163521':['C3047', 'F', 'Metopic', '2'], '163522':['C3040', 'M', 'Sagittal', '2'], '163523':['C1806', 'M', 'Sagittal', '2'], 
		'163524':['C3016', 'M', 'Sagittal', '2'], '163525':['C1707', 'F', 'Coronal Synostosis-L', '2'], '163526':['C3034', 'F', 'Coronal', '2'], '163527':['C1785', 'M', 'Sagittal', '2'], 
		'163528':['C1597', 'M', 'Sagittal', '2'], '163529':['C3006', 'F', 'Metopic', '2'], '163530':['C3054', 'F', 'Coronal', '2'], '163531':['C3012', 'M', 'Metopic', '2'], 
		'163532':['C1926', 'M', 'Sagittal', '2'], '163533':['C1757', 'M', 'Sagittal', '2'], '163534':['C3020', 'F', 'Coronal', '2'], '163535':['C1766', 'M', 'Sagittal', '2'], 
		'163536':['STL44', 'M', 'Control', '2'], '163537':['C1686', 'M', 'Sagittal', '2'], '163538':['C3002', 'M', 'Metopic', '2'], '163539':['C1793', 'M', 'Sagittal', '2'], 
		'163540':['C1760', 'M', 'Sagittal', '2'], '163541':['C1586', 'F', 'Sagittal', '2'], '163542':['AUT31', 'F', 'Control', '2'], '163543':['C3021', 'M', 'Sagittal', '2'], 
		'163544':['C1923', 'M', 'Metopic', '2'], '163545':['C1769', 'M', 'Sagittal', '2'], '163546':['C1704', 'M', 'Metopic', '2'], '163547':['C1551', 'M', 'Sagittal', '2'], 
		'163548':['C1915', 'M', 'Coronal Synostosis-R', '2'], '163549':['C1830', 'M', 'Metopic', '2'], '163550':['AUT36', 'F', 'Control', '2'], '163551':['STL39', 'M', 'Control', '2'], 
		'163552':['STL35', 'M', 'Control', '2'], '163553':['C3044', 'M', 'Metopic', '2'], '163554':['AUT25', 'M', 'Control', '2'], '163555':['C1612', 'F', 'Sagittal', '2'], 
		'163556':['C1831', 'M', 'Metopic', '2'], '163557':['C1616', 'F', 'Coronal Synostosis-L', '2'], '163558':['C1790', 'F', 'Lambdoid-R', '2'], '163559':['C1692', 'M', 'Coronal Synostosis-R', '2'], 
		'163560':['C3049', 'F', 'Sagittal', '2'], '163561':['C1827', 'F', 'Coronal Synostosis-R', '2'], '163562':['C1695', 'M', 'Coronal Synostosis-R', '2'], '163563':['C1812', 'F', 'Coronal Synostosis-L', '2'], 
		'163564':['C1856', 'M', 'Metopic', '2'], '163565':['C1781', 'F', 'Sagittal', '2'], '163566':['C1426', 'M', 'Metopic', '2'], '163567':['C1594', 'M', 'Sagittal', '2'], '163568':['C1665', 'F', 'Sagittal', '2'], 
		'163569':['C3059', 'M', 'Metopic', '2'], '163570':['OST79', 'M', 'Control', '2'], '163571':['C1580', 'M', 'Sagittal', '2'], '163572':['C3057', 'M', 'Sagittal', '2'], '163573':['C3036', 'F', 'Sagittal', '2'], 
		'163574':['AUT35', 'F', 'Control', '2'], '163575':['C1589', 'F', 'Metopic', '2'], '163576':['C1619', 'M', 'Sagittal', '2'], '163577':['C1701', 'M', 'Metopic', '2'], '163578':['C3014', 'M', 'Sagittal', '2'], 
		'163579':['C1803', 'M', 'Sagittal', '2'], '163580':['C1853', 'M', 'Metopic', '2'], '163581':['C3032', 'M', 'Metopic', '2'], '163582':['C1483', 'F', 'Coronal Synostosis-R', '2'], 
		'163583':['STL37', 'M', 'Control', '2'], '163584':['C3046', 'F', 'Sagittal', '2'], '163585':['C3045', 'M', 'Sagittal', '2'], '163586':['C3060', 'F', 'Coronal', '2'], 
		'163587':['C1632', 'M', 'Lambdoid-L', '2'], '163588':['C1725', 'F', 'Metopic', '2'], '163589':['STL36', 'F', 'Control', '2'], '163590':['C1577', 'M', 'Metopic', '2'], 
		'163591':['C3007', 'M', 'Sagittal', '2'], '163592':['C3056', 'M', 'Sagittal', '2'], '163593':['C3055', 'M', 'Metopic', '2'], '163594':['C1647', 'M', 'Coronal Synostosis-R', '2'], 
		'163595':['C3038', 'M', 'Sagittal', '2'], '163596':['C1627', 'M', 'Sagittal', '2'], '163597':['C1486', 'F', 'Sagittal', '2'], '163598':['C1570', 'M', 'Coronal Synostosis-L', '2'], 
		'163599':['STL20', 'M', 'Control', '2'], '163600':['STL41', 'M', 'Control', '2'], '163601':['C1787', 'F', 'Sagittal', '2'], '163602':['C1453', 'F', 'Sagittal', '2'], 
		'163603':['C1751', 'F', 'Coronal Synostosis-R', '2'], '163604':['C1859', 'M', 'Metopic', '2'], '163605':['C1833', 'M', 'Sagittal', '2'], '163606':['C1722', 'M', 'Sagittal', '2'], 
		'163607':['C3043', 'M', 'Sagittal', '2'], '163608':['STL45', 'M', 'Control', '2'], '163609':['C3050', 'M', 'Lambdoid', '2'], '163610':['C1531', 'M', 'Sagittal', '2'], 
		'163611':['C1782', 'M', 'Sagittal', '2'], '163612':['AUT34', 'M', 'Control', '2'], '163613':['STL22', 'M', 'Control', '2'], '163614':['C3019', 'M', 'Lambdoid', '2'], 
		'163615':['C1583', 'M', 'Sagittal', '2'], '163616':['C1713', 'F', 'Coronal Synostosis-L', '2'], '163617':['STL32', 'M', 'Control', '2'], '163618':['OST76', 'M', 'Control', '2'], 
		'163619':['C1456', 'M', 'Metopic', '2'], '163620':['C3022', 'M', 'Metopic', '2'], '163621':['C1950', 'M', 'Sagittal', '2'], '163622':['C1775', 'F', 'Sagittal', '2'], 
		'163623':['C3026', 'M', 'Coronal', '2'], '163624':['C1406', 'F', 'Sagittal', '2'], '163625':['C1641', 'M', 'Sagittal', '2'], '163626':['C3003', 'M', 'Metopic', '2'], 
		'163627':['C1650', 'M', 'Sagittal', '2'], '163628':['C1606', 'F', 'Coronal Synostosis-R', '2'], '163629':['C1609', 'M', 'Metopic', '2'], '163630':['C3028', 'M', 'Sagittal', '2'], 
		'163631':['C3039', 'F', 'Lambdoid', '2'], '163632':['C3011', 'M', 'Sagittal', '2'], '163633':['C3031', 'M', 'Coronal ', '2'], '163634':['C1576', 'F', 'Metopic', '2'], 
		'163635':['C1622', 'M', 'Lambdoid-R', '2'], '163636':['C1653', 'M', 'Sagittal', '2'], '163637':['C3017', 'M', 'Metopic', '2'], '163638':['STL46', 'M', 'Control', '2'], 
		'163639':['C1970', 'M', 'Sagittal', '2'], '163640':['C1710', 'M', 'Sagittal', '2'], '163641':['C1625', 'M', 'Sagittal', '2'], '163642':['C1591', 'F', 'Sagittal', '2'], 
		'163643':['C1603', 'M', 'Metopic', '2'], '163644':['C1850', 'F', 'Sagittal', '2'], '163645':['C2063', 'M', 'Sagittal', '2'], '163646':['C2038', 'F', 'Coronal Synostosis-R', '2'], 
		'163647':['C1869', 'M', 'Metopic', '2'], '163648':['STL17', 'F', 'Control', '2'], '163649':['C1957', 'F', 'Sagittal', '2'], '163650':['AUT28', 'M', 'Control', '2'], 
		'163651':['C3018', 'M', 'Sagittal', '2'], '163652':['C3023', 'F', 'Metopic', '2'], '163653':['C3041', 'M', 'Metopic', '2'], '163654':['C1466', 'M', 'Sagittal', '2'], 
		'163655':['C2030', 'F', 'Coronal Synostosis-R', '2'], '163656':['C1662', 'M', 'Metopic', '2'], '163657':['C3013', 'F', 'Sagittal', '2'], '163658':['C1459', 'F', 'Coronal Synostosis-R', '2'], 
		'163659':['C3064', 'M', 'Sagittal', '2'], '163660':['C1534', 'M', 'Sagittal', '2'], '163661':['C1875', 'F', 'Sagittal', '2'], '163662':['C1954', 'F', 'Coronal Synostosis-R', '2'], 
		'163663':['C2017', 'F', 'Coronal Synostosis-R', '2'], '163664':['C2033', 'M', 'Sagittal', '2'], '163665':['C3062', 'F', 'Sagittal', '2'], '163666':['C1680', 'F', 'Sagittal', '2'], 
		'163667':['C1432', 'M', 'Sagittal', '2'], '163668':['C3061', 'F', 'Sagittal', '2'], '163669':['C1677', 'M', 'Sagittal', '2'], '163670':['STL60', 'M', 'Control', '2'], '163671':['STL64', 'M', 'Control', '2'], 
		'163672':['C3042', 'M', 'Sagittal', '2'], '163673':['STL51', 'F', 'Control', '2'], '163674':['C1728', 'M', 'Sagittal', '2'], '163675':['STL59', 'M', 'Control', '2'], '163676':['AUT20', 'M', 'Control', '2'], 
		'163677':['STL28', 'M', 'Control', '2'], '163678':['C3051', 'M', 'Metopic', '2'], '163679':['C2035', 'F', 'Sagittal', '2'], '163680':['C1847', 'M', 'Metopic', '2'], '163681':['C1573', 'M', 'Metopic', '2'], 
		'163682':['C3033', 'M', 'Sagittal', '2'], '163683':['C1778', 'M', 'Sagittal', '2'], '163684':['C3005', 'M', 'Sagittal', '2'], '163685':['C3015', 'M', 'Metopic', '2'], '163686':['C1844', 'F', 'Lambdoid-L', '2'], 
		'163687':['C3000', 'M', 'Metopic', '2'], '163688':['C1754', 'M', 'Sagittal', '2'], '163689':['STL21', 'M', 'Control', '2'], '163690':['C3029', 'M', 'Metopic', '2'], '163691':['C1742', 'M', 'Sagittal', '2'], 
		'163692':['C1736', 'M', 'Sagittal', '2'], '163693':['AUT32', 'F', 'Control', '2'], '163694':['C1809', 'M', 'Metopic', '2'], '163695':['C3009', 'M', 'Sagittal', '2'], '163696':['C1892', 'M', 'Sagittal', '2'], 
		'163697':['C2060', 'M', 'Sagittal', '2'], '163698':['C1878', 'M', 'Lambdoid-R', '2'], '163699':['C3037', 'M', 'Sagittal', '2'], '163700':['C1674', 'F', 'Coronal Synostosis-R', '2'], 
		'163701':['C1668', 'M', 'Sagittal', '2'], '163702':['C3025', 'F', 'Sagittal', '2'], '163703':['C3053', 'M', 'Sagittal', '2'], '163704':['C1733', 'F', 'Sagittal', '2'], '163705':['STL54', 'F', 'Control', '2'], 
		'163706':['C1689', 'M', 'Sagittal', '2'], '163707':['C1799', 'M', 'Sagittal', '2'], '163708':['STL58', 'M', 'Control', '2'], '163709':['C1403', 'M', 'Sagittal', '2'], '163710':['C3048', 'F', 'Sagittal', '2'], 
		'163711':['C1748', 'M', 'Sagittal', '2'], '163712':['C1796', 'F', 'Lambdoid-L', '2'], '163713':['C1815', 'M', 'Sagittal', '2'], '163714':['C2084', 'M', 'Sagittal', '2'], '163715':['OST82', 'F', 'Control', '2'], 
		'163716':['OST85', 'M', 'Control', '2'], '163717':['C2072', 'M', 'Sagittal', '2'], '163718':['C2088', 'M', 'Metopic', '2'], '163719':['C2080', 'M', 'Sagittal', '2'], '163720':['C3065', 'M', 'Metopic', '2'], 
		'163721':['C3066', 'F', 'Coronal', '2']}


#annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/data/atimms/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,dbnsfp30a,popfreq_all_20150413,kaviar_20150923,hrcr1,cosmic70,dbscsnv11,dbnsfp31a_interpro,avsnp147']
av_operation = ['-operation', 'g,r,r,f,f,f,f,f,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput']

##variables for filtering

##definitions for determining var function
exonic_definitions = ['exonic', 'splicing']
disruptive_definitions = ['','.', 'stopgain', 'stoploss','frameshift insertion', 'frameshift deletion', 'frameshift_insertion', 'frameshift_deletion']
protein_changing_definitions =['','.', 'stopgain', 'stoploss', 'nonsynonymous SNV','frameshift insertion', 'frameshift deletion', 'nonsynonymous_SNV','frameshift_insertion', 'frameshift_deletion']
nosyn_definitions =['nonsynonymous SNV','nonsynonymous_SNV']
pp2_score = [0.9, 0.9]
cadd_gerp_score = [15,3]
damaging_definitions = [0.9,0.9,15,3]
no_record_definition = ['', '.']
##list of columns for both methods (all start at 1):
##for filtering variants: refgene function, frequency dbs, refgene exonic function, pp2 x2/cadd/gerp, refgene geneame, quality and coverage
##using max_aaf
# filtering_cols = [15, 56, 18, [24,26,40,50], 16, [7,8]]
##using eac all
filtering_cols = [15, 63, 18, [24,26,40,50], 16, [7,8]]

frequencies = [0.01, 0.02, 0.04, 0.001]



##methods
##simple list maker -- give string and dictionary to get length from
def make_list(string, dict):
	l = []
	for i in range(len(dict)):
		l.append(string)
	return l



def filter_vcf_file(in_vcf, out_vcf):
	bcftools_filter = subprocess.Popen([bcftools, 'view', '-f', 'PASS', '-o', out_vcf, '-O', 'z', in_vcf])
	bcftools_filter.wait()

##convert vcf file to individual avinput file
def convert_to_annovar(vcf, project_prefix):
	av_prefix = project_prefix
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', av_prefix])
	con_ann.wait()

##annotate vcf file
def table_annovar_vcf(vcf, project_prefix):
		out_prefix = project_prefix
		command = [table_annovar] + av_buildver + [vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', out_prefix]
		annovar = subprocess.Popen(command)
		annovar.wait()

def split_info_field(info_list):
	indices = [i for i, s in enumerate(info_list) if 'ANNOVAR_DATE' in s]
	# print indices
	i_count = 0
	final_list = []
	for info in info_list:
		# print info
		if i_count > indices[0] and info != 'ALLELE_END':
			info2 = info.split('=')[1]
			#print info2
			final_list.append(info2)
		i_count += 1
	return final_list
	
def format_avinput(project_prefix):
	avinputs = glob.glob(project_prefix + '.*.avinput')
	print avinputs
	for avinput in avinputs:
		head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Zygosity', 'Qual', 'Coverage', 'Filter', 'Pos', 'Ref2', 'Alt2', 'Format', 'Info', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 'genomicSuperDups', 'SIFT_score', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score', 'Polyphen2_HVAR_pred', 'LRT_score', 'LRT_pred', 'MutationTaster_score', 'MutationTaster_pred', 'MutationAssessor_score', 'MutationAssessor_pred', 'FATHMM_score', 'FATHMM_pred', 'PROVEAN_score', 'PROVEAN_pred', 'VEST3_score', 'CADD_raw', 'CADD_phred', 'DANN_score', 'fathmm-MKL_coding_score', 'fathmm-MKL_coding_pred', 'MetaSVM_score', 'MetaSVM_pred', 'MetaLR_score', 'MetaLR_pred', 'integrated_fitCons_score', 'integrated_confidence_value', 'GERP++_RS', 'phyloP7way_vertebrate', 'phyloP20way_mammalian', 'phastCons7way_vertebrate', 'phastCons20way_mammalian', 'SiPhy_29way_logOdds', 'PopFreqMax', '1000G_ALL', '1000G_AFR', '1000G_AMR', '1000G_EAS', '1000G_EUR', '1000G_SAS', 'ExAC_ALL', 'ExAC_AFR', 'ExAC_AMR', 'ExAC_EAS', 'ExAC_FIN', 'ExAC_NFE', 'ExAC_OTH', 'ExAC_SAS', 'ESP6500siv2_ALL', 'ESP6500siv2_AA', 'ESP6500siv2_EA', 'CG46', 'Kaviar_AF', 'Kaviar_AC', 'Kaviar_AN', 'HRC_AF', 'HRC_AC', 'HRC_AN', 'HRC_non1000G_AF', 'HRC_non1000G_AC', 'HRC_non1000G_AN', 'cosmic70', 'dbscSNV_ADA_SCORE', 'dbscSNV_RF_SCORE', 'Interpro_domain', 'avsnp147']
		head_out = delim.join(head + ['\n'])
		sample = avinput.split('.')[2]
		outfile = project_prefix + '.' + sample + '.annotated.txt'
		with open(avinput, "r") as av, open(outfile, "w") as final:
			final.write(head_out)
			for line in av:
				line = line.strip('\n').split(delim)
				stuff = line[0:8] + [line[14]] + [line[9]] + line[11:13]
				info = line[15].split(';')
				info_list = split_info_field(info)
				other_stuff = line[16:]
				line_out = delim.join(stuff + other_stuff + info_list +['\n'])
				final.write(line_out)

def format_avinput_exac(project_prefix):
	avinput = project_prefix
	print project_prefix
	head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Zygosity', 'Qual', 'Coverage', 'Filter', 'Pos', 'Ref2', 'Alt2', 'Format', 'Info', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 'genomicSuperDups', 'SIFT_score', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score', 'Polyphen2_HVAR_pred', 'LRT_score', 'LRT_pred', 'MutationTaster_score', 'MutationTaster_pred', 'MutationAssessor_score', 'MutationAssessor_pred', 'FATHMM_score', 'FATHMM_pred', 'PROVEAN_score', 'PROVEAN_pred', 'VEST3_score', 'CADD_raw', 'CADD_phred', 'DANN_score', 'fathmm-MKL_coding_score', 'fathmm-MKL_coding_pred', 'MetaSVM_score', 'MetaSVM_pred', 'MetaLR_score', 'MetaLR_pred', 'integrated_fitCons_score', 'integrated_confidence_value', 'GERP++_RS', 'phyloP7way_vertebrate', 'phyloP20way_mammalian', 'phastCons7way_vertebrate', 'phastCons20way_mammalian', 'SiPhy_29way_logOdds', 'PopFreqMax', '1000G_ALL', '1000G_AFR', '1000G_AMR', '1000G_EAS', '1000G_EUR', '1000G_SAS', 'ExAC_ALL', 'ExAC_AFR', 'ExAC_AMR', 'ExAC_EAS', 'ExAC_FIN', 'ExAC_NFE', 'ExAC_OTH', 'ExAC_SAS', 'ESP6500siv2_ALL', 'ESP6500siv2_AA', 'ESP6500siv2_EA', 'CG46', 'Kaviar_AF', 'Kaviar_AC', 'Kaviar_AN', 'HRC_AF', 'HRC_AC', 'HRC_AN', 'HRC_non1000G_AF', 'HRC_non1000G_AC', 'HRC_non1000G_AN', 'cosmic70', 'dbscSNV_ADA_SCORE', 'dbscSNV_RF_SCORE', 'Interpro_domain', 'avsnp147', 'exac_info']
	head_out = delim.join(head + ['\n'])
	outfile = project_prefix + '.annotated.txt'
	with open(avinput, "r") as av, open(outfile, "w") as final:
		final.write(head_out)
		for line in av:
			line = line.strip('\n').split(delim)
			stuff = line[0:8] + [line[14]] + [line[9]] + line[11:13]
			info = line[15].split(';')
			info_list = split_info_field(info)
			other_stuff = line[16:]
			exac_info = line[15].split('ANNOVAR_DATE')[0]
			line_out = delim.join(stuff + other_stuff + info_list +[exac_info, '\n'])
			final.write(line_out)

##calls all annovar methods
def annotate_vcf_file(vcf, project_prefix):
	table_annovar_vcf(vcf, project_prefix)
	post_annovar_vcf = project_prefix + '.' + av_genome + '_multianno.vcf'
	convert_to_annovar(post_annovar_vcf, project_prefix)
	if project_prefix.split('.')[0] == 'exac':
		format_avinput_exac(project_prefix)
	else:
		format_avinput(project_prefix)

def combine_files_and_add_count(infiles, outfile, pos_to_insert):
	var_dict = {}

	for infile in infiles:
		sample = infile.split('.')[2]
		file_count = 0
		with open(infile, "r") as in_fh:
			line_count = 0
			file_count += 1
			for line in in_fh:
				line = line.rstrip().split(delim)
				line_count += 1
				if line_count == 1 and file_count ==1:
					header = delim.join(line + ['count', 'samples', '\n'])
				else:
					if line_count > 1:
						var = '.'.join(line[:5])
						# print var
						if var in var_dict:
							var_dict[var][0].append(sample)
						else:
							var_dict[var] = [[sample], line]
	with open(outfile, "w") as out_fh:
		out_fh.write(header)
		for variant in var_dict:
			# var_dict[variant][0] = str(var_dict[variant][0])
			inline = var_dict[variant][1]
			count_samples = [str(len(var_dict[variant][0])), ', '.join(var_dict[variant][0])]
			out_line = inline[:pos_to_insert] + count_samples + inline[pos_to_insert:] + ['\n']
			out_fh.write(delim.join(out_line))

def filter_ann_txt_files(infile, cols_to_filter, freq_req, outfile_prefix):
	##filter rnaseq data
	for freq in freq_req:
		##only 'rare' (<=1%)
		filtering_annotated.filter(working_dir, "and", infile, "1.temp", [cols_to_filter[1]], ['<='], [freq])
		##q>=30 and coverage >=5
		filtering_annotated.filter(working_dir, "and", "1.temp", "2.temp", cols_to_filter[5], ['>=', '>='], [30,5])
		##exonic_variants in refGene
		filtering_annotated.filter(working_dir, "or", "2.temp", outfile_prefix + '.' + str(freq) +  ".exonic.xls", [cols_to_filter[0], cols_to_filter[0]], ['==','=='], [exonic_definitions[0],exonic_definitions[1]])
		##get dispuptive
		filtering_annotated.filter(working_dir, "or", outfile_prefix + '.' + str(freq) +  ".exonic.xls", outfile_prefix + '.' + str(freq) +  ".disruptive.xls", make_list(cols_to_filter[2], disruptive_definitions), make_list('==', disruptive_definitions), disruptive_definitions)
		##get all protein changing
		filtering_annotated.filter(working_dir, "or", outfile_prefix + '.' + str(freq) +  ".exonic.xls", outfile_prefix + '.' + str(freq) +  ".protein_changing.xls", make_list(cols_to_filter[2], protein_changing_definitions), make_list('==', protein_changing_definitions), protein_changing_definitions)
		##get damaging - pp2_hdiv, pp2_hvar, cadd_phred, gerp - in all or in any
		#get non synonymous snps
		filtering_annotated.filter(working_dir, "or", outfile_prefix + '.' + str(freq) +  ".exonic.xls", "3.temp", [cols_to_filter[2], cols_to_filter[2]], ['==', '=='], nosyn_definitions)
		#get if all are positive
		# filtering_annotated.filter(working_dir, "and", group + '.' + str(freq) +  ".3.temp", group + '.' + str(freq) +  ".4.temp", case_controls_col[3], make_list('>=', case_controls_col[3]), damaging_definitions)
		filtering_annotated.filter(working_dir, "or", "3.temp",  "4a.temp", cols_to_filter[3][:2], make_list('>=', cols_to_filter[3][:2]), pp2_score)
		filtering_annotated.filter(working_dir, "and", "4a.temp", outfile_prefix + '.' + str(freq) + ".damaging_all.xls", cols_to_filter[3][2:], make_list('>=', cols_to_filter[3][2:]), cadd_gerp_score)
		# remove_rows_with_no_data(group + '.' + str(freq) +  ".4b.temp", group + '.' + str(freq) +  ".damaging_all.xls", case_controls_col[3])
		#get if any are positive
		filtering_annotated.filter(working_dir, "or", "3.temp", outfile_prefix + '.' + str(freq) +  ".damaging_any.xls", cols_to_filter[3], make_list('>=', cols_to_filter[3]), damaging_definitions)
		# remove_rows_with_no_data(group + '.' + str(freq) +  ".5.temp", group + '.' + str(freq) +  ".damaging_any.xls", case_controls_col[3])

def convert_list_to_dict_with_zero(gene_list):
	##convert list to dict to track var_count
	gene_dict = {}
	for gene in gene_list:
		gene_dict[gene] = 0
	return gene_dict

def filter_var_by_genename(genelist, input_file,output_file, gene_col):
	##convert list to dict to track var_count
	genedict = convert_list_to_dict_with_zero(genelist)
# 	print genedict
# 	print len(genedict), len(genelist)
	with open(output_file, "w") as outfh, open(input_file, "U") as infh:
		line_count = 0
		for line in infh:
			line_count += 1
			if line_count == 1:
				outfh.write(line)
			else:
				line = line.strip('\n').strip('\r').rstrip().split(delim)
				genes = re.sub(r'\(.+?\)', '', line[gene_col])
				genes = genes.split(',')
				for gene in genes:
					if gene in genedict:
						outfh.write(delim.join(line) + '\n')
						genedict[gene] += 1
						break
		#print genedict
		print 'of the variants checked %i are in the genes specified'% (sum(genedict.values()))
		print '\n'
# 		for g in genedict:
# 			print 'in file %s we have %s variants in gene %s'%(output_file, genedict[g], g)

def combine_ann_txt(sample_list, prefix, suffix, outfile):
	with open(outfile, "w") as final_file:
		sample_count = 0
		for sample in sample_list:
			file = prefix + sample + suffix
			sample_count += 1
			with open(file, "r") as open_file:
				line_count = 0
				if sample_count == 1:
					for line in open_file:
						line_count += 1
						if line_count == 1:
							final_file.write('Proband' + delim + line)
						else:
							final_file.write(sample + delim + line)
				else:
					for line in open_file:
						line_count += 1
						if line_count > 1:
							final_file.write(sample + delim + line)


##run methods

##filter for passed by gatk for cohort 1 variants
# filter_vcf_file(c1_case_vcf, c1_case_passed_vcf)
# filter_vcf_file(c1_controls_vcf, c1_controls_passed_vcf)
# filter_vcf_file(exac_vcf, exac_passed_vcf)

##annotate variants
# annotate_vcf_file(c1_case_passed_vcf, 'cohort1.cases')
# annotate_vcf_file(c1_controls_passed_vcf, 'cohort1.ctls')
# annotate_vcf_file(c2_case_passed_vcf, 'cohort2.cases')
# annotate_vcf_file(c2_controls_passed_vcf, 'cohort2.ctls')
# annotate_vcf_file(exac_vcf, 'exac.all')
# annotate_vcf_file(exac_passed_vcf, 'exac.passed')

##combine annoated.txt files and add count and sample names
# all_ann_cases = glob.glob('*cases*annotated.txt')
# all_ann_controls = glob.glob('*ctls*annotated.txt')
# print len(all_ann_cases)
# print len(all_ann_controls)
# combine_files_and_add_count(all_ann_cases, comb_annotated_cases, 88)
# combine_files_and_add_count(all_ann_controls, comb_annotated_ctls, 88)

##filter variants for gene and variant level analysis
# filter_ann_txt_files(comb_annotated_cases, filtering_cols, frequencies, 'cs_rnaseq_1116.c1_and_2.cases')
# filter_ann_txt_files(comb_annotated_ctls, filtering_cols, frequencies, 'cs_rnaseq_1116.c1_and_2.ctls')
# filter_ann_txt_files(exac_annotated, filtering_cols, frequencies, 'exac.all')
# filter_ann_txt_files(exac_annotated_passed, filtering_cols, frequencies, 'exac.passed')



##all vars in a gene
# cases_list = ['cohort1.cases.95445', 'cohort1.cases.95446', 'cohort1.cases.95447', 'cohort1.cases.95448', 'cohort1.cases.95450', 'cohort1.cases.95452', 'cohort1.cases.95453', 'cohort1.cases.95454', 'cohort1.cases.95456', 'cohort1.cases.95457', 'cohort1.cases.95458', 'cohort1.cases.95459', 'cohort1.cases.95460', 'cohort1.cases.95462', 'cohort1.cases.95463', 'cohort1.cases.95464', 'cohort1.cases.95466', 'cohort1.cases.95467', 'cohort1.cases.95468', 'cohort1.cases.95469', 'cohort1.cases.95470', 'cohort1.cases.95471', 'cohort1.cases.95472', 'cohort1.cases.95474', 'cohort1.cases.95475', 'cohort1.cases.95476', 'cohort1.cases.95477', 'cohort1.cases.95479', 'cohort1.cases.95481', 'cohort1.cases.95482', 'cohort1.cases.95484', 'cohort1.cases.95485', 'cohort1.cases.95486', 'cohort1.cases.95487', 'cohort1.cases.95488', 'cohort1.cases.95490', 'cohort1.cases.95491', 'cohort1.cases.95492', 'cohort1.cases.95493', 'cohort1.cases.95494', 'cohort1.cases.95495', 'cohort1.cases.95496', 'cohort1.cases.95498', 'cohort1.cases.95500', 'cohort1.cases.95501', 'cohort1.cases.95502', 'cohort1.cases.95505', 'cohort1.cases.95506', 'cohort1.cases.95507', 'cohort1.cases.95509', 'cohort1.cases.95510', 'cohort1.cases.95511', 'cohort1.cases.95512', 'cohort1.cases.95514', 'cohort1.cases.95515', 'cohort1.cases.95516', 'cohort1.cases.95518', 'cohort1.cases.95519', 'cohort1.cases.95520', 'cohort1.cases.95521', 'cohort1.cases.95522', 'cohort1.cases.95524', 'cohort1.cases.95525', 'cohort1.cases.95526', 'cohort1.cases.95527', 'cohort1.cases.95528', 'cohort1.cases.95529', 'cohort1.cases.95530', 'cohort1.cases.95531', 'cohort1.cases.95534', 'cohort1.cases.95535', 'cohort1.cases.95536', 'cohort1.cases.95537', 'cohort1.cases.95540', 'cohort1.cases.95542', 'cohort1.cases.95543', 'cohort1.cases.95544', 'cohort1.cases.95545', 'cohort1.cases.95546', 'cohort1.cases.95547', 'cohort1.cases.95549', 'cohort1.cases.95550', 'cohort1.cases.95551', 'cohort1.cases.95552', 'cohort1.cases.95555', 'cohort1.cases.95556', 'cohort1.cases.95557', 'cohort1.cases.95558', 'cohort1.cases.95559', 'cohort1.cases.95560', 'cohort1.cases.95562', 'cohort1.cases.95563', 'cohort1.cases.95564', 'cohort1.cases.95565', 'cohort1.cases.95566', 'cohort1.cases.95567', 'cohort1.cases.95569', 'cohort1.cases.95570', 'cohort1.cases.95571', 'cohort1.cases.95572', 'cohort1.cases.95573', 'cohort1.cases.95575', 'cohort1.cases.95576', 'cohort1.cases.95578', 'cohort1.cases.95579', 'cohort1.cases.95580', 'cohort1.cases.95581', 'cohort1.cases.95583', 'cohort1.cases.95584', 'cohort1.cases.95585', 'cohort1.cases.95586', 'cohort1.cases.95587', 'cohort1.cases.95588', 'cohort1.cases.95589', 'cohort1.cases.95591', 'cohort1.cases.95592', 'cohort1.cases.95593', 'cohort1.cases.95594', 'cohort1.cases.95595', 'cohort1.cases.95596', 'cohort1.cases.95597', 'cohort1.cases.95598', 'cohort1.cases.95599', 'cohort1.cases.95600', 'cohort1.cases.95601', 'cohort1.cases.95602', 'cohort1.cases.95603', 'cohort1.cases.95604', 'cohort1.cases.95605', 'cohort1.cases.95606', 'cohort1.cases.95609', 'cohort1.cases.95610', 'cohort1.cases.95611', 'cohort1.cases.95612', 'cohort1.cases.95614', 'cohort1.cases.95615', 'cohort1.cases.95616', 'cohort1.cases.95617', 'cohort1.cases.95619', 'cohort1.cases.95620', 'cohort1.cases.95622', 'cohort1.cases.95623', 'cohort1.cases.95624', 'cohort1.cases.95626', 'cohort1.cases.95627', 'cohort1.cases.95628', 'cohort1.cases.95629', 'cohort1.cases.95630', 'cohort1.cases.95631', 'cohort1.cases.95633', 'cohort1.cases.95634', 'cohort1.cases.95635', 'cohort1.cases.95636', 'cohort1.cases.95637', 'cohort1.cases.95638', 'cohort1.cases.95639', 'cohort1.cases.95640', 'cohort1.cases.95642', 'cohort1.cases.95644', 'cohort1.cases.95645', 'cohort1.cases.95646', 'cohort1.cases.95649', 'cohort1.cases.95650', 'cohort1.cases.95651', 'cohort1.cases.95652', 'cohort1.cases.95653', 'cohort1.cases.95654', 'cohort1.cases.95655', 'cohort1.cases.95656', 'cohort1.cases.95657', 'cohort1.cases.95658', 'cohort1.cases.95661', 'cohort1.cases.95662', 'cohort1.cases.95663', 'cohort1.cases.95665', 'cohort1.cases.95668', 'cohort1.cases.95669', 'cohort1.cases.95670', 'cohort1.cases.95671', 'cohort1.cases.95673', 'cohort1.cases.95675', 'cohort1.cases.95676', 'cohort1.cases.95677', 'cohort1.cases.95678', 'cohort1.cases.95679', 'cohort1.cases.95680', 'cohort1.cases.95681', 'cohort1.cases.95682', 'cohort1.cases.95683', 'cohort1.cases.95684', 'cohort1.cases.95685', 'cohort1.cases.95686', 'cohort1.cases.95688', 'cohort1.cases.95689', 'cohort1.cases.95690', 'cohort1.cases.95691', 'cohort1.cases.95693', 'cohort1.cases.95694', 'cohort1.cases.95696', 'cohort1.cases.95697', 'cohort1.cases.95698', 'cohort1.cases.95699', 'cohort1.cases.95700', 'cohort1.cases.95702', 'cohort1.cases.95703', 'cohort2.cases.163499', 'cohort2.cases.163502', 'cohort2.cases.163503', 'cohort2.cases.163504', 'cohort2.cases.163505', 'cohort2.cases.163506', 'cohort2.cases.163507', 'cohort2.cases.163508', 'cohort2.cases.163509', 'cohort2.cases.163510', 'cohort2.cases.163511', 'cohort2.cases.163513', 'cohort2.cases.163516', 'cohort2.cases.163517', 'cohort2.cases.163518', 'cohort2.cases.163519', 'cohort2.cases.163520', 'cohort2.cases.163521', 'cohort2.cases.163522', 'cohort2.cases.163523', 'cohort2.cases.163524', 'cohort2.cases.163525', 'cohort2.cases.163526', 'cohort2.cases.163527', 'cohort2.cases.163528', 'cohort2.cases.163529', 'cohort2.cases.163530', 'cohort2.cases.163531', 'cohort2.cases.163532', 'cohort2.cases.163533', 'cohort2.cases.163534', 'cohort2.cases.163535', 'cohort2.cases.163537', 'cohort2.cases.163538', 'cohort2.cases.163539', 'cohort2.cases.163540', 'cohort2.cases.163541', 'cohort2.cases.163543', 'cohort2.cases.163544', 'cohort2.cases.163545', 'cohort2.cases.163546', 'cohort2.cases.163547', 'cohort2.cases.163548', 'cohort2.cases.163549', 'cohort2.cases.163553', 'cohort2.cases.163555', 'cohort2.cases.163556', 'cohort2.cases.163557', 'cohort2.cases.163558', 'cohort2.cases.163559', 'cohort2.cases.163560', 'cohort2.cases.163561', 'cohort2.cases.163562', 'cohort2.cases.163563', 'cohort2.cases.163564', 'cohort2.cases.163565', 'cohort2.cases.163566', 'cohort2.cases.163567', 'cohort2.cases.163568', 'cohort2.cases.163569', 'cohort2.cases.163571', 'cohort2.cases.163572', 'cohort2.cases.163573', 'cohort2.cases.163575', 'cohort2.cases.163576', 'cohort2.cases.163577', 'cohort2.cases.163578', 'cohort2.cases.163579', 'cohort2.cases.163580', 'cohort2.cases.163581', 'cohort2.cases.163582', 'cohort2.cases.163584', 'cohort2.cases.163585', 'cohort2.cases.163586', 'cohort2.cases.163587', 'cohort2.cases.163588', 'cohort2.cases.163590', 'cohort2.cases.163591', 'cohort2.cases.163592', 'cohort2.cases.163593', 'cohort2.cases.163594', 'cohort2.cases.163595', 'cohort2.cases.163596', 'cohort2.cases.163597', 'cohort2.cases.163598', 'cohort2.cases.163601', 'cohort2.cases.163602', 'cohort2.cases.163603', 'cohort2.cases.163604', 'cohort2.cases.163605', 'cohort2.cases.163606', 'cohort2.cases.163607', 'cohort2.cases.163609', 'cohort2.cases.163610', 'cohort2.cases.163611', 'cohort2.cases.163614', 'cohort2.cases.163615', 'cohort2.cases.163616', 'cohort2.cases.163619', 'cohort2.cases.163620', 'cohort2.cases.163621', 'cohort2.cases.163622', 'cohort2.cases.163623', 'cohort2.cases.163624', 'cohort2.cases.163625', 'cohort2.cases.163626', 'cohort2.cases.163627', 'cohort2.cases.163628', 'cohort2.cases.163629', 'cohort2.cases.163630', 'cohort2.cases.163631', 'cohort2.cases.163632', 'cohort2.cases.163633', 'cohort2.cases.163634', 'cohort2.cases.163635', 'cohort2.cases.163636', 'cohort2.cases.163637', 'cohort2.cases.163639', 'cohort2.cases.163640', 'cohort2.cases.163641', 'cohort2.cases.163642', 'cohort2.cases.163643', 'cohort2.cases.163644', 'cohort2.cases.163645', 'cohort2.cases.163646', 'cohort2.cases.163647', 'cohort2.cases.163649', 'cohort2.cases.163651', 'cohort2.cases.163652', 'cohort2.cases.163653', 'cohort2.cases.163654', 'cohort2.cases.163655', 'cohort2.cases.163656', 'cohort2.cases.163657', 'cohort2.cases.163658', 'cohort2.cases.163659', 'cohort2.cases.163660', 'cohort2.cases.163661', 'cohort2.cases.163662', 'cohort2.cases.163663', 'cohort2.cases.163664', 'cohort2.cases.163665', 'cohort2.cases.163666', 'cohort2.cases.163667', 'cohort2.cases.163668', 'cohort2.cases.163669', 'cohort2.cases.163672', 'cohort2.cases.163674', 'cohort2.cases.163678', 'cohort2.cases.163679', 'cohort2.cases.163680', 'cohort2.cases.163681', 'cohort2.cases.163682', 'cohort2.cases.163683', 'cohort2.cases.163684', 'cohort2.cases.163685', 'cohort2.cases.163686', 'cohort2.cases.163687', 'cohort2.cases.163688', 'cohort2.cases.163690', 'cohort2.cases.163691', 'cohort2.cases.163692', 'cohort2.cases.163694', 'cohort2.cases.163695', 'cohort2.cases.163696', 'cohort2.cases.163697', 'cohort2.cases.163698', 'cohort2.cases.163699', 'cohort2.cases.163700', 'cohort2.cases.163701', 'cohort2.cases.163702', 'cohort2.cases.163703', 'cohort2.cases.163704', 'cohort2.cases.163706', 'cohort2.cases.163707', 'cohort2.cases.163709', 'cohort2.cases.163710', 'cohort2.cases.163711', 'cohort2.cases.163712', 'cohort2.cases.163713', 'cohort2.cases.163714', 'cohort2.cases.163717', 'cohort2.cases.163718', 'cohort2.cases.163719', 'cohort2.cases.163720', 'cohort2.cases.163721']
# control_list = ['cohort1.ctls.95449', 'cohort1.ctls.95451', 'cohort1.ctls.95455', 'cohort1.ctls.95461', 'cohort1.ctls.95465', 'cohort1.ctls.95473', 'cohort1.ctls.95478', 'cohort1.ctls.95480', 'cohort1.ctls.95483', 'cohort1.ctls.95489', 'cohort1.ctls.95497', 'cohort1.ctls.95499', 'cohort1.ctls.95504', 'cohort1.ctls.95508', 'cohort1.ctls.95513', 'cohort1.ctls.95517', 'cohort1.ctls.95523', 'cohort1.ctls.95532', 'cohort1.ctls.95533', 'cohort1.ctls.95538', 'cohort1.ctls.95539', 'cohort1.ctls.95541', 'cohort1.ctls.95548', 'cohort1.ctls.95553', 'cohort1.ctls.95554', 'cohort1.ctls.95561', 'cohort1.ctls.95568', 'cohort1.ctls.95574', 'cohort1.ctls.95577', 'cohort1.ctls.95582', 'cohort1.ctls.95590', 'cohort1.ctls.95607', 'cohort1.ctls.95608', 'cohort1.ctls.95613', 'cohort1.ctls.95618', 'cohort1.ctls.95621', 'cohort1.ctls.95625', 'cohort1.ctls.95647', 'cohort1.ctls.95648', 'cohort1.ctls.95659', 'cohort1.ctls.95660', 'cohort1.ctls.95666', 'cohort1.ctls.95667', 'cohort1.ctls.95672', 'cohort1.ctls.95674', 'cohort1.ctls.95692', 'cohort1.ctls.95695', 'cohort1.ctls.95704', 'cohort2.ctls.163500', 'cohort2.ctls.163501', 'cohort2.ctls.163512', 'cohort2.ctls.163514', 'cohort2.ctls.163515', 'cohort2.ctls.163536', 'cohort2.ctls.163542', 'cohort2.ctls.163550', 'cohort2.ctls.163551', 'cohort2.ctls.163552', 'cohort2.ctls.163554', 'cohort2.ctls.163570', 'cohort2.ctls.163574', 'cohort2.ctls.163583', 'cohort2.ctls.163589', 'cohort2.ctls.163599', 'cohort2.ctls.163600', 'cohort2.ctls.163608', 'cohort2.ctls.163612', 'cohort2.ctls.163613', 'cohort2.ctls.163617', 'cohort2.ctls.163618', 'cohort2.ctls.163638', 'cohort2.ctls.163648', 'cohort2.ctls.163650', 'cohort2.ctls.163670', 'cohort2.ctls.163671', 'cohort2.ctls.163673', 'cohort2.ctls.163675', 'cohort2.ctls.163676', 'cohort2.ctls.163677', 'cohort2.ctls.163689', 'cohort2.ctls.163693', 'cohort2.ctls.163705', 'cohort2.ctls.163708', 'cohort2.ctls.163715', 'cohort2.ctls.163716']
##BMP2_SMAD6
# gene_list = ['BMP2', 'SMAD6']
# final_suffix = '.all_BMP2_SMAD6_0717.xls'
##BMP7
# gene_list = ['BMP7']
# final_suffix = '.all_BMP7_0717.xls'

##check 0917
# gene_list = ['FGFR1', 'FGFR2', 'FGFR3', 'TWIST1', 'TCF12', 'EFNB1']
# final_suffix = '.0917.xls'
# cases_list = ['cohort1.cases.95521', 'cohort1.cases.95493', 'cohort1.cases.95682', 'cohort1.cases.95599', 'cohort2.cases.163658', 'cohort2.cases.163582', 'cohort2.cases.163502', 'cohort2.cases.163598', 'cohort2.cases.163628', 'cohort2.cases.163557', 'cohort2.cases.163594', 'cohort2.cases.163700', 'cohort2.cases.163559', 'cohort2.cases.163562', 'cohort2.cases.163505', 'cohort2.cases.163525', 'cohort2.cases.163616', 'cohort2.cases.163603', 'cohort2.cases.163563', 'cohort2.cases.163561', 'cohort2.cases.163509', 'cohort2.cases.163508', 'cohort2.cases.163548', 'cohort2.cases.163662', 'cohort2.cases.163663', 'cohort2.cases.163655', 'cohort2.cases.163646']

##get variants in those genes
# for case in cases_list:
# 	filter_var_by_genename(gene_list, case + '.annotated.txt', case + '.temp.xls', 15)
# combine_ann_txt(cases_list, '', '.temp.xls', 'cases' + final_suffix)
# for control in control_list:
# 	filter_var_by_genename(gene_list, control + '.annotated.txt', control + '.temp.xls', 15)
# combine_ann_txt(control_list, '', '.temp.xls', 'controls' + final_suffix)

##deema 1017 -- get all variants 
final_suffix = '.FLNA_associated_1017.xls'
cases_list = ['cohort1.cases.95445', 'cohort1.cases.95446', 'cohort1.cases.95447', 'cohort1.cases.95448', 'cohort1.cases.95450', 'cohort1.cases.95452', 'cohort1.cases.95453', 'cohort1.cases.95454', 'cohort1.cases.95456', 'cohort1.cases.95457', 'cohort1.cases.95458', 'cohort1.cases.95459', 'cohort1.cases.95460', 'cohort1.cases.95462', 'cohort1.cases.95463', 'cohort1.cases.95464', 'cohort1.cases.95466', 'cohort1.cases.95467', 'cohort1.cases.95468', 'cohort1.cases.95469', 'cohort1.cases.95470', 'cohort1.cases.95471', 'cohort1.cases.95472', 'cohort1.cases.95474', 'cohort1.cases.95475', 'cohort1.cases.95476', 'cohort1.cases.95477', 'cohort1.cases.95479', 'cohort1.cases.95481', 'cohort1.cases.95482', 'cohort1.cases.95484', 'cohort1.cases.95485', 'cohort1.cases.95486', 'cohort1.cases.95487', 'cohort1.cases.95488', 'cohort1.cases.95490', 'cohort1.cases.95491', 'cohort1.cases.95492', 'cohort1.cases.95493', 'cohort1.cases.95494', 'cohort1.cases.95495', 'cohort1.cases.95496', 'cohort1.cases.95498', 'cohort1.cases.95500', 'cohort1.cases.95501', 'cohort1.cases.95502', 'cohort1.cases.95505', 'cohort1.cases.95506', 'cohort1.cases.95507', 'cohort1.cases.95509', 'cohort1.cases.95510', 'cohort1.cases.95511', 'cohort1.cases.95512', 'cohort1.cases.95514', 'cohort1.cases.95515', 'cohort1.cases.95516', 'cohort1.cases.95518', 'cohort1.cases.95519', 'cohort1.cases.95520', 'cohort1.cases.95521', 'cohort1.cases.95522', 'cohort1.cases.95524', 'cohort1.cases.95525', 'cohort1.cases.95526', 'cohort1.cases.95527', 'cohort1.cases.95528', 'cohort1.cases.95529', 'cohort1.cases.95530', 'cohort1.cases.95531', 'cohort1.cases.95534', 'cohort1.cases.95535', 'cohort1.cases.95536', 'cohort1.cases.95537', 'cohort1.cases.95540', 'cohort1.cases.95542', 'cohort1.cases.95543', 'cohort1.cases.95544', 'cohort1.cases.95545', 'cohort1.cases.95546', 'cohort1.cases.95547', 'cohort1.cases.95549', 'cohort1.cases.95550', 'cohort1.cases.95551', 'cohort1.cases.95552', 'cohort1.cases.95555', 'cohort1.cases.95556', 'cohort1.cases.95557', 'cohort1.cases.95558', 'cohort1.cases.95559', 'cohort1.cases.95560', 'cohort1.cases.95562', 'cohort1.cases.95563', 'cohort1.cases.95564', 'cohort1.cases.95565', 'cohort1.cases.95566', 'cohort1.cases.95567', 'cohort1.cases.95569', 'cohort1.cases.95570', 'cohort1.cases.95571', 'cohort1.cases.95572', 'cohort1.cases.95573', 'cohort1.cases.95575', 'cohort1.cases.95576', 'cohort1.cases.95578', 'cohort1.cases.95579', 'cohort1.cases.95580', 'cohort1.cases.95581', 'cohort1.cases.95583', 'cohort1.cases.95584', 'cohort1.cases.95585', 'cohort1.cases.95586', 'cohort1.cases.95587', 'cohort1.cases.95588', 'cohort1.cases.95589', 'cohort1.cases.95591', 'cohort1.cases.95592', 'cohort1.cases.95593', 'cohort1.cases.95594', 'cohort1.cases.95595', 'cohort1.cases.95596', 'cohort1.cases.95597', 'cohort1.cases.95598', 'cohort1.cases.95599', 'cohort1.cases.95600', 'cohort1.cases.95601', 'cohort1.cases.95602', 'cohort1.cases.95603', 'cohort1.cases.95604', 'cohort1.cases.95605', 'cohort1.cases.95606', 'cohort1.cases.95609', 'cohort1.cases.95610', 'cohort1.cases.95611', 'cohort1.cases.95612', 'cohort1.cases.95614', 'cohort1.cases.95615', 'cohort1.cases.95616', 'cohort1.cases.95617', 'cohort1.cases.95619', 'cohort1.cases.95620', 'cohort1.cases.95622', 'cohort1.cases.95623', 'cohort1.cases.95624', 'cohort1.cases.95626', 'cohort1.cases.95627', 'cohort1.cases.95628', 'cohort1.cases.95629', 'cohort1.cases.95630', 'cohort1.cases.95631', 'cohort1.cases.95633', 'cohort1.cases.95634', 'cohort1.cases.95635', 'cohort1.cases.95636', 'cohort1.cases.95637', 'cohort1.cases.95638', 'cohort1.cases.95639', 'cohort1.cases.95640', 'cohort1.cases.95642', 'cohort1.cases.95644', 'cohort1.cases.95645', 'cohort1.cases.95646', 'cohort1.cases.95649', 'cohort1.cases.95650', 'cohort1.cases.95651', 'cohort1.cases.95652', 'cohort1.cases.95653', 'cohort1.cases.95654', 'cohort1.cases.95655', 'cohort1.cases.95656', 'cohort1.cases.95657', 'cohort1.cases.95658', 'cohort1.cases.95661', 'cohort1.cases.95662', 'cohort1.cases.95663', 'cohort1.cases.95665', 'cohort1.cases.95668', 'cohort1.cases.95669', 'cohort1.cases.95670', 'cohort1.cases.95671', 'cohort1.cases.95673', 'cohort1.cases.95675', 'cohort1.cases.95676', 'cohort1.cases.95677', 'cohort1.cases.95678', 'cohort1.cases.95679', 'cohort1.cases.95680', 'cohort1.cases.95681', 'cohort1.cases.95682', 'cohort1.cases.95683', 'cohort1.cases.95684', 'cohort1.cases.95685', 'cohort1.cases.95686', 'cohort1.cases.95688', 'cohort1.cases.95689', 'cohort1.cases.95690', 'cohort1.cases.95691', 'cohort1.cases.95693', 'cohort1.cases.95694', 'cohort1.cases.95696', 'cohort1.cases.95697', 'cohort1.cases.95698', 'cohort1.cases.95699', 'cohort1.cases.95700', 'cohort1.cases.95702', 'cohort1.cases.95703', 'cohort2.cases.163499', 'cohort2.cases.163502', 'cohort2.cases.163503', 'cohort2.cases.163504', 'cohort2.cases.163505', 'cohort2.cases.163506', 'cohort2.cases.163507', 'cohort2.cases.163508', 'cohort2.cases.163509', 'cohort2.cases.163510', 'cohort2.cases.163511', 'cohort2.cases.163513', 'cohort2.cases.163516', 'cohort2.cases.163517', 'cohort2.cases.163518', 'cohort2.cases.163519', 'cohort2.cases.163520', 'cohort2.cases.163521', 'cohort2.cases.163522', 'cohort2.cases.163523', 'cohort2.cases.163524', 'cohort2.cases.163525', 'cohort2.cases.163526', 'cohort2.cases.163527', 'cohort2.cases.163528', 'cohort2.cases.163529', 'cohort2.cases.163530', 'cohort2.cases.163531', 'cohort2.cases.163532', 'cohort2.cases.163533', 'cohort2.cases.163534', 'cohort2.cases.163535', 'cohort2.cases.163537', 'cohort2.cases.163538', 'cohort2.cases.163539', 'cohort2.cases.163540', 'cohort2.cases.163541', 'cohort2.cases.163543', 'cohort2.cases.163544', 'cohort2.cases.163545', 'cohort2.cases.163546', 'cohort2.cases.163547', 'cohort2.cases.163548', 'cohort2.cases.163549', 'cohort2.cases.163553', 'cohort2.cases.163555', 'cohort2.cases.163556', 'cohort2.cases.163557', 'cohort2.cases.163558', 'cohort2.cases.163559', 'cohort2.cases.163560', 'cohort2.cases.163561', 'cohort2.cases.163562', 'cohort2.cases.163563', 'cohort2.cases.163564', 'cohort2.cases.163565', 'cohort2.cases.163566', 'cohort2.cases.163567', 'cohort2.cases.163568', 'cohort2.cases.163569', 'cohort2.cases.163571', 'cohort2.cases.163572', 'cohort2.cases.163573', 'cohort2.cases.163575', 'cohort2.cases.163576', 'cohort2.cases.163577', 'cohort2.cases.163578', 'cohort2.cases.163579', 'cohort2.cases.163580', 'cohort2.cases.163581', 'cohort2.cases.163582', 'cohort2.cases.163584', 'cohort2.cases.163585', 'cohort2.cases.163586', 'cohort2.cases.163587', 'cohort2.cases.163588', 'cohort2.cases.163590', 'cohort2.cases.163591', 'cohort2.cases.163592', 'cohort2.cases.163593', 'cohort2.cases.163594', 'cohort2.cases.163595', 'cohort2.cases.163596', 'cohort2.cases.163597', 'cohort2.cases.163598', 'cohort2.cases.163601', 'cohort2.cases.163602', 'cohort2.cases.163603', 'cohort2.cases.163604', 'cohort2.cases.163605', 'cohort2.cases.163606', 'cohort2.cases.163607', 'cohort2.cases.163609', 'cohort2.cases.163610', 'cohort2.cases.163611', 'cohort2.cases.163614', 'cohort2.cases.163615', 'cohort2.cases.163616', 'cohort2.cases.163619', 'cohort2.cases.163620', 'cohort2.cases.163621', 'cohort2.cases.163622', 'cohort2.cases.163623', 'cohort2.cases.163624', 'cohort2.cases.163625', 'cohort2.cases.163626', 'cohort2.cases.163627', 'cohort2.cases.163628', 'cohort2.cases.163629', 'cohort2.cases.163630', 'cohort2.cases.163631', 'cohort2.cases.163632', 'cohort2.cases.163633', 'cohort2.cases.163634', 'cohort2.cases.163635', 'cohort2.cases.163636', 'cohort2.cases.163637', 'cohort2.cases.163639', 'cohort2.cases.163640', 'cohort2.cases.163641', 'cohort2.cases.163642', 'cohort2.cases.163643', 'cohort2.cases.163644', 'cohort2.cases.163645', 'cohort2.cases.163646', 'cohort2.cases.163647', 'cohort2.cases.163649', 'cohort2.cases.163651', 'cohort2.cases.163652', 'cohort2.cases.163653', 'cohort2.cases.163654', 'cohort2.cases.163655', 'cohort2.cases.163656', 'cohort2.cases.163657', 'cohort2.cases.163658', 'cohort2.cases.163659', 'cohort2.cases.163660', 'cohort2.cases.163661', 'cohort2.cases.163662', 'cohort2.cases.163663', 'cohort2.cases.163664', 'cohort2.cases.163665', 'cohort2.cases.163666', 'cohort2.cases.163667', 'cohort2.cases.163668', 'cohort2.cases.163669', 'cohort2.cases.163672', 'cohort2.cases.163674', 'cohort2.cases.163678', 'cohort2.cases.163679', 'cohort2.cases.163680', 'cohort2.cases.163681', 'cohort2.cases.163682', 'cohort2.cases.163683', 'cohort2.cases.163684', 'cohort2.cases.163685', 'cohort2.cases.163686', 'cohort2.cases.163687', 'cohort2.cases.163688', 'cohort2.cases.163690', 'cohort2.cases.163691', 'cohort2.cases.163692', 'cohort2.cases.163694', 'cohort2.cases.163695', 'cohort2.cases.163696', 'cohort2.cases.163697', 'cohort2.cases.163698', 'cohort2.cases.163699', 'cohort2.cases.163700', 'cohort2.cases.163701', 'cohort2.cases.163702', 'cohort2.cases.163703', 'cohort2.cases.163704', 'cohort2.cases.163706', 'cohort2.cases.163707', 'cohort2.cases.163709', 'cohort2.cases.163710', 'cohort2.cases.163711', 'cohort2.cases.163712', 'cohort2.cases.163713', 'cohort2.cases.163714', 'cohort2.cases.163717', 'cohort2.cases.163718', 'cohort2.cases.163719', 'cohort2.cases.163720', 'cohort2.cases.163721']
cases_list = ['cohort2.cases.163641', 'cohort2.cases.163714', 'cohort2.cases.163560', 'cohort1.cases.95479', 'cohort1.cases.95597', 'cohort2.cases.163606']
# combine_ann_txt(cases_list, '', '.annotated.txt', 'cases' + final_suffix)

##jg get all vars in a set of genesgene
cases_list = ['cohort1.cases.95445', 'cohort1.cases.95446', 'cohort1.cases.95447', 'cohort1.cases.95448', 'cohort1.cases.95450', 'cohort1.cases.95452', 'cohort1.cases.95453', 'cohort1.cases.95454', 'cohort1.cases.95456', 'cohort1.cases.95457', 'cohort1.cases.95458', 'cohort1.cases.95459', 'cohort1.cases.95460', 'cohort1.cases.95462', 'cohort1.cases.95463', 'cohort1.cases.95464', 'cohort1.cases.95466', 'cohort1.cases.95467', 'cohort1.cases.95468', 'cohort1.cases.95469', 'cohort1.cases.95470', 'cohort1.cases.95471', 'cohort1.cases.95472', 'cohort1.cases.95474', 'cohort1.cases.95475', 'cohort1.cases.95476', 'cohort1.cases.95477', 'cohort1.cases.95479', 'cohort1.cases.95481', 'cohort1.cases.95482', 'cohort1.cases.95484', 'cohort1.cases.95485', 'cohort1.cases.95486', 'cohort1.cases.95487', 'cohort1.cases.95488', 'cohort1.cases.95490', 'cohort1.cases.95491', 'cohort1.cases.95492', 'cohort1.cases.95493', 'cohort1.cases.95494', 'cohort1.cases.95495', 'cohort1.cases.95496', 'cohort1.cases.95498', 'cohort1.cases.95500', 'cohort1.cases.95501', 'cohort1.cases.95502', 'cohort1.cases.95505', 'cohort1.cases.95506', 'cohort1.cases.95507', 'cohort1.cases.95509', 'cohort1.cases.95510', 'cohort1.cases.95511', 'cohort1.cases.95512', 'cohort1.cases.95514', 'cohort1.cases.95515', 'cohort1.cases.95516', 'cohort1.cases.95518', 'cohort1.cases.95519', 'cohort1.cases.95520', 'cohort1.cases.95521', 'cohort1.cases.95522', 'cohort1.cases.95524', 'cohort1.cases.95525', 'cohort1.cases.95526', 'cohort1.cases.95527', 'cohort1.cases.95528', 'cohort1.cases.95529', 'cohort1.cases.95530', 'cohort1.cases.95531', 'cohort1.cases.95534', 'cohort1.cases.95535', 'cohort1.cases.95536', 'cohort1.cases.95537', 'cohort1.cases.95540', 'cohort1.cases.95542', 'cohort1.cases.95543', 'cohort1.cases.95544', 'cohort1.cases.95545', 'cohort1.cases.95546', 'cohort1.cases.95547', 'cohort1.cases.95549', 'cohort1.cases.95550', 'cohort1.cases.95551', 'cohort1.cases.95552', 'cohort1.cases.95555', 'cohort1.cases.95556', 'cohort1.cases.95557', 'cohort1.cases.95558', 'cohort1.cases.95559', 'cohort1.cases.95560', 'cohort1.cases.95562', 'cohort1.cases.95563', 'cohort1.cases.95564', 'cohort1.cases.95565', 'cohort1.cases.95566', 'cohort1.cases.95567', 'cohort1.cases.95569', 'cohort1.cases.95570', 'cohort1.cases.95571', 'cohort1.cases.95572', 'cohort1.cases.95573', 'cohort1.cases.95575', 'cohort1.cases.95576', 'cohort1.cases.95578', 'cohort1.cases.95579', 'cohort1.cases.95580', 'cohort1.cases.95581', 'cohort1.cases.95583', 'cohort1.cases.95584', 'cohort1.cases.95585', 'cohort1.cases.95586', 'cohort1.cases.95587', 'cohort1.cases.95588', 'cohort1.cases.95589', 'cohort1.cases.95591', 'cohort1.cases.95592', 'cohort1.cases.95593', 'cohort1.cases.95594', 'cohort1.cases.95595', 'cohort1.cases.95596', 'cohort1.cases.95597', 'cohort1.cases.95598', 'cohort1.cases.95599', 'cohort1.cases.95600', 'cohort1.cases.95601', 'cohort1.cases.95602', 'cohort1.cases.95603', 'cohort1.cases.95604', 'cohort1.cases.95605', 'cohort1.cases.95606', 'cohort1.cases.95609', 'cohort1.cases.95610', 'cohort1.cases.95611', 'cohort1.cases.95612', 'cohort1.cases.95614', 'cohort1.cases.95615', 'cohort1.cases.95616', 'cohort1.cases.95617', 'cohort1.cases.95619', 'cohort1.cases.95620', 'cohort1.cases.95622', 'cohort1.cases.95623', 'cohort1.cases.95624', 'cohort1.cases.95626', 'cohort1.cases.95627', 'cohort1.cases.95628', 'cohort1.cases.95629', 'cohort1.cases.95630', 'cohort1.cases.95631', 'cohort1.cases.95633', 'cohort1.cases.95634', 'cohort1.cases.95635', 'cohort1.cases.95636', 'cohort1.cases.95637', 'cohort1.cases.95638', 'cohort1.cases.95639', 'cohort1.cases.95640', 'cohort1.cases.95642', 'cohort1.cases.95644', 'cohort1.cases.95645', 'cohort1.cases.95646', 'cohort1.cases.95649', 'cohort1.cases.95650', 'cohort1.cases.95651', 'cohort1.cases.95652', 'cohort1.cases.95653', 'cohort1.cases.95654', 'cohort1.cases.95655', 'cohort1.cases.95656', 'cohort1.cases.95657', 'cohort1.cases.95658', 'cohort1.cases.95661', 'cohort1.cases.95662', 'cohort1.cases.95663', 'cohort1.cases.95665', 'cohort1.cases.95668', 'cohort1.cases.95669', 'cohort1.cases.95670', 'cohort1.cases.95671', 'cohort1.cases.95673', 'cohort1.cases.95675', 'cohort1.cases.95676', 'cohort1.cases.95677', 'cohort1.cases.95678', 'cohort1.cases.95679', 'cohort1.cases.95680', 'cohort1.cases.95681', 'cohort1.cases.95682', 'cohort1.cases.95683', 'cohort1.cases.95684', 'cohort1.cases.95685', 'cohort1.cases.95686', 'cohort1.cases.95688', 'cohort1.cases.95689', 'cohort1.cases.95690', 'cohort1.cases.95691', 'cohort1.cases.95693', 'cohort1.cases.95694', 'cohort1.cases.95696', 'cohort1.cases.95697', 'cohort1.cases.95698', 'cohort1.cases.95699', 'cohort1.cases.95700', 'cohort1.cases.95702', 'cohort1.cases.95703', 'cohort2.cases.163499', 'cohort2.cases.163502', 'cohort2.cases.163503', 'cohort2.cases.163504', 'cohort2.cases.163505', 'cohort2.cases.163506', 'cohort2.cases.163507', 'cohort2.cases.163508', 'cohort2.cases.163509', 'cohort2.cases.163510', 'cohort2.cases.163511', 'cohort2.cases.163513', 'cohort2.cases.163516', 'cohort2.cases.163517', 'cohort2.cases.163518', 'cohort2.cases.163519', 'cohort2.cases.163520', 'cohort2.cases.163521', 'cohort2.cases.163522', 'cohort2.cases.163523', 'cohort2.cases.163524', 'cohort2.cases.163525', 'cohort2.cases.163526', 'cohort2.cases.163527', 'cohort2.cases.163528', 'cohort2.cases.163529', 'cohort2.cases.163530', 'cohort2.cases.163531', 'cohort2.cases.163532', 'cohort2.cases.163533', 'cohort2.cases.163534', 'cohort2.cases.163535', 'cohort2.cases.163537', 'cohort2.cases.163538', 'cohort2.cases.163539', 'cohort2.cases.163540', 'cohort2.cases.163541', 'cohort2.cases.163543', 'cohort2.cases.163544', 'cohort2.cases.163545', 'cohort2.cases.163546', 'cohort2.cases.163547', 'cohort2.cases.163548', 'cohort2.cases.163549', 'cohort2.cases.163553', 'cohort2.cases.163555', 'cohort2.cases.163556', 'cohort2.cases.163557', 'cohort2.cases.163558', 'cohort2.cases.163559', 'cohort2.cases.163560', 'cohort2.cases.163561', 'cohort2.cases.163562', 'cohort2.cases.163563', 'cohort2.cases.163564', 'cohort2.cases.163565', 'cohort2.cases.163566', 'cohort2.cases.163567', 'cohort2.cases.163568', 'cohort2.cases.163569', 'cohort2.cases.163571', 'cohort2.cases.163572', 'cohort2.cases.163573', 'cohort2.cases.163575', 'cohort2.cases.163576', 'cohort2.cases.163577', 'cohort2.cases.163578', 'cohort2.cases.163579', 'cohort2.cases.163580', 'cohort2.cases.163581', 'cohort2.cases.163582', 'cohort2.cases.163584', 'cohort2.cases.163585', 'cohort2.cases.163586', 'cohort2.cases.163587', 'cohort2.cases.163588', 'cohort2.cases.163590', 'cohort2.cases.163591', 'cohort2.cases.163592', 'cohort2.cases.163593', 'cohort2.cases.163594', 'cohort2.cases.163595', 'cohort2.cases.163596', 'cohort2.cases.163597', 'cohort2.cases.163598', 'cohort2.cases.163601', 'cohort2.cases.163602', 'cohort2.cases.163603', 'cohort2.cases.163604', 'cohort2.cases.163605', 'cohort2.cases.163606', 'cohort2.cases.163607', 'cohort2.cases.163609', 'cohort2.cases.163610', 'cohort2.cases.163611', 'cohort2.cases.163614', 'cohort2.cases.163615', 'cohort2.cases.163616', 'cohort2.cases.163619', 'cohort2.cases.163620', 'cohort2.cases.163621', 'cohort2.cases.163622', 'cohort2.cases.163623', 'cohort2.cases.163624', 'cohort2.cases.163625', 'cohort2.cases.163626', 'cohort2.cases.163627', 'cohort2.cases.163628', 'cohort2.cases.163629', 'cohort2.cases.163630', 'cohort2.cases.163631', 'cohort2.cases.163632', 'cohort2.cases.163633', 'cohort2.cases.163634', 'cohort2.cases.163635', 'cohort2.cases.163636', 'cohort2.cases.163637', 'cohort2.cases.163639', 'cohort2.cases.163640', 'cohort2.cases.163641', 'cohort2.cases.163642', 'cohort2.cases.163643', 'cohort2.cases.163644', 'cohort2.cases.163645', 'cohort2.cases.163646', 'cohort2.cases.163647', 'cohort2.cases.163649', 'cohort2.cases.163651', 'cohort2.cases.163652', 'cohort2.cases.163653', 'cohort2.cases.163654', 'cohort2.cases.163655', 'cohort2.cases.163656', 'cohort2.cases.163657', 'cohort2.cases.163658', 'cohort2.cases.163659', 'cohort2.cases.163660', 'cohort2.cases.163661', 'cohort2.cases.163662', 'cohort2.cases.163663', 'cohort2.cases.163664', 'cohort2.cases.163665', 'cohort2.cases.163666', 'cohort2.cases.163667', 'cohort2.cases.163668', 'cohort2.cases.163669', 'cohort2.cases.163672', 'cohort2.cases.163674', 'cohort2.cases.163678', 'cohort2.cases.163679', 'cohort2.cases.163680', 'cohort2.cases.163681', 'cohort2.cases.163682', 'cohort2.cases.163683', 'cohort2.cases.163684', 'cohort2.cases.163685', 'cohort2.cases.163686', 'cohort2.cases.163687', 'cohort2.cases.163688', 'cohort2.cases.163690', 'cohort2.cases.163691', 'cohort2.cases.163692', 'cohort2.cases.163694', 'cohort2.cases.163695', 'cohort2.cases.163696', 'cohort2.cases.163697', 'cohort2.cases.163698', 'cohort2.cases.163699', 'cohort2.cases.163700', 'cohort2.cases.163701', 'cohort2.cases.163702', 'cohort2.cases.163703', 'cohort2.cases.163704', 'cohort2.cases.163706', 'cohort2.cases.163707', 'cohort2.cases.163709', 'cohort2.cases.163710', 'cohort2.cases.163711', 'cohort2.cases.163712', 'cohort2.cases.163713', 'cohort2.cases.163714', 'cohort2.cases.163717', 'cohort2.cases.163718', 'cohort2.cases.163719', 'cohort2.cases.163720', 'cohort2.cases.163721']
control_list = ['cohort1.ctls.95449', 'cohort1.ctls.95451', 'cohort1.ctls.95455', 'cohort1.ctls.95461', 'cohort1.ctls.95465', 'cohort1.ctls.95473', 'cohort1.ctls.95478', 'cohort1.ctls.95480', 'cohort1.ctls.95483', 'cohort1.ctls.95489', 'cohort1.ctls.95497', 'cohort1.ctls.95499', 'cohort1.ctls.95504', 'cohort1.ctls.95508', 'cohort1.ctls.95513', 'cohort1.ctls.95517', 'cohort1.ctls.95523', 'cohort1.ctls.95532', 'cohort1.ctls.95533', 'cohort1.ctls.95538', 'cohort1.ctls.95539', 'cohort1.ctls.95541', 'cohort1.ctls.95548', 'cohort1.ctls.95553', 'cohort1.ctls.95554', 'cohort1.ctls.95561', 'cohort1.ctls.95568', 'cohort1.ctls.95574', 'cohort1.ctls.95577', 'cohort1.ctls.95582', 'cohort1.ctls.95590', 'cohort1.ctls.95607', 'cohort1.ctls.95608', 'cohort1.ctls.95613', 'cohort1.ctls.95618', 'cohort1.ctls.95621', 'cohort1.ctls.95625', 'cohort1.ctls.95647', 'cohort1.ctls.95648', 'cohort1.ctls.95659', 'cohort1.ctls.95660', 'cohort1.ctls.95666', 'cohort1.ctls.95667', 'cohort1.ctls.95672', 'cohort1.ctls.95674', 'cohort1.ctls.95692', 'cohort1.ctls.95695', 'cohort1.ctls.95704', 'cohort2.ctls.163500', 'cohort2.ctls.163501', 'cohort2.ctls.163512', 'cohort2.ctls.163514', 'cohort2.ctls.163515', 'cohort2.ctls.163536', 'cohort2.ctls.163542', 'cohort2.ctls.163550', 'cohort2.ctls.163551', 'cohort2.ctls.163552', 'cohort2.ctls.163554', 'cohort2.ctls.163570', 'cohort2.ctls.163574', 'cohort2.ctls.163583', 'cohort2.ctls.163589', 'cohort2.ctls.163599', 'cohort2.ctls.163600', 'cohort2.ctls.163608', 'cohort2.ctls.163612', 'cohort2.ctls.163613', 'cohort2.ctls.163617', 'cohort2.ctls.163618', 'cohort2.ctls.163638', 'cohort2.ctls.163648', 'cohort2.ctls.163650', 'cohort2.ctls.163670', 'cohort2.ctls.163671', 'cohort2.ctls.163673', 'cohort2.ctls.163675', 'cohort2.ctls.163676', 'cohort2.ctls.163677', 'cohort2.ctls.163689', 'cohort2.ctls.163693', 'cohort2.ctls.163705', 'cohort2.ctls.163708', 'cohort2.ctls.163715', 'cohort2.ctls.163716']
gene_list = ['IGF1', 'IGF1R', 'GSK3b', 'TWIST1', 'RUNX2']
final_suffix = '.genes_0618.xls'
# for case in cases_list:
# 	filter_var_by_genename(gene_list, case + '.annotated.txt', case + '.temp.xls', 15)
# combine_ann_txt(cases_list, '', '.temp.xls', 'cases' + final_suffix)
# for control in control_list:
# 	filter_var_by_genename(gene_list, control + '.annotated.txt', control + '.temp.xls', 15)
# combine_ann_txt(control_list, '', '.temp.xls', 'controls' + final_suffix)


##jg get all vars in 3 genesets
cases_list = ['cohort1.cases.95445', 'cohort1.cases.95446', 'cohort1.cases.95447', 'cohort1.cases.95448', 'cohort1.cases.95450', 'cohort1.cases.95452', 'cohort1.cases.95453', 'cohort1.cases.95454', 'cohort1.cases.95456', 'cohort1.cases.95457', 'cohort1.cases.95458', 'cohort1.cases.95459', 'cohort1.cases.95460', 'cohort1.cases.95462', 'cohort1.cases.95463', 'cohort1.cases.95464', 'cohort1.cases.95466', 'cohort1.cases.95467', 'cohort1.cases.95468', 'cohort1.cases.95469', 'cohort1.cases.95470', 'cohort1.cases.95471', 'cohort1.cases.95472', 'cohort1.cases.95474', 'cohort1.cases.95475', 'cohort1.cases.95476', 'cohort1.cases.95477', 'cohort1.cases.95479', 'cohort1.cases.95481', 'cohort1.cases.95482', 'cohort1.cases.95484', 'cohort1.cases.95485', 'cohort1.cases.95486', 'cohort1.cases.95487', 'cohort1.cases.95488', 'cohort1.cases.95490', 'cohort1.cases.95491', 'cohort1.cases.95492', 'cohort1.cases.95493', 'cohort1.cases.95494', 'cohort1.cases.95495', 'cohort1.cases.95496', 'cohort1.cases.95498', 'cohort1.cases.95500', 'cohort1.cases.95501', 'cohort1.cases.95502', 'cohort1.cases.95505', 'cohort1.cases.95506', 'cohort1.cases.95507', 'cohort1.cases.95509', 'cohort1.cases.95510', 'cohort1.cases.95511', 'cohort1.cases.95512', 'cohort1.cases.95514', 'cohort1.cases.95515', 'cohort1.cases.95516', 'cohort1.cases.95518', 'cohort1.cases.95519', 'cohort1.cases.95520', 'cohort1.cases.95521', 'cohort1.cases.95522', 'cohort1.cases.95524', 'cohort1.cases.95525', 'cohort1.cases.95526', 'cohort1.cases.95527', 'cohort1.cases.95528', 'cohort1.cases.95529', 'cohort1.cases.95530', 'cohort1.cases.95531', 'cohort1.cases.95534', 'cohort1.cases.95535', 'cohort1.cases.95536', 'cohort1.cases.95537', 'cohort1.cases.95540', 'cohort1.cases.95542', 'cohort1.cases.95543', 'cohort1.cases.95544', 'cohort1.cases.95545', 'cohort1.cases.95546', 'cohort1.cases.95547', 'cohort1.cases.95549', 'cohort1.cases.95550', 'cohort1.cases.95551', 'cohort1.cases.95552', 'cohort1.cases.95555', 'cohort1.cases.95556', 'cohort1.cases.95557', 'cohort1.cases.95558', 'cohort1.cases.95559', 'cohort1.cases.95560', 'cohort1.cases.95562', 'cohort1.cases.95563', 'cohort1.cases.95564', 'cohort1.cases.95565', 'cohort1.cases.95566', 'cohort1.cases.95567', 'cohort1.cases.95569', 'cohort1.cases.95570', 'cohort1.cases.95571', 'cohort1.cases.95572', 'cohort1.cases.95573', 'cohort1.cases.95575', 'cohort1.cases.95576', 'cohort1.cases.95578', 'cohort1.cases.95579', 'cohort1.cases.95580', 'cohort1.cases.95581', 'cohort1.cases.95583', 'cohort1.cases.95584', 'cohort1.cases.95585', 'cohort1.cases.95586', 'cohort1.cases.95587', 'cohort1.cases.95588', 'cohort1.cases.95589', 'cohort1.cases.95591', 'cohort1.cases.95592', 'cohort1.cases.95593', 'cohort1.cases.95594', 'cohort1.cases.95595', 'cohort1.cases.95596', 'cohort1.cases.95597', 'cohort1.cases.95598', 'cohort1.cases.95599', 'cohort1.cases.95600', 'cohort1.cases.95601', 'cohort1.cases.95602', 'cohort1.cases.95603', 'cohort1.cases.95604', 'cohort1.cases.95605', 'cohort1.cases.95606', 'cohort1.cases.95609', 'cohort1.cases.95610', 'cohort1.cases.95611', 'cohort1.cases.95612', 'cohort1.cases.95614', 'cohort1.cases.95615', 'cohort1.cases.95616', 'cohort1.cases.95617', 'cohort1.cases.95619', 'cohort1.cases.95620', 'cohort1.cases.95622', 'cohort1.cases.95623', 'cohort1.cases.95624', 'cohort1.cases.95626', 'cohort1.cases.95627', 'cohort1.cases.95628', 'cohort1.cases.95629', 'cohort1.cases.95630', 'cohort1.cases.95631', 'cohort1.cases.95633', 'cohort1.cases.95634', 'cohort1.cases.95635', 'cohort1.cases.95636', 'cohort1.cases.95637', 'cohort1.cases.95638', 'cohort1.cases.95639', 'cohort1.cases.95640', 'cohort1.cases.95642', 'cohort1.cases.95644', 'cohort1.cases.95645', 'cohort1.cases.95646', 'cohort1.cases.95649', 'cohort1.cases.95650', 'cohort1.cases.95651', 'cohort1.cases.95652', 'cohort1.cases.95653', 'cohort1.cases.95654', 'cohort1.cases.95655', 'cohort1.cases.95656', 'cohort1.cases.95657', 'cohort1.cases.95658', 'cohort1.cases.95661', 'cohort1.cases.95662', 'cohort1.cases.95663', 'cohort1.cases.95665', 'cohort1.cases.95668', 'cohort1.cases.95669', 'cohort1.cases.95670', 'cohort1.cases.95671', 'cohort1.cases.95673', 'cohort1.cases.95675', 'cohort1.cases.95676', 'cohort1.cases.95677', 'cohort1.cases.95678', 'cohort1.cases.95679', 'cohort1.cases.95680', 'cohort1.cases.95681', 'cohort1.cases.95682', 'cohort1.cases.95683', 'cohort1.cases.95684', 'cohort1.cases.95685', 'cohort1.cases.95686', 'cohort1.cases.95688', 'cohort1.cases.95689', 'cohort1.cases.95690', 'cohort1.cases.95691', 'cohort1.cases.95693', 'cohort1.cases.95694', 'cohort1.cases.95696', 'cohort1.cases.95697', 'cohort1.cases.95698', 'cohort1.cases.95699', 'cohort1.cases.95700', 'cohort1.cases.95702', 'cohort1.cases.95703', 'cohort2.cases.163499', 'cohort2.cases.163502', 'cohort2.cases.163503', 'cohort2.cases.163504', 'cohort2.cases.163505', 'cohort2.cases.163506', 'cohort2.cases.163507', 'cohort2.cases.163508', 'cohort2.cases.163509', 'cohort2.cases.163510', 'cohort2.cases.163511', 'cohort2.cases.163513', 'cohort2.cases.163516', 'cohort2.cases.163517', 'cohort2.cases.163518', 'cohort2.cases.163519', 'cohort2.cases.163520', 'cohort2.cases.163521', 'cohort2.cases.163522', 'cohort2.cases.163523', 'cohort2.cases.163524', 'cohort2.cases.163525', 'cohort2.cases.163526', 'cohort2.cases.163527', 'cohort2.cases.163528', 'cohort2.cases.163529', 'cohort2.cases.163530', 'cohort2.cases.163531', 'cohort2.cases.163532', 'cohort2.cases.163533', 'cohort2.cases.163534', 'cohort2.cases.163535', 'cohort2.cases.163537', 'cohort2.cases.163538', 'cohort2.cases.163539', 'cohort2.cases.163540', 'cohort2.cases.163541', 'cohort2.cases.163543', 'cohort2.cases.163544', 'cohort2.cases.163545', 'cohort2.cases.163546', 'cohort2.cases.163547', 'cohort2.cases.163548', 'cohort2.cases.163549', 'cohort2.cases.163553', 'cohort2.cases.163555', 'cohort2.cases.163556', 'cohort2.cases.163557', 'cohort2.cases.163558', 'cohort2.cases.163559', 'cohort2.cases.163560', 'cohort2.cases.163561', 'cohort2.cases.163562', 'cohort2.cases.163563', 'cohort2.cases.163564', 'cohort2.cases.163565', 'cohort2.cases.163566', 'cohort2.cases.163567', 'cohort2.cases.163568', 'cohort2.cases.163569', 'cohort2.cases.163571', 'cohort2.cases.163572', 'cohort2.cases.163573', 'cohort2.cases.163575', 'cohort2.cases.163576', 'cohort2.cases.163577', 'cohort2.cases.163578', 'cohort2.cases.163579', 'cohort2.cases.163580', 'cohort2.cases.163581', 'cohort2.cases.163582', 'cohort2.cases.163584', 'cohort2.cases.163585', 'cohort2.cases.163586', 'cohort2.cases.163587', 'cohort2.cases.163588', 'cohort2.cases.163590', 'cohort2.cases.163591', 'cohort2.cases.163592', 'cohort2.cases.163593', 'cohort2.cases.163594', 'cohort2.cases.163595', 'cohort2.cases.163596', 'cohort2.cases.163597', 'cohort2.cases.163598', 'cohort2.cases.163601', 'cohort2.cases.163602', 'cohort2.cases.163603', 'cohort2.cases.163604', 'cohort2.cases.163605', 'cohort2.cases.163606', 'cohort2.cases.163607', 'cohort2.cases.163609', 'cohort2.cases.163610', 'cohort2.cases.163611', 'cohort2.cases.163614', 'cohort2.cases.163615', 'cohort2.cases.163616', 'cohort2.cases.163619', 'cohort2.cases.163620', 'cohort2.cases.163621', 'cohort2.cases.163622', 'cohort2.cases.163623', 'cohort2.cases.163624', 'cohort2.cases.163625', 'cohort2.cases.163626', 'cohort2.cases.163627', 'cohort2.cases.163628', 'cohort2.cases.163629', 'cohort2.cases.163630', 'cohort2.cases.163631', 'cohort2.cases.163632', 'cohort2.cases.163633', 'cohort2.cases.163634', 'cohort2.cases.163635', 'cohort2.cases.163636', 'cohort2.cases.163637', 'cohort2.cases.163639', 'cohort2.cases.163640', 'cohort2.cases.163641', 'cohort2.cases.163642', 'cohort2.cases.163643', 'cohort2.cases.163644', 'cohort2.cases.163645', 'cohort2.cases.163646', 'cohort2.cases.163647', 'cohort2.cases.163649', 'cohort2.cases.163651', 'cohort2.cases.163652', 'cohort2.cases.163653', 'cohort2.cases.163654', 'cohort2.cases.163655', 'cohort2.cases.163656', 'cohort2.cases.163657', 'cohort2.cases.163658', 'cohort2.cases.163659', 'cohort2.cases.163660', 'cohort2.cases.163661', 'cohort2.cases.163662', 'cohort2.cases.163663', 'cohort2.cases.163664', 'cohort2.cases.163665', 'cohort2.cases.163666', 'cohort2.cases.163667', 'cohort2.cases.163668', 'cohort2.cases.163669', 'cohort2.cases.163672', 'cohort2.cases.163674', 'cohort2.cases.163678', 'cohort2.cases.163679', 'cohort2.cases.163680', 'cohort2.cases.163681', 'cohort2.cases.163682', 'cohort2.cases.163683', 'cohort2.cases.163684', 'cohort2.cases.163685', 'cohort2.cases.163686', 'cohort2.cases.163687', 'cohort2.cases.163688', 'cohort2.cases.163690', 'cohort2.cases.163691', 'cohort2.cases.163692', 'cohort2.cases.163694', 'cohort2.cases.163695', 'cohort2.cases.163696', 'cohort2.cases.163697', 'cohort2.cases.163698', 'cohort2.cases.163699', 'cohort2.cases.163700', 'cohort2.cases.163701', 'cohort2.cases.163702', 'cohort2.cases.163703', 'cohort2.cases.163704', 'cohort2.cases.163706', 'cohort2.cases.163707', 'cohort2.cases.163709', 'cohort2.cases.163710', 'cohort2.cases.163711', 'cohort2.cases.163712', 'cohort2.cases.163713', 'cohort2.cases.163714', 'cohort2.cases.163717', 'cohort2.cases.163718', 'cohort2.cases.163719', 'cohort2.cases.163720', 'cohort2.cases.163721']
control_list = ['cohort1.ctls.95449', 'cohort1.ctls.95451', 'cohort1.ctls.95455', 'cohort1.ctls.95461', 'cohort1.ctls.95465', 'cohort1.ctls.95473', 'cohort1.ctls.95478', 'cohort1.ctls.95480', 'cohort1.ctls.95483', 'cohort1.ctls.95489', 'cohort1.ctls.95497', 'cohort1.ctls.95499', 'cohort1.ctls.95504', 'cohort1.ctls.95508', 'cohort1.ctls.95513', 'cohort1.ctls.95517', 'cohort1.ctls.95523', 'cohort1.ctls.95532', 'cohort1.ctls.95533', 'cohort1.ctls.95538', 'cohort1.ctls.95539', 'cohort1.ctls.95541', 'cohort1.ctls.95548', 'cohort1.ctls.95553', 'cohort1.ctls.95554', 'cohort1.ctls.95561', 'cohort1.ctls.95568', 'cohort1.ctls.95574', 'cohort1.ctls.95577', 'cohort1.ctls.95582', 'cohort1.ctls.95590', 'cohort1.ctls.95607', 'cohort1.ctls.95608', 'cohort1.ctls.95613', 'cohort1.ctls.95618', 'cohort1.ctls.95621', 'cohort1.ctls.95625', 'cohort1.ctls.95647', 'cohort1.ctls.95648', 'cohort1.ctls.95659', 'cohort1.ctls.95660', 'cohort1.ctls.95666', 'cohort1.ctls.95667', 'cohort1.ctls.95672', 'cohort1.ctls.95674', 'cohort1.ctls.95692', 'cohort1.ctls.95695', 'cohort1.ctls.95704', 'cohort2.ctls.163500', 'cohort2.ctls.163501', 'cohort2.ctls.163512', 'cohort2.ctls.163514', 'cohort2.ctls.163515', 'cohort2.ctls.163536', 'cohort2.ctls.163542', 'cohort2.ctls.163550', 'cohort2.ctls.163551', 'cohort2.ctls.163552', 'cohort2.ctls.163554', 'cohort2.ctls.163570', 'cohort2.ctls.163574', 'cohort2.ctls.163583', 'cohort2.ctls.163589', 'cohort2.ctls.163599', 'cohort2.ctls.163600', 'cohort2.ctls.163608', 'cohort2.ctls.163612', 'cohort2.ctls.163613', 'cohort2.ctls.163617', 'cohort2.ctls.163618', 'cohort2.ctls.163638', 'cohort2.ctls.163648', 'cohort2.ctls.163650', 'cohort2.ctls.163670', 'cohort2.ctls.163671', 'cohort2.ctls.163673', 'cohort2.ctls.163675', 'cohort2.ctls.163676', 'cohort2.ctls.163677', 'cohort2.ctls.163689', 'cohort2.ctls.163693', 'cohort2.ctls.163705', 'cohort2.ctls.163708', 'cohort2.ctls.163715', 'cohort2.ctls.163716']
'''
##set1
gene_list = ['SEMA3F', 'MCUB', 'ITGA3', 'BAIAP2L1', 'CYB561', 'DCN', 'SEMA3B', 'IGF1', 'NRXN3', 'TNFRSF1B', 'APBA2', 'DAPK2', 'RAB27B', 'CAPG', 'HDAC9', 'LAMA3', 'SPAG4', 'EYA2', 'SLC9A7', 'PRKCZ', 'TGFBR3', 'CYBRD1', 'MRVI1', 'CLCN4', 'JADE1', 'FBLN1', 'SLC1A3', 'TCF7', 'COBLL1', 'KCNK2', 'FER1L4', 'OAS1', 'SLC7A8', 'TGFB2', 'CDC6', 'MAP3K1', 'NEFH', 'HMOX1', 'SYNGR1', 'LGMN', 'BDKRB1', 'MID1', 'CRISPLD2', 'WISP1', 'NDRG1', 'FSD1', 'ZNF175', 'ISYNA1', 'PTN', 'TSPAN12', 'PCOLCE', 'TSPAN13', 'GLI3', 'C5', 'APBA1', 'ATRNL1', 'KAZALD1', 'SYNGR2', 'CCDC34', 'TSPAN11', 'KRT18', 'OAS3', 'COL12A1', 'NEDD9', 'KCTD20', 'RBM24', 'ARRDC3', 'PDE4D', 'FGF1', 'C3orf52', 'PODXL2', 'SLC4A3', 'ITGA4', 'IL1R1', 'ID2', 'PRRX1', 'RGS2', 'OLFML3', 'AL591845.1', 'RGS4', 'SLC2A1', 'MAN1C1', 'AVPI1', 'PLS1', 'TMPO', 'INHBA', 'TWIST1', 'SRGN', 'P4HA1', 'CDKN2C', 'ATP7B', 'NR4A1', 'DAW1', 'RUNX2', 'HS3ST3B1', 'GRIA3', 'FLRT3', 'LAMP5', 'AMOT', 'WNK4', 'SIX1', 'RHOJ', 'POM121L9P', 'DOCK4', 'SNRPN', 'SAT1', 'CNN1', 'ULBP2', 'POPDC3', 'MATN2', 'DCLK1', 'COL4A2', 'MDFIC', 'SLC26A10', 'ALPK3', 'HS6ST1', 'KLF4', 'ANGPTL2', 'ENPP2', 'LRRC32', 'UACA', 'SMAD6', 'IFI44L', 'PCDH10', 'BMPR1B', 'SHROOM3', 'KIAA1644', 'SLC38A4', 'AMIGO2', 'LUM', 'MFGE8', 'IMPA2', 'GALNT1', 'SLC16A3', 'SECTM1', 'PMAIP1', 'PADI1', 'EPHA2', 'CTSK', 'SELENBP1', 'WNT9A', 'STAC', 'BOC', 'OCIAD2', 'SCD5', 'SLC2A12', 'EGFR', 'RSPO2', 'IDNK', 'ANKRD1', 'CCDC81', 'ADAMTS12', 'DPYSL4', 'RNF144A', 'TMEM56', 'HSPB8', 'ZFP36L2', 'PDLIM3', 'AFF2', 'PCDH1', 'SH3RF2', 'FZD1', 'AFAP1L1', 'AUTS2', 'HIST1H4H', 'SLC45A3', 'RCAN1', 'C1R', 'ACE', 'MPZL3', 'CXCL16', 'FGF11', 'PLPP3', 'PRKAA2', 'AK4', 'FBLIM1', 'SNED1', 'FBXO41', 'ACTG2', 'SLC16A14', 'FZD5', 'GPR155', 'GPR85', 'C8orf34', 'PGM2L1', 'CLMP', 'DCHS1', 'SAMD14', 'PRRX2', 'PRRT2', 'FILIP1L', 'BDKRB2', 'MN1', 'LDB2', 'RAC3', 'PCDH7', 'ALCAM', 'CDH2', 'KRT19', 'PDE7B', 'MAP6', 'ENC1', 'CTPS1', 'SPTLC3', 'IL17D', 'SLFN11', 'RARG', 'EGFL7', 'MRGPRF', 'BNC2', 'ABLIM3', 'STOX2', 'OLR1', 'TNFRSF10D', 'CNTNAP2', 'FZD4', 'YPEL2', 'BNIP3', 'IRX5', 'KCTD12', 'C5orf46', 'PCED1B', 'GAS1', 'GREM2', 'OXTR', 'NQO1', 'PENK', 'NXPH4', 'NXPH3', 'COL18A1', 'TMEM119', 'OVCH2', 'KREMEN1', 'EFNA5', 'AP001062.1', 'PROS1', '43348', 'SLC24A3', 'FLRT2', 'NDUFA4L2', 'NTF3', 'C11orf87', 'CCDC190', 'INSIG1', 'COL4A1', 'SAMD11', 'COL4A5', 'NHS', 'C15orf52', 'DPYD', 'PCDH18', 'CACNA1H', 'WNK3', 'SCN8A', 'BNIP3P1', 'PDCD1LG2', 'C2orf27A', 'SH3BGRL2', 'SLC5A3', 'ALPK2', 'MFAP3L', 'GP1BB', 'LGR4', 'LCMT1', 'DUXAP8', 'STK38L', 'GPSM3', 'SH3D21', 'ARHGEF28', 'MTMR9LP', 'CEBPD', 'AL645608.1', 'AC012462.1', 'ELFN1', 'AC046143.1', 'CCDC144NL-AS1', 'RNF217-AS1', 'LINC00888', 'DUXAP10', 'NEAT1', 'MIR210HG', 'NR2F2-AS1', 'CARMN', 'AC105383.1', 'TMEM200B', 'PCDHGB7', 'LINC02407', 'MGAM', 'ITGB3', 'AL590004.4']
final_suffix = '.TWIST1_IGF1.xls'
for case in cases_list:
	filter_var_by_genename(gene_list, case + '.annotated.txt', case + '.temp.xls', 15)
combine_ann_txt(cases_list, '', '.temp.xls', 'cases' + final_suffix)

for control in control_list:
	filter_var_by_genename(gene_list, control + '.annotated.txt', control + '.temp.xls', 15)
combine_ann_txt(control_list, '', '.temp.xls', 'controls' + final_suffix)
##set2
gene_list = ['SEMA3F', 'TFPI', 'MCUB', 'ITGA3', 'CRLF1', 'CYB561', 'IGF1', 'NRXN3', 'DAPK2', 'DNAH5', 'RAB27B', 'CAPG', 'HDAC9', 'PTPRU', 'EYA2', 'NGEF', 'FGFR2', 'TGFBR3', 'MRVI1', 'CLCN4', 'MCAM', 'JADE1', 'EDN1', 'BRINP1', 'SLC1A3', 'TCF7', 'COBLL1', 'KCNK2', 'MMP2', 'ASAP3', 'FER1L4', 'MAP3K1', 'TIMP3', 'HMOX1', 'SYNGR1', 'LGMN', 'SYNDIG1', 'MID1', 'CRISPLD2', 'HOMER2', 'WISP1', 'NDRG1', 'SH2D4A', 'APLP1', 'HAS1', 'ISYNA1', 'PTN', 'TSPAN12', 'PCOLCE', 'TSPAN13', 'GLI3', 'APBA1', 'ATRNL1', 'SYNGR2', 'CCDC34', 'VWA5A', 'TSPAN11', 'ARHGDIB', 'COL12A1', 'NEDD9', 'KCTD20', 'BACH2', 'RBM24', 'SEMA5A', 'ARRDC3', 'FGF1', 'PODXL2', 'ADAM23', 'ITGA4', 'PLCD4', 'IL1R1', 'EPHA4', 'MARK1', 'RGS2', 'OLFML3', 'AL591845.1', 'RGS4', 'SGIP1', 'TMPO', 'GLT8D2', 'INHBA', 'TWIST1', 'SRGN', 'BICC1', 'CDKN2C', 'ATP7B', 'RUNX2', 'HS3ST3B1', 'FLRT3', 'AMOT', 'RHOJ', 'HSPA2', 'PLEKHG3', 'TNFRSF19', 'SNRPN', 'SAT1', 'CNN1', 'LSP1', 'SH3BP5', 'KHDRBS3', 'PODNL1', 'SERPINF1', 'NES', 'DCLK1', 'COL4A2', 'MDFIC', 'FAIM2', 'SLC26A10', 'THSD1', 'ALPK3', 'HS6ST1', 'ENPP2', 'LRRC32', 'TRPC6', 'UACA', 'PCDH10', 'BMPR1B', 'SHROOM3', 'KIAA1644', 'FGF7', 'MFGE8', 'PADI1', 'CRABP2', 'FLG', 'WNT9A', 'MARCH4', 'STAC', 'COL8A1', 'BOC', 'SCD5', 'PDGFC', 'EGFR', 'ITGB1BP2', 'RSPO2', 'ANKRD1', 'INA', 'ZNF214', 'CCDC81', 'ADAMTS12', 'RNF144A', 'TMEM56', 'HSPB8', 'TCF7L1', 'LURAP1L', 'PDLIM3', 'AFF2', 'KCNMA1', 'ADAMTSL3', 'PCDH1', 'NRG1', 'FZD1', 'AFAP1L1', 'AUTS2', 'HIST1H4H', 'RCAN1', 'IGF2BP1', 'ACE', 'CXCL16', 'FGF11', 'PLPP3', 'PRKAA2', 'AK4', 'FBLIM1', 'PDPN', 'SNED1', 'ACTG2', 'GPR155', 'GPR85', 'CTHRC1', 'C8orf34', 'PGM2L1', 'PACSIN3', 'CLMP', 'PLEKHF1', 'SAMD14', 'FILIP1L', 'BDKRB2', 'MN1', 'SH3TC2', 'BNC1', 'LDB2', 'PCDH7', 'TCAF2', 'KRT19', 'PDE7B', 'CTPS1', 'SPTLC3', 'IL16', 'RASGRP1', 'RARG', 'EGFL7', 'GXYLT2', 'BNC2', 'ABLIM3', 'STOX2', 'NDNF', 'OLR1', 'FZD4', 'YPEL2', 'CSRP2', 'IRX5', 'GCNT4', 'FAM20C', 'KCTD12', 'C5orf46', 'EGR3', 'ARSJ', 'GREM2', 'OXTR', 'NQO1', 'PENK', 'C1S', 'NXPH4', 'NXPH3', 'COL18A1', 'GPC6', 'TMEM119', 'OVCH2', 'TNFAIP8L3', 'AP001062.1', 'PROS1', '43348', 'CA13', 'NDUFA4L2', 'SHC4', 'CCDC190', 'INSIG1', 'CEACAM19', 'GCNT1', 'COL4A1', 'SAMD11', 'SPRY4', 'COL4A5', 'NHS', 'C15orf52', 'AKR1C3', 'CACNA1H', 'SULF2', 'WNK3', 'BNIP3P1', 'MAP3K5', 'C2orf27A', 'SH3BGRL2', 'SMOC1', 'SLC5A3', 'ZNF521', 'GRK5', 'MFAP3L', 'GP1BB', 'LGR4', 'DUXAP8', 'STK38L', 'GPSM3', 'SH3D21', 'PLXNA4', 'CEBPD', 'AL645608.1', 'AC012462.1', 'SFTA1P', 'ELFN1', 'AC046143.1', 'CCDC144NL-AS1', 'RNF217-AS1', 'LINC00888', 'HLA-L', 'IFITM10', 'DUXAP10', 'MIR210HG', 'CARMN', 'AC105383.1', 'TMEM200B', 'BBOX1-AS1', 'MGAM', 'AL590004.4', 'ADAMTS7P3']
final_suffix = '.RUNX2_IGF1.xls'
for case in cases_list:
	filter_var_by_genename(gene_list, case + '.annotated.txt', case + '.temp.xls', 15)
combine_ann_txt(cases_list, '', '.temp.xls', 'cases' + final_suffix)
for control in control_list:
	filter_var_by_genename(gene_list, control + '.annotated.txt', control + '.temp.xls', 15)
combine_ann_txt(control_list, '', '.temp.xls', 'controls' + final_suffix)
##set3
gene_list = ['COLEC12', 'FLRT3', 'EYA2', 'PREX1', 'STMN3', 'CAPG', 'IL1R1', 'HS6ST1', 'CCDC74B', 'SNED1', 'DAPK2', 'ANPEP', 'SLC1A3', 'ARRDC3', 'TCF7', 'BNC2', 'IRX3', 'IRX5', 'NQO1', 'CRISPLD2', 'DCLK1', 'MID1', 'NHS', 'BX119917.1', 'AMOT', 'AFF2', 'ELFN1', 'TWIST1', 'CDCA7L', 'GLI3', 'IGFBP3', 'EGFR', 'AUTS2', 'FZD1', 'MCM7', 'PCOLCE', 'MDFIC', 'PTN', 'INSIG1', 'AL645608.1', 'DHRS3', 'TMEM200B', 'CDKN2C', 'PLPP3', 'TGFBR3', 'TMEM56', 'OLFML3', 'KCNK2', 'WNT9A', 'GREM2', 'MN1', 'HMOX1', 'KIAA1644', 'TSPAN11', 'RARG', 'PRIM1', 'TMPO', 'IGF1', 'TMEM119', 'CXCL16', 'EFNB3', 'HS3ST3A1', 'HS3ST3B1', 'KRT19', 'ACE', 'MBOAT1', 'TCF19', 'RUNX2', 'COL12A1', 'LAMA4', 'PDE7B', 'BOC', 'LDB2', 'MCUB', 'HMGB2', 'BDKRB2', 'SCARA3', 'CEBPD', 'PENK', 'RSPO2', 'ENPP2', 'WISP1', 'NDRG1', 'SLC5A3', 'CBR3', 'CCDC34', 'LGR4', 'SLC43A3', 'CLMP', 'RNF144A', 'SAMD11', 'AL136452.1', 'CELF2', 'C10orf10', 'SRGN', 'ANKRD1', 'PDLIM1', 'CNNM2', 'ATRNL1', 'FOXS1', 'KLHL29', 'ACTG2', 'C2orf27A', 'COBLL1', 'SCN9A', 'GPR155', 'MAP2', 'AC012462.1', 'SNRPN', 'C15orf52', 'UACA', 'ALPK3', 'MFGE8', 'ZNF774', 'ARRDC4', 'ADAMTS12', 'MAP3K1', 'PCDHGA4', 'PCDH1', 'FGF1', 'KCTD16', 'C5orf46', 'ABLIM3', 'AFAP1L1', 'GPX3', 'CNN1', 'C5AR1', 'APBA1', 'EGFL7', 'ATP7B', 'COL4A1', 'COL4A2', 'SAT1', 'DMD', 'WNK3', 'COL4A5', 'PDGFA', 'TSPAN13', 'HDAC9', 'INHBA', 'TSPAN12', 'TMEM178B', 'MGAM', 'AGRN', 'FBLIM1', 'PADI1', 'SH3D21', 'PRKAA2', 'AK4', 'FAM212B-AS1', 'CCDC190', 'RGS4', 'DUXAP8', 'DGCR6', 'GP1BB', 'CCND2', 'PTPRO', 'STK38L', 'NXPH4', 'NDUFA4L2', 'NUAK1', 'HSPB8', 'CCDC144NL-AS1', 'NXPH3', 'ITGA3', 'YPEL2', 'SYNGR2', 'AL590004.4', 'NEDD9', 'RBM24', 'CMAHP', 'HIST1H2BK', 'DDR1', 'GPSM3', 'KCTD20', 'ADGRF5', 'SH3BGRL2', 'ANKRD6', 'OXTR', 'STAC', 'FILIP1L', 'RAB6B', 'LINC00888', 'AC046143.1', 'PCDH7', 'ATP10D', 'IGFBP7', 'SHROOM3', 'MYOZ2', 'AC108062.1', 'JADE1', 'AC105383.1', 'PCDH10', 'FNIP2', 'MFAP3L', 'PDLIM3', 'DUXAP9', 'DUXAP10', 'RHOJ', 'NRXN3', 'LGMN', 'RCAN1', 'AP001062.1', 'COL18A1', 'MIR210HG', 'MRVI1', 'LUZP2', 'PGM2L1', 'LRRC32', 'CCDC81', 'FZD4', 'RBM20', 'PCDHGA1', 'RAB33A', 'OLR1', 'FGF11', 'HIST1H2BD', 'HIST1H4H', 'RNF217-AS1', 'SLC6A1', 'POPDC2', 'STOX2', 'C8orf34', 'OVCH2', 'RAB27B', 'SPTLC3', 'KCTD12', 'CLCN4', 'GPR85', 'STRIP2', 'CYB561', 'SEMA3F', 'PROS1', 'BMPR1B', 'FER1L4', 'ITGA4', 'CARMN', 'ADGRL1', 'ISYNA1', 'HSPB6', 'CACNA1H', 'PDK3', 'AC004656.1', 'AL591845.1', 'CTPS1', 'RGS2', 'SEPT5', 'SYNGR1', 'ITGA7', 'SLC26A10', 'SAMD14', 'PODXL2', 'BST1', 'SCD5', 'BNIP3P1', 'PNMA2', 'RAB11FIP1', 'JAM2', 'ANO1', 'CITED4']
final_suffix = '.TWIST1_RUNX2.xls'
for case in cases_list:
	filter_var_by_genename(gene_list, case + '.annotated.txt', case + '.temp.xls', 15)
combine_ann_txt(cases_list, '', '.temp.xls', 'cases' + final_suffix)
for control in control_list:
	filter_var_by_genename(gene_list, control + '.annotated.txt', control + '.temp.xls', 15)
combine_ann_txt(control_list, '', '.temp.xls', 'controls' + final_suffix)
'''

##take filter vars and filter by gene list 0718
'''
gene_list = ['SEMA3F', 'MCUB', 'ITGA3', 'BAIAP2L1', 'CYB561', 'DCN', 'SEMA3B', 'IGF1', 'NRXN3', 'TNFRSF1B', 'APBA2', 'DAPK2', 'RAB27B', 'CAPG', 'HDAC9', 'LAMA3', 'SPAG4', 'EYA2', 'SLC9A7', 'PRKCZ', 'TGFBR3', 'CYBRD1', 'MRVI1', 'CLCN4', 'JADE1', 'FBLN1', 'SLC1A3', 'TCF7', 'COBLL1', 'KCNK2', 'FER1L4', 'OAS1', 'SLC7A8', 'TGFB2', 'CDC6', 'MAP3K1', 'NEFH', 'HMOX1', 'SYNGR1', 'LGMN', 'BDKRB1', 'MID1', 'CRISPLD2', 'WISP1', 'NDRG1', 'FSD1', 'ZNF175', 'ISYNA1', 'PTN', 'TSPAN12', 'PCOLCE', 'TSPAN13', 'GLI3', 'C5', 'APBA1', 'ATRNL1', 'KAZALD1', 'SYNGR2', 'CCDC34', 'TSPAN11', 'KRT18', 'OAS3', 'COL12A1', 'NEDD9', 'KCTD20', 'RBM24', 'ARRDC3', 'PDE4D', 'FGF1', 'C3orf52', 'PODXL2', 'SLC4A3', 'ITGA4', 'IL1R1', 'ID2', 'PRRX1', 'RGS2', 'OLFML3', 'AL591845.1', 'RGS4', 'SLC2A1', 'MAN1C1', 'AVPI1', 'PLS1', 'TMPO', 'INHBA', 'TWIST1', 'SRGN', 'P4HA1', 'CDKN2C', 'ATP7B', 'NR4A1', 'DAW1', 'RUNX2', 'HS3ST3B1', 'GRIA3', 'FLRT3', 'LAMP5', 'AMOT', 'WNK4', 'SIX1', 'RHOJ', 'POM121L9P', 'DOCK4', 'SNRPN', 'SAT1', 'CNN1', 'ULBP2', 'POPDC3', 'MATN2', 'DCLK1', 'COL4A2', 'MDFIC', 'SLC26A10', 'ALPK3', 'HS6ST1', 'KLF4', 'ANGPTL2', 'ENPP2', 'LRRC32', 'UACA', 'SMAD6', 'IFI44L', 'PCDH10', 'BMPR1B', 'SHROOM3', 'KIAA1644', 'SLC38A4', 'AMIGO2', 'LUM', 'MFGE8', 'IMPA2', 'GALNT1', 'SLC16A3', 'SECTM1', 'PMAIP1', 'PADI1', 'EPHA2', 'CTSK', 'SELENBP1', 'WNT9A', 'STAC', 'BOC', 'OCIAD2', 'SCD5', 'SLC2A12', 'EGFR', 'RSPO2', 'IDNK', 'ANKRD1', 'CCDC81', 'ADAMTS12', 'DPYSL4', 'RNF144A', 'TMEM56', 'HSPB8', 'ZFP36L2', 'PDLIM3', 'AFF2', 'PCDH1', 'SH3RF2', 'FZD1', 'AFAP1L1', 'AUTS2', 'HIST1H4H', 'SLC45A3', 'RCAN1', 'C1R', 'ACE', 'MPZL3', 'CXCL16', 'FGF11', 'PLPP3', 'PRKAA2', 'AK4', 'FBLIM1', 'SNED1', 'FBXO41', 'ACTG2', 'SLC16A14', 'FZD5', 'GPR155', 'GPR85', 'C8orf34', 'PGM2L1', 'CLMP', 'DCHS1', 'SAMD14', 'PRRX2', 'PRRT2', 'FILIP1L', 'BDKRB2', 'MN1', 'LDB2', 'RAC3', 'PCDH7', 'ALCAM', 'CDH2', 'KRT19', 'PDE7B', 'MAP6', 'ENC1', 'CTPS1', 'SPTLC3', 'IL17D', 'SLFN11', 'RARG', 'EGFL7', 'MRGPRF', 'BNC2', 'ABLIM3', 'STOX2', 'OLR1', 'TNFRSF10D', 'CNTNAP2', 'FZD4', 'YPEL2', 'BNIP3', 'IRX5', 'KCTD12', 'C5orf46', 'PCED1B', 'GAS1', 'GREM2', 'OXTR', 'NQO1', 'PENK', 'NXPH4', 'NXPH3', 'COL18A1', 'TMEM119', 'OVCH2', 'KREMEN1', 'EFNA5', 'AP001062.1', 'PROS1', '43348', 'SLC24A3', 'FLRT2', 'NDUFA4L2', 'NTF3', 'C11orf87', 'CCDC190', 'INSIG1', 'COL4A1', 'SAMD11', 'COL4A5', 'NHS', 'C15orf52', 'DPYD', 'PCDH18', 'CACNA1H', 'WNK3', 'SCN8A', 'BNIP3P1', 'PDCD1LG2', 'C2orf27A', 'SH3BGRL2', 'SLC5A3', 'ALPK2', 'MFAP3L', 'GP1BB', 'LGR4', 'LCMT1', 'DUXAP8', 'STK38L', 'GPSM3', 'SH3D21', 'ARHGEF28', 'MTMR9LP', 'CEBPD', 'AL645608.1', 'AC012462.1', 'ELFN1', 'AC046143.1', 'CCDC144NL-AS1', 'RNF217-AS1', 'LINC00888', 'DUXAP10', 'NEAT1', 'MIR210HG', 'NR2F2-AS1', 'CARMN', 'AC105383.1', 'TMEM200B', 'PCDHGB7', 'LINC02407', 'MGAM', 'ITGB3', 'AL590004.4']
final_suffix = '.TWIST1_IGF1.xls'
for freq in frequencies:
	for file_suffix in ['.damaging_all', '.damaging_any', '.protein_changing']:
		infile = 'cs_rnaseq_1116.c1_and_2.cases.' + str(freq) + file_suffix + '.xls'
		outfile = 'cs_rnaseq_1116.c1_and_2.cases.' + str(freq) + file_suffix + final_suffix
		filter_var_by_genename(gene_list, infile, outfile, 15)
		infile = 'cs_rnaseq_1116.c1_and_2.ctls.' + str(freq) + file_suffix + '.xls'
		outfile = 'cs_rnaseq_1116.c1_and_2.ctls.' + str(freq) + file_suffix + final_suffix
		filter_var_by_genename(gene_list, infile, outfile, 15)
gene_list = ['SEMA3F', 'TFPI', 'MCUB', 'ITGA3', 'CRLF1', 'CYB561', 'IGF1', 'NRXN3', 'DAPK2', 'DNAH5', 'RAB27B', 'CAPG', 'HDAC9', 'PTPRU', 'EYA2', 'NGEF', 'FGFR2', 'TGFBR3', 'MRVI1', 'CLCN4', 'MCAM', 'JADE1', 'EDN1', 'BRINP1', 'SLC1A3', 'TCF7', 'COBLL1', 'KCNK2', 'MMP2', 'ASAP3', 'FER1L4', 'MAP3K1', 'TIMP3', 'HMOX1', 'SYNGR1', 'LGMN', 'SYNDIG1', 'MID1', 'CRISPLD2', 'HOMER2', 'WISP1', 'NDRG1', 'SH2D4A', 'APLP1', 'HAS1', 'ISYNA1', 'PTN', 'TSPAN12', 'PCOLCE', 'TSPAN13', 'GLI3', 'APBA1', 'ATRNL1', 'SYNGR2', 'CCDC34', 'VWA5A', 'TSPAN11', 'ARHGDIB', 'COL12A1', 'NEDD9', 'KCTD20', 'BACH2', 'RBM24', 'SEMA5A', 'ARRDC3', 'FGF1', 'PODXL2', 'ADAM23', 'ITGA4', 'PLCD4', 'IL1R1', 'EPHA4', 'MARK1', 'RGS2', 'OLFML3', 'AL591845.1', 'RGS4', 'SGIP1', 'TMPO', 'GLT8D2', 'INHBA', 'TWIST1', 'SRGN', 'BICC1', 'CDKN2C', 'ATP7B', 'RUNX2', 'HS3ST3B1', 'FLRT3', 'AMOT', 'RHOJ', 'HSPA2', 'PLEKHG3', 'TNFRSF19', 'SNRPN', 'SAT1', 'CNN1', 'LSP1', 'SH3BP5', 'KHDRBS3', 'PODNL1', 'SERPINF1', 'NES', 'DCLK1', 'COL4A2', 'MDFIC', 'FAIM2', 'SLC26A10', 'THSD1', 'ALPK3', 'HS6ST1', 'ENPP2', 'LRRC32', 'TRPC6', 'UACA', 'PCDH10', 'BMPR1B', 'SHROOM3', 'KIAA1644', 'FGF7', 'MFGE8', 'PADI1', 'CRABP2', 'FLG', 'WNT9A', 'MARCH4', 'STAC', 'COL8A1', 'BOC', 'SCD5', 'PDGFC', 'EGFR', 'ITGB1BP2', 'RSPO2', 'ANKRD1', 'INA', 'ZNF214', 'CCDC81', 'ADAMTS12', 'RNF144A', 'TMEM56', 'HSPB8', 'TCF7L1', 'LURAP1L', 'PDLIM3', 'AFF2', 'KCNMA1', 'ADAMTSL3', 'PCDH1', 'NRG1', 'FZD1', 'AFAP1L1', 'AUTS2', 'HIST1H4H', 'RCAN1', 'IGF2BP1', 'ACE', 'CXCL16', 'FGF11', 'PLPP3', 'PRKAA2', 'AK4', 'FBLIM1', 'PDPN', 'SNED1', 'ACTG2', 'GPR155', 'GPR85', 'CTHRC1', 'C8orf34', 'PGM2L1', 'PACSIN3', 'CLMP', 'PLEKHF1', 'SAMD14', 'FILIP1L', 'BDKRB2', 'MN1', 'SH3TC2', 'BNC1', 'LDB2', 'PCDH7', 'TCAF2', 'KRT19', 'PDE7B', 'CTPS1', 'SPTLC3', 'IL16', 'RASGRP1', 'RARG', 'EGFL7', 'GXYLT2', 'BNC2', 'ABLIM3', 'STOX2', 'NDNF', 'OLR1', 'FZD4', 'YPEL2', 'CSRP2', 'IRX5', 'GCNT4', 'FAM20C', 'KCTD12', 'C5orf46', 'EGR3', 'ARSJ', 'GREM2', 'OXTR', 'NQO1', 'PENK', 'C1S', 'NXPH4', 'NXPH3', 'COL18A1', 'GPC6', 'TMEM119', 'OVCH2', 'TNFAIP8L3', 'AP001062.1', 'PROS1', '43348', 'CA13', 'NDUFA4L2', 'SHC4', 'CCDC190', 'INSIG1', 'CEACAM19', 'GCNT1', 'COL4A1', 'SAMD11', 'SPRY4', 'COL4A5', 'NHS', 'C15orf52', 'AKR1C3', 'CACNA1H', 'SULF2', 'WNK3', 'BNIP3P1', 'MAP3K5', 'C2orf27A', 'SH3BGRL2', 'SMOC1', 'SLC5A3', 'ZNF521', 'GRK5', 'MFAP3L', 'GP1BB', 'LGR4', 'DUXAP8', 'STK38L', 'GPSM3', 'SH3D21', 'PLXNA4', 'CEBPD', 'AL645608.1', 'AC012462.1', 'SFTA1P', 'ELFN1', 'AC046143.1', 'CCDC144NL-AS1', 'RNF217-AS1', 'LINC00888', 'HLA-L', 'IFITM10', 'DUXAP10', 'MIR210HG', 'CARMN', 'AC105383.1', 'TMEM200B', 'BBOX1-AS1', 'MGAM', 'AL590004.4', 'ADAMTS7P3']
final_suffix = '.RUNX2_IGF1.xls'
for freq in frequencies:
	for file_suffix in ['.damaging_all', '.damaging_any', '.protein_changing']:
		infile = 'cs_rnaseq_1116.c1_and_2.cases.' + str(freq) + file_suffix + '.xls'
		outfile = 'cs_rnaseq_1116.c1_and_2.cases.' + str(freq) + file_suffix + final_suffix
		filter_var_by_genename(gene_list, infile, outfile, 15)
		infile = 'cs_rnaseq_1116.c1_and_2.ctls.' + str(freq) + file_suffix + '.xls'
		outfile = 'cs_rnaseq_1116.c1_and_2.ctls.' + str(freq) + file_suffix + final_suffix
		filter_var_by_genename(gene_list, infile, outfile, 15)
gene_list = ['COLEC12', 'FLRT3', 'EYA2', 'PREX1', 'STMN3', 'CAPG', 'IL1R1', 'HS6ST1', 'CCDC74B', 'SNED1', 'DAPK2', 'ANPEP', 'SLC1A3', 'ARRDC3', 'TCF7', 'BNC2', 'IRX3', 'IRX5', 'NQO1', 'CRISPLD2', 'DCLK1', 'MID1', 'NHS', 'BX119917.1', 'AMOT', 'AFF2', 'ELFN1', 'TWIST1', 'CDCA7L', 'GLI3', 'IGFBP3', 'EGFR', 'AUTS2', 'FZD1', 'MCM7', 'PCOLCE', 'MDFIC', 'PTN', 'INSIG1', 'AL645608.1', 'DHRS3', 'TMEM200B', 'CDKN2C', 'PLPP3', 'TGFBR3', 'TMEM56', 'OLFML3', 'KCNK2', 'WNT9A', 'GREM2', 'MN1', 'HMOX1', 'KIAA1644', 'TSPAN11', 'RARG', 'PRIM1', 'TMPO', 'IGF1', 'TMEM119', 'CXCL16', 'EFNB3', 'HS3ST3A1', 'HS3ST3B1', 'KRT19', 'ACE', 'MBOAT1', 'TCF19', 'RUNX2', 'COL12A1', 'LAMA4', 'PDE7B', 'BOC', 'LDB2', 'MCUB', 'HMGB2', 'BDKRB2', 'SCARA3', 'CEBPD', 'PENK', 'RSPO2', 'ENPP2', 'WISP1', 'NDRG1', 'SLC5A3', 'CBR3', 'CCDC34', 'LGR4', 'SLC43A3', 'CLMP', 'RNF144A', 'SAMD11', 'AL136452.1', 'CELF2', 'C10orf10', 'SRGN', 'ANKRD1', 'PDLIM1', 'CNNM2', 'ATRNL1', 'FOXS1', 'KLHL29', 'ACTG2', 'C2orf27A', 'COBLL1', 'SCN9A', 'GPR155', 'MAP2', 'AC012462.1', 'SNRPN', 'C15orf52', 'UACA', 'ALPK3', 'MFGE8', 'ZNF774', 'ARRDC4', 'ADAMTS12', 'MAP3K1', 'PCDHGA4', 'PCDH1', 'FGF1', 'KCTD16', 'C5orf46', 'ABLIM3', 'AFAP1L1', 'GPX3', 'CNN1', 'C5AR1', 'APBA1', 'EGFL7', 'ATP7B', 'COL4A1', 'COL4A2', 'SAT1', 'DMD', 'WNK3', 'COL4A5', 'PDGFA', 'TSPAN13', 'HDAC9', 'INHBA', 'TSPAN12', 'TMEM178B', 'MGAM', 'AGRN', 'FBLIM1', 'PADI1', 'SH3D21', 'PRKAA2', 'AK4', 'FAM212B-AS1', 'CCDC190', 'RGS4', 'DUXAP8', 'DGCR6', 'GP1BB', 'CCND2', 'PTPRO', 'STK38L', 'NXPH4', 'NDUFA4L2', 'NUAK1', 'HSPB8', 'CCDC144NL-AS1', 'NXPH3', 'ITGA3', 'YPEL2', 'SYNGR2', 'AL590004.4', 'NEDD9', 'RBM24', 'CMAHP', 'HIST1H2BK', 'DDR1', 'GPSM3', 'KCTD20', 'ADGRF5', 'SH3BGRL2', 'ANKRD6', 'OXTR', 'STAC', 'FILIP1L', 'RAB6B', 'LINC00888', 'AC046143.1', 'PCDH7', 'ATP10D', 'IGFBP7', 'SHROOM3', 'MYOZ2', 'AC108062.1', 'JADE1', 'AC105383.1', 'PCDH10', 'FNIP2', 'MFAP3L', 'PDLIM3', 'DUXAP9', 'DUXAP10', 'RHOJ', 'NRXN3', 'LGMN', 'RCAN1', 'AP001062.1', 'COL18A1', 'MIR210HG', 'MRVI1', 'LUZP2', 'PGM2L1', 'LRRC32', 'CCDC81', 'FZD4', 'RBM20', 'PCDHGA1', 'RAB33A', 'OLR1', 'FGF11', 'HIST1H2BD', 'HIST1H4H', 'RNF217-AS1', 'SLC6A1', 'POPDC2', 'STOX2', 'C8orf34', 'OVCH2', 'RAB27B', 'SPTLC3', 'KCTD12', 'CLCN4', 'GPR85', 'STRIP2', 'CYB561', 'SEMA3F', 'PROS1', 'BMPR1B', 'FER1L4', 'ITGA4', 'CARMN', 'ADGRL1', 'ISYNA1', 'HSPB6', 'CACNA1H', 'PDK3', 'AC004656.1', 'AL591845.1', 'CTPS1', 'RGS2', 'SEPT5', 'SYNGR1', 'ITGA7', 'SLC26A10', 'SAMD14', 'PODXL2', 'BST1', 'SCD5', 'BNIP3P1', 'PNMA2', 'RAB11FIP1', 'JAM2', 'ANO1', 'CITED4']
final_suffix = '.TWIST1_RUNX2.xls'
for freq in frequencies:
	for file_suffix in ['.damaging_all', '.damaging_any', '.protein_changing']:
		infile = 'cs_rnaseq_1116.c1_and_2.cases.' + str(freq) + file_suffix + '.xls'
		outfile = 'cs_rnaseq_1116.c1_and_2.cases.' + str(freq) + file_suffix + final_suffix
		filter_var_by_genename(gene_list, infile, outfile, 15)
		infile = 'cs_rnaseq_1116.c1_and_2.ctls.' + str(freq) + file_suffix + '.xls'
		outfile = 'cs_rnaseq_1116.c1_and_2.ctls.' + str(freq) + file_suffix + final_suffix
		filter_var_by_genename(gene_list, infile, outfile, 15)
'''


final_suffix = '.ADAM_0219.xls'
gene_list = ['ADAM10', 'ADAM11', 'ADAM11,GJC1', 'ADAM11\x3bADAM11', 'ADAM12', 'ADAM12,LINC00601', 'ADAM15', 'ADAM15\x3bADAM15', 'ADAM17', 'ADAM17,IAH1', 'ADAM17,YWHAQ', 'ADAM18', 'ADAM19', 'ADAM1A', 'ADAM2', 'ADAM20', 'ADAM20,MED6', 'ADAM20P1', 'ADAM21', 'ADAM21P1', 'ADAM22', 'ADAM22\x3bADAM22', 'ADAM23', 'ADAM28', 'ADAM29', 'ADAM30', 'ADAM30,NOTCH2', 'ADAM32', 'ADAM33', 'ADAM33\x3bADAM33', 'ADAM3A', 'ADAM6', 'ADAM6,LINC00226', 'ADAM7', 'ADAM8', 'ADAM8\x3bADAM8', 'ADAM9', 'ADAM9,TM2D2', 'ADAMDEC1', 'ADAMTS1', 'ADAMTS10', 'ADAMTS12', 'ADAMTS12\x3bADAMTS12', 'ADAMTS13', 'ADAMTS14', 'ADAMTS14\x3bADAMTS14', 'ADAMTS15', 'ADAMTS15,MIR8052', 'ADAMTS16', 'ADAMTS16,ICE1', 'ADAMTS17', 'ADAMTS18', 'ADAMTS19', 'ADAMTS19-AS1', 'ADAMTS19-AS1\x3bADAMTS19', 'ADAMTS1,ADAMTS5', 'ADAMTS2', 'ADAMTS20', 'ADAMTS20,PUS7L', 'ADAMTS2,RUFY1', 'ADAMTS3', 'ADAMTS3,COX18', 'ADAMTS4', 'ADAMTS5', 'ADAMTS5,LINC00113', 'ADAMTS6', 'ADAMTS6,CENPK', 'ADAMTS7', 'ADAMTS7,MORF4L1', 'ADAMTS7P1', 'ADAMTS7P1,GOLGA6L17P', 'ADAMTS8', 'ADAMTS9', 'ADAMTS9-AS1', 'ADAMTS9-AS1\x3bADAMTS9', 'ADAMTS9-AS2', 'ADAMTS9-AS2,MAGI1', 'ADAMTSL1', 'ADAMTSL1\x3bADAMTSL1', 'ADAMTSL2', 'ADAMTSL3', 'ADAMTSL3,EFTUD1P1', 'ADAMTSL4', 'ADAMTSL4-AS1', 'ADAMTSL5', 'CACFD1\x3bADAMTS13', 'CWC27,ADAMTS6', 'FALEC,ADAMTSL4', 'FAM30A,ADAM6', 'FANK1,ADAM12', 'LIPC,ADAM10', 'LOC105369739,ADAMTS20', 'LOC400464,ADAMTS17', 'LOC646938,ADAMTS7', 'MAPKAPK5,ADAM1A', 'PRICKLE2,ADAMTS9', 'SYCE1L,ADAMTS18', 'ZBTB44,ADAMTS8', 'ZNF354C,ADAMTS2']
# for case in cases_list:
# 	filter_var_by_genename(gene_list, case + '.annotated.txt', case + '.temp.xls', 15)
# combine_ann_txt(cases_list, '', '.temp.xls', 'cases' + final_suffix)
# for control in control_list:
# 	filter_var_by_genename(gene_list, control + '.annotated.txt', control + '.temp.xls', 15)
# combine_ann_txt(control_list, '', '.temp.xls', 'controls' + final_suffix)

##reduced set of cases and controls -- get all vars and then filter
'''
new_list = ['cohort1.cases.95445', 'cohort1.cases.95471', 'cohort1.cases.95476', 'cohort1.cases.95482', 'cohort1.cases.95484', 'cohort1.cases.95495', 'cohort1.cases.95509', 'cohort1.cases.95547', 'cohort1.cases.95558', 'cohort1.cases.95563', 'cohort1.cases.95565', 'cohort1.cases.95571', 'cohort1.cases.95609', 'cohort1.cases.95684', 'cohort1.cases.95699', 'cohort2.cases.163507', 'cohort2.cases.163509', 'cohort2.cases.163517', 'cohort2.cases.163522', 'cohort2.cases.163524', 'cohort2.cases.163525', 'cohort2.cases.163539', 'cohort2.cases.163544', 'cohort2.cases.163548', 'cohort2.cases.163558', 'cohort2.cases.163559', 'cohort2.cases.163566', 'cohort2.cases.163575', 'cohort2.cases.163585', 'cohort2.cases.163588', 'cohort2.cases.163596', 'cohort2.cases.163598', 'cohort2.cases.163603', 'cohort2.cases.163607', 'cohort2.cases.163611', 'cohort2.cases.163620', 'cohort2.cases.163621', 'cohort2.cases.163625', 'cohort2.cases.163627', 'cohort2.cases.163630', 'cohort2.cases.163635', 'cohort2.cases.163637', 'cohort2.cases.163640', 'cohort2.cases.163654', 'cohort2.cases.163656', 'cohort2.cases.163658', 'cohort2.cases.163660', 'cohort2.cases.163663', 'cohort2.cases.163664', 'cohort2.cases.163678', 'cohort2.cases.163692', 'cohort2.cases.163695', 'cohort2.cases.163700', 'cohort2.cases.163706', 'cohort2.cases.163710', 'cohort2.cases.163712', 'cohort2.cases.163718', 'cohort1.ctls.95449', 'cohort1.ctls.95461', 'cohort1.ctls.95480', 'cohort1.ctls.95499', 'cohort1.ctls.95541', 'cohort1.ctls.95666', 'cohort1.ctls.95695', 'cohort2.ctls.163612', 'cohort2.ctls.163613']
combine_ann_txt(new_list, '', '.annotated.txt', 'subset_samples_0219.annotated.txt')
##get all protein changing
filtering_cols_new =[16, 64, 19]
filtering_annotated.filter(working_dir, "or", 'subset_samples_0219.annotated.txt', "subset_samples_0219.exonic.xls", [filtering_cols_new[0], filtering_cols_new[0]], ['==','=='], [exonic_definitions[0],exonic_definitions[1]])
filtering_annotated.filter(working_dir, "or", 'subset_samples_0219.exonic.xls', "subset_samples_0219.protein_changing.xls", make_list(filtering_cols_new[2], protein_changing_definitions), make_list('==', protein_changing_definitions), protein_changing_definitions)
'''

##different set of cases and controls (022719) -- get all vars and then filter
'''
new_list = ['cohort1.cases.95665', 'cohort1.cases.95663', 'cohort1.cases.95610', 'cohort1.cases.95462', 'cohort1.ctls.95451', 'cohort1.cases.95468', 'cohort1.cases.95570', 'cohort1.cases.95498', 'cohort1.cases.95694', 'cohort1.cases.95572', 'cohort1.cases.95565', 'cohort1.ctls.95538', 'cohort1.cases.95454', 'cohort1.ctls.95647', 'cohort1.cases.95676', 'cohort1.cases.95668', 'cohort1.cases.95696', 'cohort1.cases.95486', 'cohort1.cases.95446', 'cohort1.cases.95697', 'cohort1.cases.95502', 'cohort1.ctls.95590', 'cohort1.cases.95475', 'cohort1.cases.95579', 'cohort1.cases.95463', 'cohort1.cases.95680', 'cohort1.cases.95495', 'cohort1.cases.95491', 'cohort1.cases.95534', 'cohort1.cases.95622', 'cohort1.cases.95587', 'cohort1.cases.95481', 'cohort1.ctls.95621', 'cohort1.cases.95490', 'cohort1.ctls.95497', 'cohort1.cases.95469', 'cohort1.ctls.95582', 'cohort1.ctls.95449', 'cohort1.cases.95447', 'cohort1.ctls.95704', 'cohort1.cases.95703', 'cohort1.cases.95645', 'cohort1.cases.95634', 'cohort1.ctls.95574', 'cohort1.cases.95524', 'cohort1.cases.95594']
combine_ann_txt(new_list, '', '.annotated.txt', 'subset_samples_2_0219.annotated.txt')
##get all protein changing
filtering_cols_new =[16, 64, 19]
filtering_annotated.filter(working_dir, "or", 'subset_samples_2_0219.annotated.txt', "subset_samples_2_0219.exonic.xls", [filtering_cols_new[0], filtering_cols_new[0]], ['==','=='], [exonic_definitions[0],exonic_definitions[1]])
filtering_annotated.filter(working_dir, "or", 'subset_samples_2_0219.exonic.xls', "subset_samples_2_0219.protein_changing.xls", make_list(filtering_cols_new[2], protein_changing_definitions), make_list('==', protein_changing_definitions), protein_changing_definitions)
'''

##another set
new_list = ['cohort1.cases.95446', 'cohort1.cases.95447', 'cohort1.ctls.95449', 'cohort1.ctls.95451', 'cohort1.cases.95454', 'cohort1.cases.95462', 'cohort1.cases.95463', 'cohort1.cases.95468', 'cohort1.cases.95469', 'cohort1.cases.95475', 'cohort1.cases.95481', 'cohort1.cases.95486', 'cohort1.cases.95490', 'cohort1.cases.95491', 'cohort1.cases.95495', 'cohort1.ctls.95497', 'cohort1.cases.95498', 'cohort1.cases.95502', 'cohort1.cases.95524', 'cohort1.cases.95534', 'cohort1.ctls.95538', 'cohort1.cases.95565', 'cohort1.cases.95570', 'cohort1.cases.95572', 'cohort1.ctls.95574', 'cohort1.cases.95579', 'cohort1.ctls.95582', 'cohort1.cases.95587', 'cohort1.ctls.95590', 'cohort1.cases.95594', 'cohort1.cases.95610', 'cohort1.ctls.95621', 'cohort1.cases.95622', 'cohort1.cases.95634', 'cohort1.cases.95645', 'cohort1.ctls.95647', 'cohort1.cases.95663', 'cohort1.cases.95665', 'cohort1.cases.95668', 'cohort1.cases.95676', 'cohort1.cases.95680', 'cohort1.cases.95694', 'cohort1.cases.95696', 'cohort1.cases.95697', 'cohort1.cases.95703', 'cohort1.ctls.95704', 'cohort2.cases.163507', 'cohort2.cases.163509', 'cohort2.cases.163522', 'cohort2.cases.163524', 'cohort2.cases.163525', 'cohort2.cases.163539', 'cohort2.cases.163548', 'cohort2.cases.163559', 'cohort2.cases.163566', 'cohort2.cases.163585', 'cohort2.cases.163588', 'cohort2.cases.163596', 'cohort2.cases.163598', 'cohort2.cases.163607', 'cohort2.cases.163611', 'cohort2.cases.163620', 'cohort2.cases.163621', 'cohort2.cases.163625', 'cohort2.cases.163635', 'cohort2.cases.163637', 'cohort2.cases.163640', 'cohort2.cases.163654', 'cohort2.cases.163656', 'cohort2.cases.163658', 'cohort2.cases.163660', 'cohort2.cases.163663', 'cohort2.cases.163678', 'cohort2.cases.163692', 'cohort2.cases.163700', 'cohort2.cases.163706', 'cohort2.cases.163710', 'cohort2.cases.163712', 'cohort2.cases.163718']
combine_ann_txt(new_list, '', '.annotated.txt', 'subset_samples_3_0219.annotated.txt')
##get all protein changing
filtering_cols_new =[16, 64, 19]
filtering_annotated.filter(working_dir, "or", 'subset_samples_3_0219.annotated.txt', "subset_samples_3_0219.exonic.xls", [filtering_cols_new[0], filtering_cols_new[0]], ['==','=='], [exonic_definitions[0],exonic_definitions[1]])
filtering_annotated.filter(working_dir, "or", 'subset_samples_3_0219.exonic.xls', "subset_samples_3_0219.protein_changing.xls", make_list(filtering_cols_new[2], protein_changing_definitions), make_list('==', protein_changing_definitions), protein_changing_definitions)












