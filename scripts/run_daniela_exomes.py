#!/tools/BioBuilds-2015.04/bin/python
import daniela_exome_pipeline_ub_v1


##call all methods, supply:
##working dir, file prefix (usually proband name), if starting from bam or fastq, sample dict, sample list (proband, dad, mom, sib(if available)

##templates
##fq, trio and single
#daniela_exome_pipeline_ub_v1.call_all_exome_methods(working_dir, '', 'fastq', {'':['', ''], '':['', ''], '':['', '']})
#daniela_exome_pipeline_ub_v1.call_all_exome_methods(working_dir, '', 'fastq', {'':['', '']})
##bam, trio and single
#daniela_exome_pipeline_ub_v1.call_all_exome_methods(working_dir, '', 'bam', {'':'', '':'', '':''})
#daniela_exome_pipeline_ub_v1.call_all_exome_methods(working_dir, '', 'bam', {'':''})


##batch1 -- ready
working_dir = '/data/atimms/microtia_exomes/batch1'
daniela_exome_pipeline_ub_v1.call_all_exome_methods(working_dir, '1010001', 'bam', {'101000101':'60017.bam'})
##daniela_exome_pipeline_ub_v1.call_all_exome_methods(working_dir, '1010002', 'bam', {'101000201':'60018.bam', '101000202':'60019.bam', '101000203':'60020.bam'})
##daniela_exome_pipeline_ub_v1.call_all_exome_methods(working_dir, '1010003', 'bam', {'101000301':'60021.bam', '101000302':'60022.bam'})
##daniela_exome_pipeline_ub_v1.call_all_exome_methods(working_dir, '1010004', 'bam', {'101000401':'60023.bam', '101000402':'60024.bam', '101000403':'60025.bam'})
##daniela_exome_pipeline_ub_v1.call_all_exome_methods(working_dir, '1010005', 'bam', {'101000501':'60026.bam', '101000502':'60027.bam', '101000503':'60028.bam'})
##daniela_exome_pipeline_ub_v1.call_all_exome_methods(working_dir, '1010006', 'bam', {'101000601':'60029.bam', '101000602':'60030.bam'})
##daniela_exome_pipeline_ub_v1.call_all_exome_methods(working_dir, '1010007', 'bam', {'101000701':'60031.bam', '101000702':'60032.bam', '101000703':'101000703.94638.bam'})
##daniela_exome_pipeline_ub_v1.call_all_exome_methods(working_dir, '1030001', 'bam', {'103000101':'60033.bam', '103000102':'60034.bam', '103000103':'60035.bam'})