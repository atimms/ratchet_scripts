#!/usr/bin/env python
import sys
import subprocess
import os
import glob


##note
'''
for mscentipede....
need python 2
load modules:
module load biobuilds/2017.11
module load local_python/2.7.14
install packages:
conda install -c conda-forge cvxopt
to plot profile need to alter the python scripts plot_accessibility_profile.py
add the lines before import matplotlib.pyplot as plot:
import matplotlib
matplotlib.use('agg')
'''


##parameters
delim = '\t'
working_dir = '/data/atimms/cherry_temp_1218'
os.chdir(working_dir)

##program
call_binding = '/home/atimms/programs/msCentipede/call_binding.py'
plt_accesability = '/home/atimms/programs/msCentipede/plot_accessibility_profile_adj.py'



def run_mscentipede(motif_files, bam_files):
	##python call_binding.py --task learn test/CTCF_chr10_motifs.txt.gz test/Gm12878_Rep1.bam test/Gm12878_Rep2.bam
	for motif_file in motif_files:
		print(motif_file, bam_files)
		# run_msc = subprocess.Popen(['python', call_binding, '--task', 'learn', '--protocol=ATAC_seq', motif_file] + bam_files)
		# run_msc.wait()
		run_msc_infer = subprocess.Popen(['python', call_binding, '--task', 'infer', '--protocol=ATAC_seq', motif_file] + bam_files)
		run_msc_infer.wait()
		# run_msc_plot = subprocess.Popen(['python', plt_accesability, '--protocol=ATAC_seq', motif_file])
		# run_msc_plot.wait()


# mscent_motifs = ['BORIS_all.homer_hg38_1018.txt.gz']
# mscent_motifs = ['BORIS.homer_hg38_1018.txt.gz']
# mscent_motifs = ['CRE.homer_hg38_1018.txt.gz']
mscent_motifs = ['CRX.homer_hg38_1018.txt.gz']
bams_to_analyze = glob.glob('*.bwa_mkdup_filtered.bam')

##run mscentipede
run_mscentipede(mscent_motifs, bams_to_analyze)		