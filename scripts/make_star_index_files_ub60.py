#!/usr/bin/env python
import sys
import subprocess
import os

##parameters
delim = '\t'
threads = '16'

##programs
star = '/home/atimms/programs/STAR/bin/Linux_x86_64/STAR'


##method
def make_star_index_files(star_genome_dir, genome_fas, genome_gtf, threads_to_use):
	star_index = subprocess.Popen([star, '--runMode', 'genomeGenerate', '--genomeDir', star_genome_dir, '--genomeFastaFiles', genome_fas, '--sjdbGTFfile', genome_gtf, '--runThreadN', threads_to_use])
	star_index.wait()






##run methods
genome_name = 'acomy1'
star_index_dir = '/data/atimms/references/star/' + genome_name
fa_file = '/data/atimms/references/igenomes/' + genome_name + '/genome.fa'
gtf_file = '/data/atimms/references/igenomes/' + genome_name +'/genes.gtf'

##make star index files for acomy 0618 (only need to do once)
make_star_index_files(star_index_dir, fa_file, gtf_file, threads)