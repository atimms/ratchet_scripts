#!/usr/bin/env python
import sys
import subprocess
import os
import glob

##parameters
delim = '\t'
working_dir = '/data/atimms/cox_frg_aj_1217'
os.chdir(working_dir)

##files and programs
picard = '/home/atimms/programs/picard-tools-1.139/picard.jar'
aj_fasta = 'A_J.chromosomes.unplaced.gt2k.fa'
mm10_fasta = '/data/atimms/references/mm10/mm10.fa'
fasta_dict = {'aj_fa':aj_fasta, 'mm10_fa':mm10_fasta}
svaba = '/home/atimms/programs/svaba/bin/svaba'
##methods
def convert_bam_fastq_bedtools(bamfile, r1_fastq, r2_fastq):
	read_sorted_bam = bamfile[:-4] + 'n_sorted.bam'
	st_n_sort = subprocess.Popen(['samtools', 'sort', '-nO', 'bam', '-o', read_sorted_bam, '-@', '10', '-m', '10G', '-T', 'tempy', bamfile])
	st_n_sort.wait()
	bam_fq = subprocess.Popen(['bedtools', 'bamtofastq', '-i', read_sorted_bam, '-fq', r1_fastq, '-fq2', r2_fastq])
	bam_fq.wait()

def align_with_bwa_one_at_time(sample, r1_fq, r2_fq, fasta_name):
	rg = '@RG\\tID:' + sample + '_RG\\tSM:' + sample + "\\tPL:ILLUMINA"
	post_bwa_bam = sample + '.bwa.bam'
	sort_bam = sample + '.bwa_sort.bam'
	mkdup_bam = sample + '.' + fasta_name + '.bwa_mkdup.bam'
	fasta = fasta_dict[fasta_name]
	##bwa and convert to bam
	bwa_pe = subprocess.Popen(['bwa', 'mem', '-M', '-t', '30', '-R', rg, fasta, r1_fq, r2_fq], stdout=subprocess.PIPE)
	st_sam_bam_pe = subprocess.Popen(['samtools', 'view', '-b', '-@', '5', '-o', post_bwa_bam, '-'], stdin=bwa_pe.stdout)
	st_sam_bam_pe.wait()
	##sort bam
	st_sort_pe = subprocess.Popen(['samtools', 'sort', '-O', 'bam', '-T', sample, '-o', sort_bam, '-@', '10', '-m', '10G', post_bwa_bam])
	st_sort_pe.wait()
	##mark duplicates
	picard_md = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=LENIENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + mkdup_bam, 'METRICS_FILE=' + sample + '.metrics'])
	picard_md.wait()

def run_svaba_single_sample(bam, fasta, out_name):
	#svaba run -t $GERMLINE_BAM -p $CORES -L 6 -I -a germline_run -G $REF
	st_sort_pe = subprocess.Popen([svaba, 'run', '-t', bam, '-p', '20', '-I', '-a', out_name, '-G', fasta])
	st_sort_pe.wait()

def run_svaba_paired(t_bam, n_bam, fasta, out_name, region):
	#svaba run -t $GERMLINE_BAM -p $CORES -L 6 -I -a germline_run -G $REF
	st_sort_pe = subprocess.Popen([svaba, 'run', '-t', t_bam,'-n', n_bam, '-p', '20', '-I', '-a', out_name, '-G', fasta, '-k', region])
	st_sort_pe.wait()

def make_refgene_text_info_bed(in_text, out_bed):
	with open(in_text, "r") as infh, open(out_bed, "w") as outfh:
		lc = 0
		for line in infh:
			lc += 1
			if lc > 1:
				line = line.split(delim)
				lineout = [line[2]] + line[4:6] + [line[12]] + ['\n']
				print lineout
				outfh.write(delim.join(lineout))


def add_genename_to_vcf(vcfs, refgene):
	for vcf in vcfs:
		mut_type = vcf.split('.')[-2]
		analysis = vcf.split('.')[0]
		print mut_type
		##make bed file for bedtools intersect and temp text file
		var_bed = vcf.rsplit('.', 1)[0] + '.bed'
		int_bed = vcf.rsplit('.', 1)[0] + '.refgene_int.bed'
		temp_txt = vcf.rsplit('.', 1)[0] + '.temp.txt'
		final_text = vcf.rsplit('.', 1)[0] + '.with_genes.xls'
		with open(vcf, "r") as infh, open(var_bed, "w") as outbfh, open(temp_txt, "w") as outtfh:
			for line in infh:
				if line[:2] !='##':
					line = line.rstrip().split(delim)
					if line[0][0] == '#':
						outtfh.write(delim.join(line + ['chrom', 'start', 'end', 'genes', '\n']))
					else:
						
						chrom = line[0]
						start = str(int(line[1]) - 1)
						if mut_type == 'indel':
							end = line[1]
						elif mut_type == 'sv':
							second_info = line[4].replace('[',':').replace(']',':').split(':')
							second_chr = second_info[1]
							end = second_info[2]
							if second_chr == 'chr13':
								end = second_info[2]
								if int(start) > int(end) and chrom == 'chr13':
									# print start
									newend = start
									start = end
									end = newend
							else:
								##if on different chrom just use start
								end = line[1]

								# print line
								# print second_info
						if chrom == 'chr13':
							outbfh.write(delim.join([chrom, start, end, '\n']))
							outtfh.write(delim.join(line + [chrom, start, end, '\n']))
		##run bedtools intersect
		with open(int_bed, "w") as intbed_fh: 
			bt_mw = subprocess.Popen(['bedtools', 'intersect', '-a', var_bed, '-b', refgene, '-wao'], stdout=intbed_fh)
			bt_mw.wait()
		##make dict of results
		region_gene_dict = {}
		with open(int_bed, "r") as intbed_fh:
			for line in intbed_fh:
				line = line.rstrip().split(delim)
				region = '_'.join(line[:3])
				gene = line[6]
				# print region, gene
				if region in region_gene_dict:
					if gene not in region_gene_dict[region]:
						region_gene_dict[region].append(gene)
				else:
					region_gene_dict[region] = [gene]
		##then add gene name to file 
		for r in region_gene_dict:
			print r, region_gene_dict[r]
		with open(final_text, "w") as final_fh, open(temp_txt, "r") as temp_fh:
			lc = 0
			for line in temp_fh:
				lc += 1
				if lc == 1:
					final_fh.write(line)
				else:
					line = line.rstrip().split(delim)
					if analysis == 'frg':
						reg = '_'.join(line[10:13])
					else:
						reg = '_'.join(line[11:14])
					genes = '_'.join(region_gene_dict[reg])
					final_fh.write(delim.join(line + [genes, '\n']))




##run methods
aj_original_bam = 'AJ_map_region.bam'
frg_original_bam = 'mm9_B3_fog_13_38855131_53042973.rmdup.bam'
aj = 'AJ'
frg = 'frg'
sample_dict = {aj:aj_original_bam, frg:frg_original_bam}

##make bam file for aj and frg
for sample in sample_dict:
	bam = sample_dict[sample]
	r1_fq = sample + '.r1.fastq'
	r2_fq = sample + '.r2.fastq'
	# convert_bam_fastq_bedtools(bam, r1_fq, r2_fq)
	# align_with_bwa_one_at_time(sample, r1_fq, r2_fq, 'aj_fa')
	# align_with_bwa_one_at_time(sample, r1_fq, r2_fq, 'mm10_fa')

##run svaba
# run_svaba_single_sample('frg.aj_fa.bwa_mkdup.bam', aj_fasta, 'frg.aj_fa')
# run_svaba_paired('frg.mm10_fa.bwa_mkdup.bam', 'AJ.mm10_fa.bwa_mkdup.bam', mm10_fasta, 'frg_vs_aj', 'chr13:38855131-53042973')

##add genes to vcf
vcfs_to_add_gene_names = ['frg.aj_fa.svaba.indel.vcf', 'frg.aj_fa.svaba.sv.vcf', 'frg.aj_fa.svaba.unfiltered.indel.vcf', 'frg.aj_fa.svaba.unfiltered.sv.vcf', 'frg_vs_aj.svaba.germline.indel.vcf', 'frg_vs_aj.svaba.germline.sv.vcf', 'frg_vs_aj.svaba.somatic.indel.vcf', 'frg_vs_aj.svaba.somatic.sv.vcf', 'frg_vs_aj.svaba.unfiltered.germline.indel.vcf', 'frg_vs_aj.svaba.unfiltered.germline.sv.vcf', 'frg_vs_aj.svaba.unfiltered.somatic.indel.vcf', 'frg_vs_aj.svaba.unfiltered.somatic.sv.vcf']
# vcfs_to_add_gene_names = ['frg.aj_fa.svaba.indel.vcf', 'frg.aj_fa.svaba.sv.vcf']
refgene_file = 'mm10_refGene.txt'
refgene_bed = 'mm10_refGene_tx.bed'

##make bed of all gens
# make_refgene_text_info_bed(refgene_file, refgene_bed)
##then reformat vcf file
add_genename_to_vcf(vcfs_to_add_gene_names, refgene_bed)





