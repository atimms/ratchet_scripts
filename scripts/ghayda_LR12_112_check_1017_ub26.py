#!/usr/bin/env python
import pybedtools as pbt
import subprocess
import os 


##parameters
delim = '\t'
working_dir = '/data/atimms/ghayda_ped_test_1017'
os.chdir(working_dir)

plink = '/home/atimms/programs/plink'
genename_nm = 'hg19_refGene_nm_and_name.txt'
exons_bed = 'hg19_refGene_coding_exons.plus2.no_chr.bed'

##calculate coverage from refgene genelist
##need to supply list of bams, file with nm\tgenename, exon bed file, list of genes, name of outfile
def calculate_exon_coverage_by_gene_from_genelist(bam_files, name_to_nm, exon_bed, genelist, output_file):
	##get bedfile for each gene
	nm_dict = {}        #dict contain list of nms found for each gene
	gene_bed_files = []  #list of all gene specific bed file
	for gene in genelist:
		##convert name to nms
		temp_bed = gene + '.temp.bed'
		nms = []
		with open(name_to_nm, "r") as name_file:
			for line in name_file:
				line = line.strip('\n').split(delim)
				if gene == line[1]:
					nms.append(line[0])
		print 'for gene %s we found the following %i NM names: %s'% (gene, len(nms), ','.join(nms))
		##get all exons that correspond to those nms and print to gene specific bed
		with open(exon_bed, "r") as ebed, open(temp_bed, "w") as tbed:
			exon_count = 0
			for line in ebed:
				line = line.split(delim)
				for nm in nms:
					bed_nm = line[3].split('_')[0] + '_' + line[3].split('_')[1]
					if nm == bed_nm:
						tbed.write(delim.join(line))
						exon_count += 1
		nm_dict[gene] = nms
		temp_merged_bed = gene + '.temp_merged.bed'
		print '%i total exons, which were sorted and merged and put in the file %s'% (exon_count, temp_merged_bed)
		##sort and merge bed file
		a = pbt.BedTool(temp_bed)
		b = a.sort().merge().moveto(temp_merged_bed)
		gene_bed_files.append(temp_merged_bed)
	print 'we have the bed files:', gene_bed_files
	
	##combine gene bed files
	combined_bed = 'all_temp.bed'
	comb_sorted_bed = 'all_sorted_temp.bed'
	with open(combined_bed, "w") as cbed:
		for gene_bed in gene_bed_files:
			with open(gene_bed, "r") as gbed:
				for line in gbed:
					cbed.write(line)
	c = pbt.BedTool(combined_bed)
	d = c.sort().merge().moveto(comb_sorted_bed)
	print 'bed files combined and then sorted in file:', comb_sorted_bed
			
	##parameters and header
	coverage_required = [1,5,10,20,50,100,200]
	cov_head = []
	for cr in coverage_required:
		cr1 = 'total bp >=' + str(cr)
		cr2 = 'percentage covered >=' + str(cr)
		cov_head.extend([cr1, cr2])
	header = delim.join(['gene', 'NMs', 'bam file', 'target size', 'total bp', 'coverage'] + cov_head + ['\n'])
	
	##open results file and write header
	with open(output_file, "w") as outfh:
		outfh.write(header)
		##go through bam files and get coverage histogram
		for bam in bam_files:
			print 'calculating coverge for bam file', bam
			##get temp coverage files
			e = pbt.BedTool(bam)
			f = pbt.BedTool(comb_sorted_bed)
			# f_cov = f.coverage(e, hist=True, sorted=True, output=bam + '.target_coverage_temp.txt')
			f_cov = f.coverage(e, hist=True, output=bam + '.target_coverage_temp.txt')

			# coverage_histogram = bam + '.target_coverage_temp.txt'
			# hist_fh = open(coverage_histogram, 'w')
			# output=bam + '.target_coverage_temp.txt'
			# bedtools_cov = subprocess.Popen(['bedtools', 'coverage', '-a', comb_sorted_bed, '-b', bam, '-hist', '-sorted'], stdout=hist_fh)
			# bedtools_cov.wait()
			# hist_fh.close()
			##calculate numbers from histogram for all genes
			with open(bam + '.target_coverage_temp.txt', "r") as cov_temp:
				seq_total = 0
				for line in cov_temp:
					line = line.strip('\n').split(delim)
					if line[0] == 'all':
						target_size = float(line[3])
						if line[1] != '0':
							seq_total += int(line[1]) *  int(line[2])
			cov_stats = []
			for cr in coverage_required:
				coverage_bp = 0
				with open(bam + '.target_coverage_temp.txt', "r")as cov_temp:
					for line in cov_temp:
						line = line.strip('\n').split(delim)
						if line[0] == 'all':
							target_size = float(line[3])
							if int(line[1]) >= cr:
								coverage_bp += int(line[2])
					cov_stats.extend([str(coverage_bp), str((coverage_bp / target_size) * 100)])
			line_out = delim.join(['all', '', bam, str(target_size), str(seq_total), str(seq_total / target_size)] + cov_stats + ['\n'])
			outfh.write(line_out)
		##go through data gene by gene
		for gene in genelist:
			##get list of exons and target size for gene
			exon_list = []
			target_size = 0
			gene_bed = gene + '.temp_merged.bed'
			with open(gene_bed, "r") as gbed:
				for exon in gbed:
					exon = exon.strip('\n').split(delim)
					exon_list.append(exon)
					target_size += (int(exon[2]) - float(exon[1]))
			# print gene, exon_list 
			##
			for bam in bam_files:
				seq_total = 0
				with open(bam + '.target_coverage_temp.txt', "r") as bam_cov:
					for line in bam_cov:
						line = line.strip('\n').split(delim)
						# print gene, line[:3], exon_list
						if line[:3] in exon_list:
							if line[4] != '0':
								seq_total += int(line[3]) *  int(line[4])
				if seq_total != 0:
					final_coverage =  str(seq_total / target_size)
				else:
					final_coverage = '0'
				cov_stats = []
				for cr in coverage_required:
					coverage_bp = 0
					with open(bam + '.target_coverage_temp.txt', "r") as bam_cov:
						for line in bam_cov:
							line = line.strip('\n').split(delim)
							if line[:3] in exon_list:
								if int(line[3]) >= cr:
									coverage_bp += int(line[4])
					#print bam, gene, cr, coverage_bp, target_size
					if coverage_bp != 0:
						cov_stats.extend([str(coverage_bp), str((coverage_bp / target_size) * 100)])
					else:
						cov_stats.extend([str(coverage_bp), '0'])

				line_out = delim.join([gene, ', '.join(nm_dict[gene]), bam, str(target_size), str(seq_total), final_coverage] + cov_stats + ['\n'])
				outfh.write(line_out)

def plink_roh(vcf, file_prefix):
	##bzgip vcf file
	if os.path.isfile(vcf):
		bgzip_cmd = subprocess.Popen(['bgzip', vcf])
		bgzip_cmd.wait()
	##correct filtering??
	bcftools_filter = subprocess.Popen(['bcftools', 'view', '-f', 'PASS', '-i', "QUAL>50 & MIN(DP)>50", '-V', 'indels', '-o', 'temp_plink.vcf.gz', '-O', 'z', vcf + '.gz'])
	bcftools_filter.wait()
	##generate plink file from vcf
	make_plink = subprocess.Popen([plink, '--vcf', 'temp_plink.vcf.gz', '--vcf-filter', '--const-fid', '--out', 'temp.pass_q50_dp50'])
	make_plink.wait()
	##prune for ld and then extract those variants
	plink_ld_prune = subprocess.Popen([plink, '--bfile', 'temp.pass_q50_dp50', '--indep-pairwise', '50', '5', '0.5'])
	plink_ld_prune.wait()
	plink_extract = subprocess.Popen([plink, '--bfile', 'temp.pass_q50_dp50', '--extract', 'plink.prune.in', '--make-bed', '--out', 'temp.pass_q50_dp50_prune'])
	plink_extract.wait()
	##check sex -results in .sexcheck file
	plink_roh = subprocess.Popen([plink, '--bfile', 'temp.pass_q50_dp50_prune', '--homozyg', '--out', file_prefix + '.pass_q50_dp50'])
	plink_roh.wait()

##calculate coverage from genes
bams = ['LR12-112a1.bwa_gatk.bam', 'LR12-112f.bwa_gatk.bam', 'LR12-112m.bwa_gatk.bam']
genelist = ['L1CAM', 'AP1S2', 'CCDC88C', 'B3GALTL', 'HYLS1', 'GPSM2', 'MPDZ', 'EML1', 'WDR81', 'AKT3', 'PIK3CA', 'PIK3R2', 'CCND2', 'MAP2K1', 'MAP2K2', 'NFIA']
# calculate_exon_coverage_by_gene_from_genelist(bams, genename_nm, exons_bed, genelist, 'LR12-112.gene_coverage.xls')

##roh in plink
vcf_file = '0002.vcf'
plink_roh(vcf_file, 'LR12-112')

