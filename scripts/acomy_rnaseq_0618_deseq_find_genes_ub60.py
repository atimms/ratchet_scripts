#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import shutil
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


'''
must load local python 3.5 to run
module load local_python/3.6.5
'''
##parameters
delim = '\t'
threads = '16'

##setup working directory where results will be
working_dir = '/data/atimms/acomy_rnaseq_0618'
os.chdir(working_dir)




##methods
def make_fasta_from_gtf(gtf, out_fasta, genome_fasta):
	gtf_dict = {}
	##mk dict from gtf, conatain chr, starts of exons, ends of exons and strand
	with open(gtf, "r") as gtfph:
		line_count = 0
		for line in gtfph:
			line = line.strip('\n').split(delim)
			el_type = line[2]
			gene_id = line[8].split(';')[0].split('"')[1]
			chrom = line[0]
			start = line[3]
			end = line[4]
			strand = line[6]
			if el_type == 'CDS':
				# print(gene_id, el_type)
				if gene_id in gtf_dict:
					if chrom != gtf_dict[gene_id][0]:
						print('multiple chrs for gene: ', gene_id)
					if strand != gtf_dict[gene_id][3]:
						print('different strandedness for gene: ', gene_id)
					gtf_dict[gene_id][1].append(start)
					gtf_dict[gene_id][2].append(end)
				else:
					gtf_dict[gene_id] = [chrom, [start], [end], strand]
	print(len(gtf_dict), 'entries in gtf dict')
	with open(out_fasta, "w") as out_fh:
		for g in gtf_dict:
			# print (g, gtf_dict[g])
			chrom = gtf_dict[g][0]
			starts = gtf_dict[g][1]
			ends = gtf_dict[g][2]
			strand = gtf_dict[g][3]
			temp_gene_fa = g + '_temp.fa'
			with open(temp_gene_fa, "w") as tbf_fh:
				##for each exon get the sequence
				for i, start in enumerate(starts):
					end = ends[i]
					cse = chrom + ':' + start + '-' + end
					# print(g, cse, i)
					st_faidx = subprocess.Popen(['samtools', 'faidx', genome_fasta, cse], stdout=tbf_fh)
					st_faidx.wait()
			##get the sequence from the fastq
			# '''
			sequence = ''
			with open(temp_gene_fa, "r") as tbf2_fh:
				for line in tbf2_fh:
					if line[0] != '>':
						sequence += line.rstrip()
						# print g, sequence, strand
			if strand == '-':
				out_sequence = SeqRecord(Seq(sequence, IUPAC.unambiguous_dna).reverse_complement(), id=g, description="")
			else:
				out_sequence = SeqRecord(Seq(sequence, IUPAC.unambiguous_dna), id=g, description="")
			##blast the sequence
			# result_handle = NCBIWWW.qblast("blastn", "nt", out_sequence.seq)
			# print out_sequence
			# print g, sequence, strand
			# print sequence.format("fasta")
			##write final fasta file
			SeqIO.write(out_sequence, out_fh, 'fasta')

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator) ##had to change this for python3
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

def split_fasta(in_fasta, out_fasta_prefix, no_of_records_req):
	filenames = []
	record_iter = SeqIO.parse(open(in_fasta),"fasta")
	for i, batch in enumerate(batch_iterator(record_iter, no_of_records_req)):
		filename = out_fasta_prefix + "%i.fasta" % (i + 1)
		filenames.append(filename)
		with open(filename, "w") as handle:
			count = SeqIO.write(batch, handle, "fasta")
			print("Wrote %i records to %s" % (count, filename))
	return filenames


def blastx_on_fasta_file(fasta_file, out_file):
	# record = SeqIO.read(fasta_file, format="fasta")
	fasta_string = open(fasta_file).read()
	# result_handle = NCBIWWW.qblast("blastx", "nr", fasta_string.format("fasta"), entrez_query="mouse[orgn]", format_type='Text', expect='0.00000001', hitlist_size=5)
	result_handle = NCBIWWW.qblast("blastx", "nr", fasta_string.format("fasta"), entrez_query="mouse[orgn]", expect='0.000001', hitlist_size=5)
	with open(out_file, "w") as out_handle:
		out_handle.write(result_handle.read())
	result_handle.close()

def blastx_on_multiple_fasta_files(fasta_files, out_suffix):
	for fasta_file in fasta_files:
		out_file = fasta_file.rsplit('.', 1)[0] + out_suffix
		# record = SeqIO.read(fasta_file, format="fasta")
		fasta_string = open(fasta_file).read()
		# result_handle = NCBIWWW.qblast("blastx", "nr", fasta_string.format("fasta"), entrez_query="mouse[orgn]", format_type='Text', expect='0.00000001', hitlist_size=5)
		result_handle = NCBIWWW.qblast("blastx", "nr", fasta_string.format("fasta"), entrez_query="mouse[orgn]", expect='0.000001', hitlist_size=5)
		with open(out_file, "w") as out_handle:
			out_handle.write(result_handle.read())
		result_handle.close()

def make_accession_to_genename_dict():
	acc_gene_dict = {}
	line_count = 0
	with open('mouse.gene2accession', "r") as acc_gene_file:
		for line in acc_gene_file:
			line_count += 1
			if line_count > 1:
				line = line.strip('\n').split(delim)
				accession = line[5].split('.')[0]
				# print accession
				gene = line[-1]
				if accession != '-':
					if accession not in acc_gene_dict:
						acc_gene_dict[accession] = [gene]
					else:
						acc_gene_dict[accession] = acc_gene_dict[accession] + [gene]
	#remove duplicates in values
	for i in acc_gene_dict:
		acc_gene_dict[i] = list(set(acc_gene_dict[i])) 
	return acc_gene_dict


def parse_blast_results(blast_xml_files):
	##get dict to convert accession id to gene names
	accession_dict = make_accession_to_genename_dict()
	##set up dict for results
	results_dict = {}
	for blast_xml_file in blast_xml_files:
		print(blast_xml_file)
		with open(blast_xml_file, 'r') as result_handle:
			blast_records = NCBIXML.parse(result_handle)
			for blast_record in blast_records:
				for alignment in blast_record.alignments:
					# print(alignment.description)
					for hsp in alignment.hsps:
						query_id = blast_record.query
						e_value = hsp.expect
						alignment_info = alignment.title
						##get id from 4th pos and rm the version
						alignment_id = alignment_info.split('|')[3].split('.')[0]
						species = alignment_info.split('[')[-1].replace(']', '')
						if alignment_id in accession_dict:
							gene_id = accession_dict[alignment_id][0]
						else:
							gene_id = 'na'
						if gene_id != 'na' and species == 'Mus musculus':
							# print(query_id, e_value, alignment_id, gene_id)
							if query_id in results_dict:
								if gene_id not in results_dict[query_id]:
									results_dict[query_id].append(gene_id)
							else:
								results_dict[query_id] = [gene_id]
	'''
	for acomy_id in results_dict:
		print(acomy_id, results_dict[acomy_id])
	print(len(results_dict))
	'''
	return results_dict

def parse_blast_results_only_one_gene(blast_xml_files):

	##get dict to convert accession id to gene names
	accession_dict = make_accession_to_genename_dict()
	##make dict with acomy_id: [[gene_ids],[evalues]]
	results_dict = {}
	for blast_xml_file in blast_xml_files:
		# print(blast_xml_file)
		with open(blast_xml_file, 'r') as result_handle:
			blast_records = NCBIXML.parse(result_handle)
			for blast_record in blast_records:
				for alignment in blast_record.alignments:
					# print(alignment.description)
					for hsp in alignment.hsps:
						query_id = blast_record.query
						e_value = hsp.expect
						alignment_info = alignment.title
						##get id from 4th pos and rm the version
						alignment_id = alignment_info.split('|')[3].split('.')[0]
						species = alignment_info.split('[')[-1].replace(']', '')
						if alignment_id in accession_dict:
							gene_id = accession_dict[alignment_id][0]
							# if gene_id == 'Cdh6':
							# 	print('Cdh6', query_id, alignment_id)
							# elif gene_id == 'Cdh11':
							# 	print('Cdh11', query_id, alignment_id)								
						else:
							gene_id = 'na'
						if gene_id != 'na' and species == 'Mus musculus':
							# print(query_id, e_value, alignment_id, gene_id)
							if query_id in results_dict:
								if gene_id not in results_dict[query_id][0]:
									results_dict[query_id][0].append(gene_id)
									results_dict[query_id][1].append(e_value)
							else:
								results_dict[query_id] = [[gene_id],[e_value]]
	##make new dict = gene_id:[best acomy id, evalue] 
	gene_acomy_dict = {}
	for ai in results_dict:
		gene = results_dict[ai][0][0]
		evalue = float(results_dict[ai][1][0])
		if gene in gene_acomy_dict:
			if evalue < gene_acomy_dict[gene][1]:
				gene_acomy_dict[gene] = [ai, evalue]
		else:
			gene_acomy_dict[gene] = [ai, evalue]
	##then make back into acomy:[gene_id]
	acomy_gene_dict = {}
	for g in gene_acomy_dict:
		a_id = gene_acomy_dict[g][0]
		if a_id in acomy_gene_dict:
			print(a_id, 'already in dict')
		else:
			acomy_gene_dict[a_id] = [g]


	'''
	for acomy_id in acomy_gene_dict:
		print(acomy_id, acomy_gene_dict[acomy_id])
	print(len(acomy_gene_dict))
	'''
	return acomy_gene_dict


def add_gene_id_de_results(infiles, id_dict, out_suffix):
	for infile in infiles:
		outfile = infile.rsplit('.', 1)[0] + out_suffix
		lc = 0
		with open(infile, 'r') as in_fh, open(outfile, 'w') as out_fh:
			for line in in_fh:
				lc += 1
				if lc == 1:
					out_fh.write(line)
				else:
					line = line.split(',')
					gene_id = line[0].replace('"', '')
					##just print if we have a genename
					if gene_id in id_dict:
						new_gid = id_dict[gene_id][0]
						line_out = [new_gid] + line[1:]
						out_fh.write(','.join(line_out))
					else:
						new_gid = gene_id
					
					# print(gene_id, new_gid)
					
def keep_acomy_genes_in_mouse_data(infiles, id_dict, out_suffix):
	gene_dict = {}
	for id_name in id_dict:
		gene_id = id_dict[id_name][0]
		gene_dict[gene_id] = id_name
	for infile in infiles:
		outfile = infile.rsplit('.', 1)[0] + out_suffix
		lc = 0
		with open(infile, 'r') as in_fh, open(outfile, 'w') as out_fh:
			for line in in_fh:
				lc += 1
				if lc == 1:
					out_fh.write(line)
				else:
					line = line.split(',')
					gene_id = line[0].replace('"', '')
					##just print if we have a genename
					if gene_id in gene_dict:
						new_gid = gene_id
						line_out = [new_gid] + line[1:]
						out_fh.write(','.join(line_out))
						
					
					# print(gene_id, new_gid)			
def convert_gene_id_to_acomy_id(gene_acomy_dict, infile, outfile):
	##flip the dict
	acomy_id_dict ={}
	for acomy_id in gene_acomy_dict:
		gene_id = gene_acomy_dict[acomy_id][0]
		if gene_id in acomy_id_dict:
			acomy_id_dict[gene_id].append(acomy_id)
		else:
			acomy_id_dict[gene_id] = [acomy_id]
	with open(infile, 'r') as in_fh, open(outfile, 'w') as out_fh:
		for line in in_fh:
			gene = line.rstrip()
			if gene in acomy_id_dict:
				for acomy_name in acomy_id_dict[gene]:
					print(gene, acomy_name)
					out_fh.write(acomy_name + '\n')

def get_acomy_ids_from_blast_results(blast_xml_files, genes_to_query):
	##get dict to convert accession id to gene names
	accession_dict = make_accession_to_genename_dict()
	##make dict with acomy_id: [[gene_ids],[evalues]]
	results_dict = {}
	for blast_xml_file in blast_xml_files:
		# print(blast_xml_file)
		with open(blast_xml_file, 'r') as result_handle:
			blast_records = NCBIXML.parse(result_handle)
			for blast_record in blast_records:
				for alignment in blast_record.alignments:
					# print(alignment.description)
					for hsp in alignment.hsps:
						query_id = blast_record.query
						e_value = hsp.expect
						alignment_info = alignment.title
						##get id from 4th pos and rm the version
						alignment_id = alignment_info.split('|')[3].split('.')[0]
						species = alignment_info.split('[')[-1].replace(']', '')
						if alignment_id in accession_dict:
							gene_id = accession_dict[alignment_id][0]
							for gene_to_query in genes_to_query:
								if gene_id == gene_to_query:
									print(gene_id, query_id, alignment_id)


def print_gene_conversion(in_dict, outfile):
	with open(outfile, 'w') as out_fh:
		for acomy_id in in_dict:
			gene_ids = in_dict[acomy_id]
			print(acomy_id, gene_ids)
			out_fh.write(delim.join([acomy_id] + gene_ids) + '\n')


##run methods
gtf_file = '/data/atimms/references/igenomes/acomy1/genes.gtf'
fa_file = 'acomy_draft_genes_0618.fa'
fa_file_split_prefix = 'acomy_draft_genes_0618.split'
split_out_suffix = '.blast.xml'
genome_fa = '/data/atimms/references/igenomes/acomy1/genome.fa'
blast_results_file = 'acomy_draft_genes_0618.blast.out'
acomy_gene_id_file = 'acomy_genes_to_id_0618.xls'
xml_blast_results = glob.glob('*' + split_out_suffix)
print(xml_blast_results)
# xml_blast_results = ['acomy_draft_genes_0618.split1.blast.xml', 'acomy_draft_genes_0618.split2.blast.xml']
acomy_de_files = ['acomy_rnaseq_0618_acomy_reduced.sham_vs_day2.csv', 'acomy_rnaseq_0618_acomy_reduced.sham_vs_day5.csv', 
		'acomy_rnaseq_0618_acomy.sham_vs_day2.csv', 'acomy_rnaseq_0618_acomy.sham_vs_day5.csv', 
		'acomy_rnaseq_0618_acomy_treatment.sham_vs_treatment.csv']
mouse_de_files = ['acomy_rnaseq_0618_mouse.sham_vs_day2.csv', 'acomy_rnaseq_0618_mouse.sham_vs_day5.csv', 
		'acomy_rnaseq_0618_mouse_treatment.sham_vs_treatment.csv']
# de_files = ['acomy_rnaseq_0618_acomy_treatment.sham_vs_treatment.csv']
de_altered_suffix = '.genename.csv'
de_altered_suffix_unique = '.genename_unique.csv'
de_acomy_genes_suffix = '.acomy_genes.csv'
dave_genelist_dna_trans = 'DNA_templated_transcription_gene_list.txt'
dave_genelist_dna_trans_acomy_id = 'DNA_templated_transcription_gene_list.acomy_id.txt'
sam_genelist_myo_tf = 'sam_myo_tf.txt'
sam_genelist_myo_tf_trans_acomy_id = 'sam_myo_tf.acomy_id.txt'
dave_genelist_bilavent = 'dave_bivalent_genes_1118.txt'
dave_genelist_bilavent_acomy_id = 'dave_bivalent_genes_1118.acomy_id.txt'
##test: 4 genes from LachesisContig2926
# gtf_file = 'test.gtf'
# fa_file = 'test.fa'
# gtf_file = 'test2.gtf'
# fa_file = 'test2.fa'
##1 20th of whole set
# fa_file_t30 = 'test_30kl.fa'
# blast_file_t30 = 'test_30kl.out'
# fa_file = 'test.fa'
# blast_results_file = 'test.out'

##get fasta formatted file for all genes
# make_fasta_from_gtf(gtf_file, fa_file, genome_fa)

##split fasta into 50 records per fasta
# split_fa_files = split_fasta(fa_file, fa_file_split_prefix, 50)
# print(split_fa_files)
##blast fa file against mouse nr database
##test set
# blastx_on_fasta_file(fa_file_t30, blast_file_t30) 
##big fasta file ##didn't work, so split into smaller files and run
#blastx_on_fasta_file(fa_file, blast_results_file)
##on split files ##if to many quesries get error: CPU usage limit was exceeded, resulting in SIGXCPU (24)
# blastx_on_multiple_fasta_files(split_fa_files, split_out_suffix)
#parse and then format
gene_id_dict = parse_blast_results(xml_blast_results)
print(len(gene_id_dict))
##add mouse gene_ids to de results
# add_gene_id_de_results(acomy_de_files, gene_id_dict, de_altered_suffix)
# keep_acomy_genes_in_mouse_data(mouse_de_files, gene_id_dict, de_acomy_genes_suffix)
##for getting heatmap from acomy data need to convert genename into acomy gene id
# convert_gene_id_to_acomy_id(gene_id_dict, dave_genelist_dna_trans, dave_genelist_dna_trans_acomy_id)
# convert_gene_id_to_acomy_id(gene_id_dict, sam_genelist_myo_tf, sam_genelist_myo_tf_trans_acomy_id)
convert_gene_id_to_acomy_id(gene_id_dict, dave_genelist_bilavent, dave_genelist_bilavent_acomy_id)
##so dave wants just one instance of each gene - change the parsing and redo the conversion
# gene_id_dict = parse_blast_results_only_one_gene(xml_blast_results)
# print(len(gene_id_dict))
##add mouse gene_ids to de results
# add_gene_id_de_results(acomy_de_files, gene_id_dict, de_altered_suffix_unique)


##get acomy id from gene id from blast results, then grep fa file (acomy_draft_genes_0618.fa)
# get_acomy_ids_from_blast_results(xml_blast_results, ['Cdh11', 'Cdh6'])

##print table of all acomy gene to id
print_gene_conversion(gene_id_dict, acomy_gene_id_file)
