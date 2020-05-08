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


def parse_blast_results(blast_xml_files):
	for blast_xml_file in blast_xml_files:
		print(blast_xml_file)
		with open(blast_xml_file, 'r') as result_handle:
			blast_records = NCBIXML.parse(result_handle)
			# blast_record = next(blast_records)
			# for alignment in blast_record.alignments:
			# 	for hsp in alignment.hsps:
			# 		print(hsp.expect)
			for blast_record in blast_records:
				# print(blast_record)
				for alignment in blast_record.alignments:
					# print(alignment)
					for hsp in alignment.hsps:
						print(hsp)
						print("sequence:", alignment.title)
						print("e value:", hsp.expect)

##run methods
de_files = ['acomy_rnaseq_0618_acomy.sham_vs_day2.csv']
gtf_file = '/data/atimms/references/igenomes/acomy1/genes.gtf'
fa_file = 'acomy_draft_genes_0618.fa'
fa_file_split_prefix = 'acomy_draft_genes_0618.split'
split_out_suffix = '.blast.xml'
genome_fa = '/data/atimms/references/igenomes/acomy1/genome.fa'
blast_results_file = 'acomy_draft_genes_0618.blast.out'
xml_blast_results = ['acomy_draft_genes_0618.split10.blast.xml', 'acomy_draft_genes_0618.split1.blast.xml', 
		'acomy_draft_genes_0618.split2.blast.xml', 'acomy_draft_genes_0618.split3.blast.xml', 
		'acomy_draft_genes_0618.split4.blast.xml', 'acomy_draft_genes_0618.split5.blast.xml', 
		'acomy_draft_genes_0618.split6.blast.xml', 'acomy_draft_genes_0618.split7.blast.xml', 
		'acomy_draft_genes_0618.split8.blast.xml', 'acomy_draft_genes_0618.split9.blast.xml']
xml_blast_results = ['acomy_draft_genes_0618.split1.blast.xml']
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

##split fasta into 500 records per fasta
split_fa_files = split_fasta(fa_file, fa_file_split_prefix, 50)
print(split_fa_files)
##blast fa file against mouse nr database
##test set
# blastx_on_fasta_file(fa_file_t30, blast_file_t30) 
##big fasta file ##didn't work, so split into smaller files and run
#blastx_on_fasta_file(fa_file, blast_results_file)
##on split files
blastx_on_multiple_fasta_files(split_fa_files, split_out_suffix)
##parse and then format
# parse_blast_results(xml_blast_results)

