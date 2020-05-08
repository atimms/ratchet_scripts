#!/usr/bin/env python
import sys
import subprocess
import os
from Bio.Seq import Seq
from Bio import SeqIO


'''
must load local python 3.5 to run
module load local_python/3.6.5
'''
##parameters
delim = '\t'

##setup working directory where results will be
working_dir = '/data/atimms/acomy_protein_seq_1018'
os.chdir(working_dir)

##programs
getorf = '/home/atimms/programs/EMBOSS-6.6.0/emboss/getorf'



def add_id_to_fasta(in_fasta, out_fasta, gene_id_file):
	acomy_dict = {}
	with open(gene_id_file, 'r') as gi_fh:
		for line in gi_fh:
			line = line.rstrip().split(delim)
			acomy_id = line[0]
			gene = line[1]
			# print(acomy_id, gene)
			if acomy_id in acomy_dict:
				print(acomy_id, 'seen twice')
			else:
				acomy_dict[acomy_id] = gene
	with open(in_fasta, 'r') as in_fh, open(out_fasta, 'w') as out_fh:
		for line in in_fh:
			if line[0] == '>':
				ai = line.rstrip()[1:]
				print(ai)
				if ai in acomy_dict:
					line_out = '>' + ai + '(' + acomy_dict[ai] + ')\n'
					out_fh.write(line_out)
				else:
					out_fh.write(line)				
			else:
				out_fh.write(line)



def convert_dna_to_protein(in_fasta, out_fasta):
	# with open(out_fasta, 'w') as out_fh:
	out_records = []
	##import fata
	for record in SeqIO.parse(in_fasta, "fasta"):
		print(record.id)
		print(record.seq)
		##convert to protein
		pseq = record.seq.translate(to_stop=True)
		##change the record to a protein seq
		record.seq = pseq
		print(record.seq)
		out_records.append(record)
		# print(pseq)
		# out_fh.write('>' + record.id + '\n')
		# out_fh.write(str(pseq) + '\n')
	print(out_records)
	SeqIO.write(out_records, out_fasta, "fasta")



##run methods
acomy_fasta = 'acomy_draft_genes_0618.fa'
acomy_fasta_with_gene_ids = 'acomy_genes_mouse_ids_0618.fa'
acomy_protein_fasta = 'acomy_protein_1018.fa'
acomy_gene_ids = 'acomy_genes_to_id_0618.xls'
##test
# acomy_fasta = 'test.fa'
# acomy_fasta_with_gene_ids = 'test_id.fa'





##modify fatsa with gene id
add_id_to_fasta(acomy_fasta, acomy_fasta_with_gene_ids, acomy_gene_ids)


##convert to protein
convert_dna_to_protein(acomy_fasta_with_gene_ids, acomy_protein_fasta)