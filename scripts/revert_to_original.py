#!/usr/bin/env python
"""
Purpose: To rewrite fasta files, replacing alias headers with original headers,
using the .table file created during pipeline.

William Brynildsen
"""

from Bio import SeqIO
import sys, argparse

if __name__ == '__main__':
	
	
	parser = argparse.ArgumentParser(description=
	'')
	parser.add_argument('-i', '--input', action='store', help='', type=argparse.FileType('r'), default = '-')
	parser.add_argument('-t', '--table', action='store', help='', type=argparse.FileType('r'), default = '-')
	args = parser.parse_args()
	
#	Loads in table generated in repeat pipeline.

	table = args.table.readlines()

#	Function that accepts the new header (like 'SEQ0') and outputs the original header

def return_original_header(newheader):
	newheader = ">" + newheader
	for line in table:
		linesplit = line.split()
		if newheader == linesplit[1]:
			print linesplit[0]

#	Writes the complete fasta file with new headers. Pipe into a file. Not formatted.

for record in SeqIO.parse(args.input, 'fasta'):
 	return_original_header(record.id)
	print record.seq
	
	