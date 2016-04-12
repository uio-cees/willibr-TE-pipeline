#!/usr/bin/env python
"""
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
	
#	Loads in table generated in repeat pipeline and a GFF file.

	table = args.table.readlines()
	gff = args.input.readlines()
	
#	Function that compares input (gff line) with name in table and prints a new line with the linked original name.

def print_line(lines):
	for line in table:
		linesplit = line.split()
		gff_compare = ">" + lines[0]
		if gff_compare == linesplit[1]:
			print linesplit[0].split(">")[1] + "\t" + lines[1] + "\t" + lines[2] +  "\t" + lines[3] + "\t" + lines[4] + "\t" + lines[5] + "\t" + lines[6] + "\t" + lines[7] + "\t" + lines[8]

#	Runs function for each line
for i in gff:
	print_line(i.split())
	

	
	
