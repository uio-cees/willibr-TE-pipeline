#!/bin/sh

from Bio import SeqIO
import sys, argparse

if __name__ == '__main__':
	
	
	parser = argparse.ArgumentParser(description=
	'Input file')
	parser.add_argument('-i', '--input', action='store', help='', type=argparse.FileType('r'), default = '-')	

	args = parser.parse_args()
print "Alias","Repeat","Count","Bp_masked","Percentage_masked"

tmp = ""

for line in args.input:
	splitline = line.split()
	if "Class" not in line:	
		if "Count" not in line:
			if "*" not in line:
				if len(line.split(":")[-1].lstrip().rstrip()) > 3:
					print line.split(".")[0], line.split(":")[-1].lstrip().rstrip()
			
			
				 
