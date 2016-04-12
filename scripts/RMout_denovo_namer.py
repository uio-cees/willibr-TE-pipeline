#!/usr/bin/env python
"""
FINISHED, William Brynildsen
"""

from Bio import SeqIO
import sys, argparse

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=
	'Renames de novo repeats from .outfile with class/family and species information')
	parser.add_argument('-i', '--input', action='store', help='', type=argparse.FileType('r'), default = '-')
	parser.add_argument('-n', '--number')

	args = parser.parse_args()
	
fish_number = args.number
		
def element_namer(outfile, number):
	str = ""
	for line in outfile:
		linesplit = line.split()
		if "rnd" in linesplit[9]:
			str += "sed -i 's/" + linesplit[9] + "/" + "fish-%s" % number + "-" + linesplit[9] + "-" + linesplit[10] + "/g' $1" + "\n"
		elif "TRIM" in linesplit[9]:
			str += "sed -i 's/" + linesplit[9] + "/" + "fish-%s" % number + "-" + linesplit[9] + "-" + linesplit[10] + "/g' $1" + "\n"
		elif "retro" in linesplit[9]:
			str += "sed -i 's/" + linesplit[9] + "/" + "fish-%s" % number + "-" + linesplit[9] + "-" + linesplit[10] + "/g' $1" + "\n"
	return str.rstrip()

sed_list_of_named_elements = element_namer(args.input, args.number)

print sed_list_of_named_elements