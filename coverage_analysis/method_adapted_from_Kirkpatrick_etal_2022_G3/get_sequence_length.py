#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"This script reads in a fasta file and outputs information about the length of each sequence"
#==============================================================================
import argparse
import sys
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("fasta", type=str,
					help="A fasta file")
parser.add_argument("outfile", type=str,
					help="Output file of the format: SequenceName \t Length \n")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def read_sequences(fasta):
	SG={}
	try:
		with open(fasta, "r") as fasta:
			for line in fasta.readlines():
				if line[0] == ">":
					name = line[1:].rstrip().split()[0]
					header = line[1:].rstrip()
					SG[name] = ["", header]
				else:
					SG[name][0] += line.rstrip()
			return SG
	except IOError:
		print "!----ERROR----!"
		print "File %s does not exit!" % fasta
		sys.exit(1)
	except KeyboardInterrupt:
		sys.exit(1)

def get_all_length(seq_dict):
	dict = {}
	for seq in seq_dict:
		dict[seq] = int(len(seq_dict[seq][0]))
	return dict
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

	seq_dict = read_sequences(args.fasta)
	print "Number of sequences = ", len(seq_dict)

	all_length_dict = get_all_length(seq_dict)
	
	with open(args.outfile, "w") as out:
		for seq in all_length_dict:
			out.write(seq)
			out.write("\t")
			out.write(str(all_length_dict[seq]))
			out.write("\n")

if __name__ == '__main__':
	main()