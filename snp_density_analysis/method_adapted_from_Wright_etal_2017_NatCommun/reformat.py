#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"This script takes alignment files and modifies them to match the formatting requirements for Bow2pro"
#==============================================================================
import argparse
import sys
import os
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infolder", type=str,
					help="A folder with snpanalysis_bowtie2_sort.sam files")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def list_folder(infolder):
    '''Returns a list of all files in a folder with the full path'''
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if f.endswith("snpanalysis_bowtie2_sort.sam")]
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
	
	infiles = list_folder(args.infolder)
	for infile in infiles:
		outfile = infile.split(".")[0]+"_reformat.sam"
		with open(outfile, "w") as out:
			with open(infile, "r") as infile:
				for line in infile:
					if not line.startswith("@"):
						chromosome = line.split()[2]
						if chromosome != "*":
							line = line.split("\t")
							for element in line[0:4]:
								out.write(element)
								out.write("\t")
							out.write(line[9])
							out.write("\t")
							out.write(line[10])
							out.write("\n")

if __name__ == '__main__':
	main()