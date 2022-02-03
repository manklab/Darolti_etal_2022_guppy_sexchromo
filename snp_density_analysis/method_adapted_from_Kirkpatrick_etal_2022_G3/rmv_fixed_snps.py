#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
#==============================================================================
import argparse
import sys
import os
from subprocess import Popen, list2cmdline
from itertools import product
from collections import defaultdict

#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infolder", type=str,
                    help="A folder containing sample vcf files")
parser.add_argument("outfolder", type=str,
                    help="A folder in which to write the output vcf files")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
'''Returns a list of all files in a folder with the full path'''
def list_folder(infolder):
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if f.endswith("vcffilt1.vcf")]

'''Read vcf file and append site and genotype information to dictionary'''
def read_file(infile, dicty):
	with open(infile, "r") as infile:
		for line in infile:
			if not line.startswith("#"):
				chromosome = line.split()[0]
				position = line.split()[1]
				chromo_pos = chromosome + "_" + position
				genotype = line.split()[-1].split(":")[0]
				dicty[chromo_pos].append(genotype)
	return dicty

def write_output(infile, sites_to_include):
	outfile = os.path.basename(infile).split(".")[0] + "_filt2.vcf"
	with open(outfile, "w") as out:
		with open(infile, "r") as infile:
			for line in infile:
				if line.startswith("#"):
					out.write(line)
				else:
					if not line.startswith("P") and not line.startswith("M"):
						chromosome = line.split()[0]
						position = line.split()[1]
						if position in sites_to_include[chromosome]:
							out.write(line)
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
	dicty = defaultdict(list)
	infiles = list_folder(args.infolder)
	for infile in infiles:
		print infile
		dicty = read_file(infile, dicty)
	
	print "Number of unfiltered sites = ", len(dicty)

	# Make dictionary of chromosomes and sites to include in the final vcf files
	# Sites are included if there are at least two heterozygous samples (thus excluding fixed sites and singletons)
	sites_to_include = defaultdict(list)
	count_filtered_sites = 0
	for site in dicty:
		if len(dicty[site]) > 1:
			count_heterozygous = 0
			for genotype in dicty[site]:
				if genotype == "0/1":
					count_heterozygous += 1
			if count_heterozygous > 1:
				count_filtered_sites += 1
				if not site.startswith("P") and not site.startswith("M"):
					chromosome = site.split("_")[0]
					position = site.split("_")[-1]
					sites_to_include[chromosome].append(position)

	print "Number of filtered sites = ", count_filtered_sites

	os.chdir(args.outfolder)
	for infile in infiles:
		write_output(infile, sites_to_include)

if __name__ == "__main__":
    main()