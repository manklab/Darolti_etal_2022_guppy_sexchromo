#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"This script ....."
#==============================================================================
import argparse
import sys
import os
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infolder", type=str,
					help="A folder containing a separate folder for each sample with coverage data for each chromosome")
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
    '''Returns a list of all sample folders'''
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if os.path.isdir(f)]

def list_files(current_dir):
    file_list = []
    for path, subdirs, files in os.walk(current_dir): # Walk directory tree
        for name in files:
            if not name.endswith(".DS_Store"):
                f = os.path.join(path, name)
                file_list.append(f)
    return file_list
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

	infolders = list_folder(args.infolder)
	for sample_folder in infolders:
		outfile = sample_folder + "_50kb_allchromosomes.txt"
		infiles = [f for f in list_files(sample_folder)]
		with open(outfile, "w") as out:
			for file in infiles:
				if file.endswith("_50kb.txt"):
					with open(file, "r") as file:
						for line in file:
							out.write(line)

if __name__ == '__main__':
	main()