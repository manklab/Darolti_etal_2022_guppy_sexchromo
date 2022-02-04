#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"This script calculates coverage in 50kb windows"
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
		infiles = [f for f in list_files(sample_folder)]
		for file in infiles:
			name = file.split("/")[-2] + "_" + file.split("/")[-1].split(".txt")[0]
			outfile = os.path.dirname(file)+"/"+name+"_50kb.txt"
			count = 0
			target = 50000 # For P. wingei analyses that use the P. wingei de novo assembly we used a window size of 10000 instead
			coverage_sum = 0 
			coverage_average = 0
			with open(outfile, "w") as out:
				with open(file, "r") as file:
					for line in file:
						line = line.rstrip()
						count += 1
						pos = int(line.split()[1])
						coverage_value = int(line.split()[2])
						coverage_sum += coverage_value
						if count == target:
							coverage_average = coverage_sum / 50000
							target += 50000
							out.write(line)
							out.write("\t")
							out.write(str(coverage_sum))
							out.write("\t")
							out.write(str(coverage_average))
							out.write("\n")
							coverage_sum = 0
							coverage_average = 0

if __name__ == '__main__':
	main()