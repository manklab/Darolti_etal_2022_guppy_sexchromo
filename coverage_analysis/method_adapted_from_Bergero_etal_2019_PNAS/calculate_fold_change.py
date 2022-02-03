#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"This script reads in female and male coverage data and calculates and outputs M:F coverage fold change for each chromosome block"
#==============================================================================
import argparse
import sys
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("females", type=str,
                    help="A file with average coverage data from females")
parser.add_argument("males", type=str,
                    help="A file with average coverage data from males")
parser.add_argument("outfile", type=str,
                    help="An output file with M:F coverage fold change values")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def extract_depth(source):
    depthdict = defaultdict(list)
    with open(source, "r") as infile:
        for line in infile:
            line = line.rstrip()
            if not line.startswith("Chromosome"):
                line = line.split(",")
                chromosome = line[0]
                windowStart = line[1]
                windowEnd = line[2]
                chromosome_block = chromosome + "_" + windowStart + "_" + windowEnd
                sumdepth = line[4]
                averagedepth = line[5]
                logaveragedepth = line[6]
                depthdict[chromosome_block].append(sumdepth)
                depthdict[chromosome_block].append(averagedepth)
                depthdict[chromosome_block].append(logaveragedepth)
    return depthdict

def combine_depth(females_depth, males_depth):
    combined_dict = defaultdict(list)
    for block in females_depth: 
        fem_sumdepth = float(females_depth[block][0])
        fem_averagedepth = float(females_depth[block][1])
        fem_logaveragedepth = float(females_depth[block][2])
        
        mal_sumdepth = float(males_depth[block][0])
        mal_averagedepth = float(males_depth[block][1])
        mal_logaveragedepth = float(males_depth[block][2])

        #combined sum and averages are identical
        combined_logaveragedepth = mal_logaveragedepth - fem_logaveragedepth

        combined_dict[block].append(combined_logaveragedepth)
        combined_dict[block].append(mal_averagedepth)
        combined_dict[block].append(mal_logaveragedepth)
        combined_dict[block].append(fem_averagedepth)
        combined_dict[block].append(fem_logaveragedepth)
    
    return combined_dict
           
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    #extract coverage in males and females
    females_depth = extract_depth(args.females)
    print "Number of chromosome blocks with coverage in females =", len(females_depth)
    males_depth = extract_depth(args.males)
    print "Number of chromosome blocks with coverage = males", len(males_depth)

    #check identical chromosome_blocks
    for block in females_depth:
        if block not in males_depth:
            print "error"
    for block in males_depth:
        if block not in females_depth:
            print "error"

    #combine coverage
    combined_depth = combine_depth(females_depth, males_depth)
    print "Number of chromosome blocks with average coverage =", len(combined_depth)

    with open(args.outfile, "w") as outfile:
        header = "Block,Chromosome,WindowStart,WindowEnd,MFLogaverage,Maverage,Mlogaverage,Faverage,Flogaverage"
        outfile.write(header)
        outfile.write("\n")
        for block in combined_depth:
            chromosome = block.split("_")[0] + "_" + block.split("_")[1] 
            windowStart = block.split("_")[2]
            windowEnd = block.split("_")[3]
            outfile.write(block)
            outfile.write(",")
            outfile.write(chromosome)
            outfile.write(",")
            outfile.write(windowStart)
            outfile.write(",")
            outfile.write(windowEnd)
            outfile.write(",")
            outfile.write(str(combined_depth[block][0]))
            outfile.write(",")
            outfile.write(str(combined_depth[block][1]))
            outfile.write(",")
            outfile.write(str(combined_depth[block][2]))
            outfile.write(",")
            outfile.write(str(combined_depth[block][3]))
            outfile.write(",")
            outfile.write(str(combined_depth[block][4]))
            outfile.write("\n")

if __name__ == '__main__':
    main()