#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"This script reads in female and male snp density data and calculates and outputs M:F snp density fold change for each chromosome block"
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("females", type=str,
                    help="A file with average snp density data from females")
parser.add_argument("males", type=str,
                    help="A file with average snp density data from males")
parser.add_argument("outfile", type=str,
                    help="Output file containing fold change values")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def extract_snps(source):
    snpsdict = defaultdict(list)
    with open(source, "r") as infile:
        next(infile)
        for line in infile:
            line = line.rstrip().split(",")
            chromosome = line[0]
            windowStart = line[1]
            windowEnd = line[2]
            chromosome_block = chromosome + "_" + windowStart + "_" + windowEnd
            sum_snpdensity = line[4]
            average_snpdensity = line[5]
            logaverage_snpdensity = line[6]
            snpsdict[chromosome_block].append(sum_snpdensity)
            snpsdict[chromosome_block].append(average_snpdensity)
            snpsdict[chromosome_block].append(logaverage_snpdensity)
    return snpsdict

def combine_density(females_snps, males_snps):
    combined_dict = defaultdict(list)
    for block in females_snps: 
        fem_sum = float(females_snps[block][0])
        fem_average = float(females_snps[block][1])
        fem_logaverage = float(females_snps[block][2])
        mal_sum = float(males_snps[block][0])
        mal_average = float(males_snps[block][1])
        mal_logaverage = float(males_snps[block][2])
        combined_logaverage = mal_logaverage - fem_logaverage
        combined_dict[block].append(combined_logaverage)
        combined_dict[block].append(mal_average)
        combined_dict[block].append(mal_logaverage)
        combined_dict[block].append(fem_average)
        combined_dict[block].append(fem_logaveraged)
    return combined_dict
           
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    #extract snps in males and females
    females_snps = extract_snps(args.females)
    print "Number of feamle chromosome blocks =", len(females_depth)
    males_snps = extract_snps(args.males)
    print "Number of male chromosome blocks = ", len(males_depth)

    #check identical chromosome_blocks
    for block in females_depth:
        if block not in males_depth:
            print "error"
    for block in males_depth:
        if block not in females_depth:
            print "error"

    #combine coverage
    combined_density = combine_density(females_snps, males_snps)
    print "Number of filtered chromosome blocks =", len(combined_density)

    with open(args.outfile, "w") as outfile:
        header = "ChromosomeBlock,Chromosome,WindowStart,WindowEnd,MFLogaverage,Maverage,Mlogaverage,Faverage,Flogaverage"
        outfile.write(header)
        outfile.write("\n")
        for block in combined_depth:
            chromosome = block.split("_")[0]
            windowStart = block.split("_")[1]
            windowEnd = block.split("_")[2]
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