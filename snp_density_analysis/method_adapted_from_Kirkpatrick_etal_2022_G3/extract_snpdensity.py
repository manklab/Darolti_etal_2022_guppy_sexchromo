#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"This script reads snp density information from multiple samples and calculates and outputs average snp density for each chromosome block"
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
import numpy as np
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infolder", type=str, nargs="+",
                    help="A folder with individual snp density files")
parser.add_argument("outfile", type=str,
                    help="")
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
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if f.endswith(".txt")]

def extract_snpdensity(source, snpdict):
    with open(source, "r") as infile:
        for line in infile:
            line = line.rstrip()
            line = line.split()
            if len(line) == 6:
                snpdensity = float(line[5])
                scaffold = line[1] + "_" + line[2] + "_" + line[3]
                snpdict[scaffold].append(snpdensity)
    return snpdict

# creates dictionary of sequence (scaff/contig) and its snp density (just as in def extract_depth) but it also calculates average snp density for each sequence
def average_snpdensity(snpdensity):
    averagedict = defaultdict(list)
    for sequence in snpdensity:
        snplist = []
        for value in snpdensity[sequence]:
            snplist.append(value)
        sum_snplist = sum(snplist)
        average = (sum_snplist/len(snplist))
        snp_values = snplist[0]
        for value in snplist[1:]:
            snp_values = str(snp_values)+"|"+str(value)
        logaverage = np.log2(average+1)
        sumsnp = sum_snplist
        averagedict[sequence].append(snp_values)
        averagedict[sequence].append(sumsnp)
        averagedict[sequence].append(average)
        averagedict[sequence].append(logaverage)
    return averagedict
           
# #==============================================================================
# #Main==========================================================================
# #==============================================================================
def main():

    #get files
    if len(args.infolder) == 1:
        if os.path.isdir(args.infolder[0]):
            infiles = list_folder(args.infolder[0])
        else:
            infiles = args.infolder
    else:
        infiles = args.infolder

    #extract snp density
    snpdensitydict = defaultdict(list)
    for infile in infiles:
        print infile 
        snpdensity = extract_snpdensity(infile, snpdensitydict)
        print "Number of chromosome blocks with snp density =", len(snpdensity)

    #calculate average snp density
    average = average_snpdensity(snpdensity)
    print "Number of chromosome blocks with average snp density =", len(average)

    with open(args.outfile, "w") as outfile:
        print args.outfile
        header = "Chromosome,WindowStart,WindowEnd,SnpDensityList,Sum,Average,Logaverage"
        outfile.write(header)
        outfile.write("\n")
        for sequence in average:
            chromosome = sequence.split("_")[0]
            windowStart = sequence.split("_")[1]
            windowEnd = sequence.split("_")[2]
            outfile.write(chromosome)
            outfile.write(",")
            outfile.write(str(windowStart))
            outfile.write(",")
            outfile.write(str(windowEnd))
            outfile.write(",")
            outfile.write(str(average[sequence][0]))
            outfile.write(",")
            outfile.write(str(average[sequence][1]))
            outfile.write(",")
            outfile.write(str(average[sequence][2]))
            outfile.write(",")
            outfile.write(str(average[sequence][3]))
            outfile.write("\n")

if __name__ == '__main__':
    main()