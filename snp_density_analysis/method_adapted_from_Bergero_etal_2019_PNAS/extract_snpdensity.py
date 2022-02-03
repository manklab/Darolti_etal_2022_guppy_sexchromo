#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
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
parser.add_argument("coverage", type=str, nargs="+",
                    help="")
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

def extract_depth(source, depthdict):
    #name, bpcovered/length, coverage, averagedepth (number of reads covering a base/coverage)
    with open(source, "r") as infile:
        for line in infile:
            line = line.rstrip()
            line = line.split()
            if len(line) == 6:
                depth = line[5]
                # depth = float(depth)
                depth = float(depth)
                scaffold = line[1] + "_" + line[2] + "_" + line[3]
                depthdict[scaffold].append(depth)
    return depthdict

# creates dictionary of sequence (scaff/contig) and its depth (just as in def extract_depth) but it also calculates 
# average depth for each sequence
def average_depth(coverage_depth):
    averagedict = defaultdict(list)
    for sequence in coverage_depth:
        coverage_list = []
        for depth in coverage_depth[sequence]:
            coverage_list.append(depth)
        sum_coverage_list = sum(coverage_list)
        average = (sum_coverage_list/len(coverage_list))

        # coverage_list[0] is the coverage of the first individual 
        #  (i think it assigns coverage_list[0] first in order to have a starting point to keep adding coverages)
        coverage_values = coverage_list[0]
        for cov in coverage_list[1:]:
            coverage_values = str(coverage_values)+"|"+str(cov)

        logaverage = np.log2(average+1)
        sumdepth = sum_coverage_list
        averagedict[sequence].append(coverage_values)
        averagedict[sequence].append(sumdepth)
        averagedict[sequence].append(average)
        averagedict[sequence].append(logaverage)
    return averagedict
           
# #==============================================================================
# #Main==========================================================================
# #==============================================================================
def main():
    #get files
    if len(args.coverage) == 1:
        if os.path.isdir(args.coverage[0]):
            infiles = list_folder(args.coverage[0])
        else:
            infiles = args.coverage
    else:
        infiles = args.coverage

    #extract coverage
    depthdict = defaultdict(list)
    for infile in infiles:
        print infile 
        coverage_depth = extract_depth(infile, depthdict)
        print "Number of chromosome blocks with coverage =", len(coverage_depth)

    #take average
    average = average_depth(coverage_depth)
    print "Number of chromosome blocks with average coverage =", len(average)

    with open(args.outfile, "w") as outfile:
        print args.outfile
        header = "Chromosome,WindowStart,WindowEnd,CovList,Sum,Average,Logaverage"
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