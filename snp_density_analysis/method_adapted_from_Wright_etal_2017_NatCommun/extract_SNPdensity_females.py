#!/usr/bin/python2.7
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
parser.add_argument("SNP", type=str,
                    help="A file of scaffolds and SNPs")
parser.add_argument("outfilepath", type=str,
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
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if not f.endswith(".DS_Store")]

def list_files(current_dir):
    file_list = []
    for path, subdirs, files in os.walk(current_dir): # Walk directory tree
        for name in files:
            if not name.endswith(".DS_Store"):
                f = os.path.join(path, name)
                file_list.append(f)
    return file_list

# def get_file_dict(infolders, popdict):
def get_file_dict(infolders):
    filedictionary = {}
    for infolder in infolders:
        name = os.path.basename(infolder)
        # print name
        infiles = list_files(infolder)
        for infile in infiles:
            # print infile
            if infile.endswith("_norm"):
                filedictionary[name] = infile
    print "No. of samples =", len(filedictionary)
    return filedictionary

def extract_depth(source, popsnpdict, genescaffolddict):
    count = 0
    with open(source, "r") as infile:
        next(infile)
        for line in infile:
            count += 1
            line = line.rstrip()
            line = line.split("\t")
            scaffold = line[1]
            gene = line[2]
            density = line[-1]
            sites = float(line[3])
            snps = float(line[4])
            row = [sites, snps, density]
            popsnpdict[gene].append(row)
            if gene in genescaffolddict:
                if scaffold != genescaffolddict[gene]:
                    print "error"
            else:
                genescaffolddict[gene] = scaffold
    return popsnpdict, genescaffolddict

def average_density(popsnpdict):
    averagedict = defaultdict(list)
    for gene in popsnpdict:
        if len(popsnpdict[gene]) == 3:
            count = 0
            sumsites = 0
            sumsnps = 0
            sumdensity = 0
            for sample in popsnpdict[gene]:
                sites = float(sample[0])
                snps = float(sample[1])
                density = float(sample[2])
                sumsites += sites
                sumsnps += snps
                sumdensity += density
                if sites > 0:
                    count += 1
            if count != 3:
                pass
            else:
                averagedensity = sumdensity/3
                logaveragedensity = np.log2(averagedensity+1)
                averagedict[gene].append(sumsites)
                averagedict[gene].append(sumsnps)
                averagedict[gene].append(averagedensity)
                averagedict[gene].append(logaveragedensity)
    return averagedict
           
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

    #get dictionary of files
    SNPfolders = list_folder(args.SNP)
    print "No. of infolders =", len(SNPfolders)
    filedictionary = get_file_dict(SNPfolders)

    #get SNP density for each sample
    snpdicty = defaultdict(list)
    genescaffolddict = {}

    for name in filedictionary:
        snpdict, genescaffolddict = extract_depth(filedictionary[name], snpdicty, genescaffolddict)
    print "Number of genes =", len(snpdict), len(genescaffolddict)

    # calculate average SNP density
    average = average_density(snpdict)

    outfilename = args.outfilepath+"/snpdensity_females.txt"
    with open(outfilename, "w") as outfile:
        header = "Scaffold,Gene,Sumsites,Sumsnps,Snpdensity,LogSnpdensity"
        outfile.write(header)
        outfile.write("\n")
        for gene in average:
            scaffold = genescaffolddict[gene]
            outfile.write(scaffold)
            outfile.write(",")
            outfile.write(gene)
            outfile.write(",")
            outfile.write(str(average[gene][0]))
            outfile.write(",")
            outfile.write(str(average[gene][1]))
            outfile.write(",")
            outfile.write(str(average[gene][2]))
            outfile.write(",")
            outfile.write(str(average[gene][3]))
            outfile.write("\n")

if __name__ == '__main__':
    main()