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
parser.add_argument("females", type=str,
                    help="A file of scaffolds and SNPs")
parser.add_argument("males", type=str,
                    help="A folder of scaffolds and SNPs")
parser.add_argument("outfilepath", type=str,
                    help="")
parser.add_argument("gene_pos", type=str,
                    help="File with genes and their position in the genome")
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

def get_file_dict(infolders, filedictionary):
    for infolder in infolders:
        name = os.path.basename(infolder)
        infiles = list_files(infolder)
        for infile in infiles:
            if infile.endswith("bygene.txt"):
                filedictionary[name] = infile
    print "No. of samples =", len(filedictionary)
    return filedictionary

def extract_snps(source):
    snpdict = defaultdict(list)
    genescaffolddict = {}
    count = 0
    with open(source, "r") as infile:
        for line in infile:
            count += 1
            line = line.rstrip()
            if not line.startswith("Scaffold"):
                line = line.split(",")
                scaffold = line[0]
                gene = line[1]
                sumsites = line[2]
                sumsnps = line[3]
                averagedensity = line[4]
                logaveragedensity = line[5]
                snpdict[gene].append(sumsites)
                snpdict[gene].append(sumsnps)
                snpdict[gene].append(averagedensity)
                snpdict[gene].append(logaveragedensity)
                genescaffolddict[gene] = scaffold
    return snpdict, genescaffolddict

def combine_density(females_snps, males_snps):
    combined_dict = defaultdict(list)
    #filter for lowly covered scaffolds ie in both male and female dictionary
    for gene in females_snps:
        if gene in males_snps:

            fem_averagedensity = float(females_snps[gene][2])
            fem_logaveragedensity = float(females_snps[gene][3])
            mal_averagedensity = float(males_snps[gene][2])
            mal_logaveragedensity = float(males_snps[gene][3])
            
            combined_logaveragedensity = mal_logaveragedensity - fem_logaveragedensity
            
            combined_dict[gene].append(combined_logaveragedensity)
            combined_dict[gene].append(mal_averagedensity)
            combined_dict[gene].append(mal_logaveragedensity)
            combined_dict[gene].append(fem_averagedensity)
            combined_dict[gene].append(fem_logaveragedensity)
    
    return combined_dict
           
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

    gene_position = {}
    gene_chromosome = {}

    with open(args.gene_pos, "r") as gene_pos:
        next(gene_pos)
        for line in gene_pos:
            gene = line.split(",")[1]
            position = line.split(",")[2]
            chromosome = line.split(",")[0]
            gene_position[gene] = position
            gene_chromosome[gene] = chromosome


    #extract snps in males and females
    females_snps, genescaffolddict  = extract_snps(args.females)
    print "Number of female scaffolds =", len(females_snps)

    males_snps, genescaffolddict = extract_snps(args.males)
    print "Number of male scaffolds =", len(males_snps)

    #combine snps
    combined_density = combine_density(females_snps, males_snps)
    print "Number of filtered scaffolds =", len(combined_density)

    outfilename = args.outfilepath+"/snpdensity_fc.txt"
    with open(outfilename, "w") as outfile:
        header = "Scaffold,Gene,MFLogaverage,Maverage,Mlogaverage,Faverage,Flogaverage,LG,Start"
        outfile.write(header)
        outfile.write("\n")
        for gene in combined_density:
            scaffold = genescaffolddict[gene]
            outfile.write(scaffold)
            outfile.write(",")
            outfile.write(gene)
            outfile.write(",")
            outfile.write(str(combined_density[gene][0]))
            outfile.write(",")
            outfile.write(str(combined_density[gene][1]))
            outfile.write(",")
            outfile.write(str(combined_density[gene][2]))
            outfile.write(",")
            outfile.write(str(combined_density[gene][3]))
            outfile.write(",")
            outfile.write(str(combined_density[gene][4]))
            outfile.write(",")
            outfile.write(gene_chromosome[gene])
            outfile.write(",")
            outfile.write(gene_position[gene])
            outfile.write("\n")

if __name__ == '__main__':
    main()