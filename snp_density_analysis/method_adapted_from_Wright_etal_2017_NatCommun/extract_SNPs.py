#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"Script for extracting SNPs within coding regions"
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("names", type=str,
                    help="A folder with a file containing sequence names")
parser.add_argument("cov", type=str,
                    help="A folder with a file containing SNPs filtered for minimum site coverage 10")
parser.add_argument("maf", type=str,
                    help="A folder with a file containing SNPs filtered for minimum site coverage 10 and SNP frequency > 0.3 x site coverage")
parser.add_argument("coordinates", type=str,
                    help="A file with gene coordinates in chromosomes")
parser.add_argument("outfilepath", type=str,
                    help="A folder for the output file")
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

def get_file_dict(infolders, filedictionary, uniqid):
    for infolder in infolders:
        name = os.path.basename(infolder)
        infiles = list_files(infolder)
        for infile in infiles:
            if infile.endswith(uniqid):
                filedictionary[name].append(infile)
    print "No. of samples =", len(filedictionary)
    return filedictionary

def get_coordinates(source):
    coord_dict = defaultdict(list)
    gene_dict = defaultdict(list)
    with open(source, "r") as coordinates:
        next(coordinates)
        for line in coordinates:
            line = line.rstrip()
            chromo = line.split()[0]
            info = line.split()[1:]
            for i in info:
                gene = i.split(",")[0]
                first = i.split(",")[1]
                second = i.split(",")[2]
                if first < second:
                    start = first
                    end = second
                else:
                    start = second
                    end = first
                row = [start,end]
                gene_dict[gene].append(row)
                coord_dict[chromo].append(gene)
    print "No. of chromosomes (from coordinates) =", len(coord_dict)
    print "No. of genes (from coordinates) =", len(gene_dict)
    return coord_dict, gene_dict

def get_names_dict(source):
    names_dict = {}
    with open(source, "r") as names:
        counter = 0
        for line in names:
            counter += 1
            line = line.rstrip()
            names_dict[counter] = line
    print "No. of chromosomes (from names) =", len(names_dict)
    return names_dict

def get_coverageormaf_dict(source, names_dict, gene_dict, coord_dict):
    coverageormaf_dict = {}
    for gene in gene_dict:
        coverageormaf_dict[gene] = 0
    with open(source, "r") as coverageormaf:
        counter = 0
        for line in coverageormaf:
            line = line.rstrip()
            if line.startswith(">"):
                counter += 1
                sys.stdout.write('%d\r' % (counter))
                sys.stdout.flush()
                scaffoldname = names_dict[counter]
                if scaffoldname in coord_dict:
                    genes = list(set(coord_dict[scaffoldname]))
                else:
                    genes = None
            else:
                if genes == None:
                    pass
                else:
                    pos = float(line.split("\t")[0])
                    for gene in genes:
                        g = gene_dict[gene]
                        include = None
                        for i in g:
                            start = float(i[0])
                            end = float(i[1])
                            if pos >= start and pos <= end:
                                include = "Yes"
                        if include == "Yes":
                            coverageormaf_dict[gene] +=1
    print "No. of chromosomes (from coverage or maf) =", len(coverageormaf_dict)
    return coverageormaf_dict
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    #order in names corresponds to same order in maf and cov
    #bow2pro renames therefore names contains proper contig names

    #get coordinates
    coord_dict, gene_dict = get_coordinates(args.coordinates)

    #get dictionary of files
    filedictionary = defaultdict(list)
    namefolder = list_folder(args.names)
    print "No. of infiles =", len(namefolder)
    filedictionary = get_file_dict(namefolder, filedictionary, ".list")

    covfolder = list_folder(args.cov)
    print "No. of infiles =", len(covfolder)
    filedictionary = get_file_dict(covfolder, filedictionary, "_sites10.pro")

    maffolder = list_folder(args.maf)
    print "No. of infiles =", len(maffolder)
    filedictionary = get_file_dict(maffolder, filedictionary, "_SNP30.pro")

    print filedictionary
    for sample in filedictionary:
        print sample
        files = filedictionary[sample]
        for f in files:
            if f.endswith(".list"):
                name = f
            elif f.endswith("_sites10.pro"):
                cov = f
            elif f.endswith("_SNP30.pro"):
                maf = f

        #make dictionary between bow2pro names and scaffold names
        names_dict = get_names_dict(name)
        
        #get number of sites
        coverage_dict = get_coverageormaf_dict(cov, names_dict, gene_dict, coord_dict)
        
        #get number of snps
        maf_dict = get_coverageormaf_dict(maf, names_dict, gene_dict, coord_dict)
        
        outfilename = args.outfilepath+"/"+sample+"/"+sample+"map_sorted_sites10_SNP30.bygene"
        print "Printing to ...", outfilename
        with open(outfilename,"w") as outfile:
            sorted_names = sorted(names_dict.keys())
            print len(sorted_names)
            for order in sorted_names:
                scaffold = names_dict[order]
                genes = list(set(coord_dict[scaffold]))
                for g in genes:
                    outfile.write(scaffold)
                    outfile.write("\t")
                    outfile.write(g)
                    outfile.write("\t")
                    outfile.write(str(coverage_dict[g]))
                    outfile.write("\t")
                    outfile.write(str(maf_dict[g]))
                    outfile.write("\n")

if __name__ == '__main__':
    main()