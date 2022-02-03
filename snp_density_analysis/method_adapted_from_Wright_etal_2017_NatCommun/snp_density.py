#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"This script calculates snp density as number of snps / number of sites for each gene and each sample"
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
parser.add_argument("snps", type=str,
                    help="A file of number of SNPs per gene")
parser.add_argument("coordinates", type=str,
                    help="A file of gene coordinates on chromosomes")
parser.add_argument("outfile", type=str,
                    help="An output file ")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()      
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    
    dicty = defaultdict(list)
    with open(args.coordinates, "r") as coordinates:
        next(coordinates)
        for line in coordinates:
            chromosome = line.split(",")[0]
            gene = line.split(",")[1]
            start = line.split(",")[2]
            dicty[gene].append(chromosome)
            dicty[gene].append(start)

    with open(args.outfile, "w") as out:
        with open(args.snps, "r") as snps:
            for line in snps:
                line = line.rstrip()
                chromosome = line.split()[0]
                gene = line.split()[1]
                sites = float(line.split()[2])
                snps = float(line.split()[3])
                if sites != 0:
                    density = snps/sites
                    logdensity = np.log2(density+1)
                    if gene in dicty:
                        if dicty[gene][0] == chromosome:
                            position = dicty[gene][1]
                            out.write(line)
                            out.write("\t")
                            out.write(str(density))
                            out.write("\t")
                            out.write(str(logdensity))
                            out.write("\t")
                            out.write(str(position))
                            out.write("\n")

if __name__ == '__main__':
    main()