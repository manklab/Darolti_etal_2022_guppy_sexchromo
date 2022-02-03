#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
#==============================================================================
import argparse
import sys
import os
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("gtf", type=str,
                    help="A gtf file")
parser.add_argument("outfile", type=str,
                    help="A gene cooridnates file")
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

    with open(args.outfile, "w") as out:
        out.write("chromosome,gene,start,end\n")
        with open(args.gtf, "r") as gtf:
            for line in gtf:
                if line.split()[2] == "gene":
                    chromosome = line.split()[0]
                    start = line.split()[3]
                    end = line.split()[4]
                    gene = line.split('gene_id "')[1].split('";')[0]
                    out.write(str(chromosome))
                    out.write(",")
                    out.write(gene)
                    out.write(",")
                    out.write(str(start))
                    out.write(",")
                    out.write(str(end))
                    out.write("\n")

if __name__ == '__main__':
    main()