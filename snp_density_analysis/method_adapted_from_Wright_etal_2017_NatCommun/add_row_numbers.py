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
parser.add_argument("snp_fc", type=str,
                    help="A file of scaffolds and SNPs")
parser.add_argument("outfile", type=str,
                    help="File with genes and their position in the genome")
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

    count = 0
    countline =0
    with open(args.outfile, "w") as out:
        with open(args.snp_fc, "r") as fc:
            for line in fc:
                if count == 0:
                    out.write("Numbers,")
                    out.write(line)
                else:
                    countline += 1
                    out.write(str(countline))
                    out.write(",")
                    out.write(line)
                count += 1

if __name__ == '__main__':
    main()