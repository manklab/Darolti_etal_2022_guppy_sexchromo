#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"This script creates a file with the start and stop positions of genes in scaffolds"
#==============================================================================
import argparse
import sys
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str,
                    help="File of M:F coverage fold change values for each chromosome block")
parser.add_argument("outfile", type=str,
                    help="Output file of M:F coverage fold change values for each chromosome block with adjusted coorinates to account for X chromosome inversions")
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

    with open(args.outfile, "w") as outfile:
        with open(args.infile, "r") as infile:
            for line in infile:
                line = line.rstrip()
                if line.startswith("Block"):
                    outfile.write(line)
                    outfile.write(",WindowStart_inversion\n")
                else:
                    start_pos = float(line.split(",")[2])
                    if line.startswith("NC_024342"):
                        if start_pos < 9900000:
                            new_start_pos = 9900000 - start_pos + 10800000
                        if start_pos >= 9900000 and start_pos < 20700000:
                            new_start_pos = start_pos - 9900000
                        if start_pos >= 20700000:
                            new_start_pos = start_pos
                        outfile.write(line)
                        outfile.write(",")
                        outfile.write(str(new_start_pos))
                        outfile.write("\n")
                    else:
                        outfile.write(line)
                        outfile.write(",")
                        outfile.write(str(start_pos))
                        outfile.write("\n")

if __name__ == '__main__':
    main()