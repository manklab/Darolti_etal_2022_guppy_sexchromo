#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"Add row numbers to file"
#==============================================================================
import argparse
import sys
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str,
                    help="An input file")
parser.add_argument("outfile", type=str,
                    help="An output file with row numbers")
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

    firstline = 0
    countline =0
    with open(args.outfile, "w") as out:
        with open(args.infile, "r") as infile:
            for line in infile:
                if firstline == 0:
                    out.write("Numbers,")
                    out.write(line)
                else:
                    countline += 1
                    out.write(str(countline))
                    out.write(",")
                    out.write(line)
                firstline += 1

if __name__ == '__main__':
    main()