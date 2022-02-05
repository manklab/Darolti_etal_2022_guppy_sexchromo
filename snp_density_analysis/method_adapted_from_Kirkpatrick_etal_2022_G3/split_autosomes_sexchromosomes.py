#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"This script is used to split the snpdensity_fc.txt file into separate files for the autosomes and the sex chromosomes"
#==============================================================================
import argparse
import sys
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str,
                    help="The snpdensity_fc.txt file")
parser.add_argument("output_autosomes", type=str,
                    help="Output file containing only the autosomes")
parser.add_argument("output_sexchromosomes", type=str,
                    help="Output file containing only the sex chromosomes")
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

    with open(args.output_autosomes, "w") as auto:
        with open(args.output_sexchromosomes, "w") as sex:
            with open(args.infile, "r") as infile:
                for line in infile:
                    if line.startswith("Chromosome"):
                        auto.write(line)
                        sex.write(line)
                    else:
                        chromosome = line.split(",")[1]
                        # The sex chromosome is labelled "NC_024342.1" for the P. reticulata genomes and "8" for the X. maculatus and P. wingei genomes.
                        if chromosome == "NC_024342.1":
                            sex.write(line)
                        else:
                            auto.write(line)

if __name__ == '__main__':
    main()