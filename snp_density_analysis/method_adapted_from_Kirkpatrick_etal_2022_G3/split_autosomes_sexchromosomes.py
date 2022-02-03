#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
from collections import OrderedDict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str,
                    help="")
parser.add_argument("out_auto", type=str,
                    help="")
parser.add_argument("out_sex", type=str,
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
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

    with open(args.out_auto, "w") as auto:
        with open(args.out_sex, "w") as sex:
            with open(args.infile, "r") as infile:
                for line in infile:
                    if line.startswith("Block"):
                        auto.write(line)
                        sex.write(line)
                    else:
                        chromo = line.split(",")[1]
                        if chromo == "NC_024342.1":
                            sex.write(line)
                        else:
                            auto.write(line)

if __name__ == '__main__':
    main()