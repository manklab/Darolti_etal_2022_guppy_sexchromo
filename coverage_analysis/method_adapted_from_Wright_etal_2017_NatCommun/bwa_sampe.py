#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
#BWA Version: 0.7.15
"This script creates and executes the commands to run bwa sampe for each sample and genome"
#==============================================================================
import argparse
import sys
import os
from subprocess import Popen
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("genomes_folder", type=str,
                    help="Infolder of genome assembly files")
parser.add_argument("reads_folder", type=str,
                    help="Infolder of fastq files and their suffix array indexes")
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
    return [os.path.join(infolder, f) for f in os.listdir(infolder)]

def exec_in_row(cmds):
    ''' Exec commands one after the other until finished.'''
    if not cmds:
        return  # empty list
    def done(p):
        return p.poll() is not None
    def success(p):
        return p.returncode == 0
    def fail():
        sys.exit(1)
    for task in cmds:
        print task
        p = Popen(task, shell=True)
        p.wait()
    if done(p):
            if success(p):
                print "done!"
            else:
                fail()
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

    #Creates list of genome files
    genomes = list_folder(args.genomes_folder)
    genomeslist = []
    for genome in genomes:
        if genome.endswith(".fa"):
            genomeslist.append(genome)
    print "Number of genome =", len(genomeslist)

    infiles = list_folder(args.reads_folder)
    #Creates dictionary of forward and reverse read pairs and the suffix array indexes for each sample
    filedictionary = defaultdict(list)
    for infile in infiles:
        name = os.path.basename(infile)
        uniqid = name.split("_")[0] # Assumes name of fastq files is of the format sample_R1.fastq or sample_R2.fastq, and that name of .sai files is of the format sample_R1_geneomeID.sai or sample_R2_genomeID.sai
        filedictionary[uniqid].append(infile)
    print "Number of samples =", len(filedictionary)

    #Cretaes list of bwa sampe commands to run
    commandstorun = []
    for genome in genomeslist:
        genomeID =  os.path.basename(genome).split("_")[0]
        for sample in filedictionary:
            for file in filedictionary[sample]: 
                # Assumes name of fastq files is of the format sample_R1.fastq or sample_R2.fastq, and that name of .sai files is of the format sample_R1_geneomeID.sai or sample_R2_genomeID.sai
                if file.endswith("R1.fastq"):
                    forward_fastq = file
                elif file.endswith("R2.fastq"):
                    reverse_fastq = file
                elif file.endswith(".sai") and genomeID in file and "R1" in file:
                    forward_sai = file
                elif file.endswith(".sai") and genomeID in file and "R2" in file:
                    reverse_sai = file
                outfile = file.split("_")[0]+"_"+genomeID+"_sampe.sam"
            command = ["bwa sampe "+genome+" "+forward_sai+" "+reverse_sai+" "+forward_fastq+" "+reverse_fastq+" -f "+outfile]
            commandstorun.append(command)

    print "Number of commands to run =", len(commandstorun)
    print "\n"
    exec_in_row(commandstorun)
    print "\n"

if __name__ == "__main__":
    main()