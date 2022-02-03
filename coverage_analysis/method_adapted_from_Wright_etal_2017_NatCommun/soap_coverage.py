#!/usr/bin/python2.7
#SOAP.coverage Version: 2.7.7
"Script used to run SOAP.coverage commands"
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
parser.add_argument("sam_folder", type=str,
                    help="Infolder of sam files filtered for uniquely mapping reads")
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

    infiles = list_folder(args.sam_folder)
    #Creates dictionary of sam files
    fileslist = []
    for infile in infiles:
        if infile.endswith("_sampe_uniq.sam"): # Assumes name of sam file is of the format sample_genomeID_sampe_uniq.sam
            fileslist.append(infile)
    print "Number of sam files =", len(fileslist)

    #Creates list of soap.coverage commands to run
    commandstorun = []
    window = 50000 # For P. wingei analyses that used the P. wingei reference genome we used a window of 10000 instead
    for genome in genomeslist:
        genomeID =  os.path.basename(genome).split("_")[0]
        for infile in fileslist:
            if genomeID in infile:
                out_window = infile.split("sampe_uniq.sam")[0]+"soapcov_window.txt"
                out_soapcov = infile.split("sampe_uniq.sam")[0]+"soapcov.txt"
                command = ["soap.coverage -sam -cvg -i "+infile+" -onlyuniq -refsingle "+genome+" -p 12 -window "+out_window+" "+window+" -o "+out_soapcov]
                commandstorun.append(command)

    print "Number of commands to run =", len(commandstorun)
    print "\n"
    exec_in_row(commandstorun)
    print "\n"

if __name__ == "__main__":
    main()