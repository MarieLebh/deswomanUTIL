#!/usr/bin/python3

import argparse
import sys
import os

"""
Fai Converter (Fai to Chrom Length bed3)
Author: Marie
Date: 30.04.25
Run like this: python3 Fai_converter.py --in FaiFile (--out Outname)
"""

#initiate the input 
parser = argparse.ArgumentParser(description='The tool will create a bed(3)file with these columns: ChromName, Start(=0), Stop(=ChromLength)',epilog="-------------------------")
parser.add_argument("--inp", help="Path to the Genome index file (.fai)", type=str)
parser.add_argument("--out", help="Path to the Output file (.bed). Default: 'PathToInput_ChromLength.bed'", type=str, default = "Test")

def index_to_chromLengthbed(path, OutPath):
    #get the chrom lengths in a bedfile from the index
    if not OutPath:
        Name = f"{path.strip('.fna.fai')}_ChromLength.bed"
    else:
        Name = OutPath
    IndexFile = open(path, "r")
    Out = open(Name, "w")
    for line in IndexFile:
        Chrom = line.split("\t")[0]
        ChromLength = line.split("\t")[1]
        Start = 0
        Out.write(f"{Chrom}\t{Start}\t{ChromLength}\n")
    IndexFile.close()
    Out.close()


if __name__ == "__main__":
    print("-------------------------\nCreate Chromosome Length Bedfiles from Fasta index(.fai)\nV.1.0\nAuthor:Marie Lebherz\n-------------------------\n")
    File = parser.parse_args().inp
    if not File:
        print("No argument supplied. Exiting...")
        sys.exit()
    if os.path.isfile(File) != True:
        print("The file you supplied doesn't exist.")
        sys.exit()
    if File[-4:] != ".fai":
        print(f"Please supply a .fai file!\nThe problematic file is: {File}")
        sys.exit()
    OutputName = parser.parse_args().out
    print("Start extracting the information...")
    index_to_chromLengthbed(File, OutputName)
    print("Finished extracting the chromosome length file.")
    print("Goodybe ;)")


