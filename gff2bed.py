#!/usr/bin/python3

import argparse
import sys
import os

"""
GFF to Bed
Author: Marie
Date: 30.04.25
Run like this: python3 GFF_converter.py --inp GFFFile (--out Outname (no suffix)) (--type Type (e.g. gene/mRNA/transcript, ...))
"""

parser = argparse.ArgumentParser(description="The tool will create a bed(6)file from a gff file",epilog="-------------------------")
parser.add_argument("--inp", help="Path to the Genome index file (.fai)", type=str)
parser.add_argument("--out", help="Path to the Output file (.bed). Default: 'PathToInput_ChromLength.bed'", type=str, default = "Test")
parser.add_argument("--type", help="Path to the Output file (.bed). Default: 'PathToInput_ChromLength.bed'", type=str, default = "mRNA")

def gff2bed6(gff, outname, seq_type):
    File = open(gff, "r")
    Bed = open(outname + ".bed", "w")
    for line in File:
        if line[0] == "#": #ignore hashtags
            continue
        else: #Get all the informaton from the different columns
            l = line.strip().split("\t") 
            if l[2] == seq_type:
                Bed.write(f"{l[0]}\t{str(int(l[3]) -1 )}\t{l[4]}\t{l[8].split(';')[0]}\t{l[5]}\t{l[6]}\n")
    Bed.close()

if __name__ == "__main__":
    print("-------------------------\nGFF/GTF to Bed6\nV.1.0\nAuthor:Marie Lebherz\n-------------------------\n")
    File = parser.parse_args().inp
    if not File:
        print("No argument supplied. Exiting...")
        sys.exit()
    if os.path.isfile(File) != True:
        print("The file you supplied doesn't exist.")
        sys.exit()
    if File[-4:] not in [".gtf", ".gff"]:
        print(f"Please supply a file with a correct suffix (gff or gtf)!\nThe problematic file is: {File}")
        sys.exit()
    OutputName = parser.parse_args().out
    Type = parser.parse_args().type
    print("Start extracting the information...")
    gff2bed6(File, OutputName, Type)
    print("Finished extracting the chromosome length file.")
    print("Goodybe ;)")
