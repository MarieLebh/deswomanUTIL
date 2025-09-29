#!/usr/bin/python3

import argparse
import sys
import os

"""
GENERATE AN OUTPUT BED FILE WITH ALL DESwoMAN NeORFs
Get a bed file with either neORF or transcript coordinates
Author: Marie
Date: 10.09.2025

usage: getBEDfromOut.py [-h] [--deswoman DESWOMAN] [--gtf GTF] [--outname OUTNAME] [--choice CHOICE]
"""

def generate_final_file_orf(DESwoMAN:str, Outname:str):
    """
    Generate a final file from the DESwoMAN output
 
    Parameters: 
    - DESwoMAN (str): Path to the DESwoMAN info file
    - Outname (str): Path to the output file

    Returns:
    - nothing
    """
    #Open all files/extract all information
    DESWOMAN = open(DESwoMAN, "r")
    Out = open(Outname + "_orf.bed", "w")
    #Now parse the DESwoMAN output file for information (all lines are saved as gff format lnes)
    DESWOMAN.readline()
    for line in DESWOMAN:
        l = line.split(",")
        Out.write(f"{l[1]}\t{l[-3]}\t{str(int(l[-2])+1)}\t{l[-4]}\t.\t{l[2]}\n")
    #Close everything again
    DESWOMAN.close()
    Out.close()

def generate_final_file_transcript(GTF:str, Outname:str, DESwoMAN:str):
    """
    Generate a final file from the DESwoMAN output 
 
    Parameters: 
    - GTF (str): Path to the transcriptome gtf file
    - Outname (str): Path to the output file
    - DESwoMAN (str): Path to the DESwoMAN info file

    Returns:
    - nothing
    """
    #Get all de novo ORFs:
    idList = []
    DESWOMAN = open(DESwoMAN, "r")
    DESWOMAN.readline() #skip the first line
    for line in DESWOMAN:
        l = line.split(",")
        idList.append(l[5])
    DESWOMAN.close()

    #Open all files/extract all information
    GTFfile = open(GTF, "r")
    Out = open(Outname + "_transcript.bed", "w")
    #Now parse the gtf output file for information 
    for line in GTFfile:
        if line[0] != "#":
            l = line.split("\t")
            if l[2] == "transcript" and l[8].split(";")[1].split(" ")[2][1:-1] in idList:
                Out.write(f"{l[0]}\t{str(int(l[3])-1)}\t{l[4]}\t{l[8].split(';')[1].split(' ')[2][1:-1]}\t.\t{l[6]}\n")
    Out.close()

def main():
    """
    Main function

    Parameter:
    - none

    Returns:
    - Nothing
    """
    #Initialize all arguments
    parser = argparse.ArgumentParser(description="Get high confidence neORFs from the DESwoMAN output",epilog="-------------------------")
    parser.add_argument("--deswoman", help="Path to the DESwoMAN info file", type=str)
    parser.add_argument("--gtf", help="Path to the transcriptome gtf file", type=str)
    parser.add_argument("--outname", help="Name of the output file (no file extension)", type=str, default = "DESwoMAN")
    parser.add_argument("--choice", help="What positions are wanted in the output bedfile (transcript, neORF or both)", type=str, default = "both")
    print("-------------------------\nGet a DESwoMAN BED file\nV.1.0\nAuthor:Marie Lebherz\n-------------------------\n")

    #Read the input    
    args = parser.parse_args()
    GTF =args.gtf
    DESwoMANPath = args.deswoman
    Outname = args.outname
    Choice = args.choice

    #Check if all is correct
    if not DESwoMANPath:
        print("[Error:] No DESwoMAN path supplied. Exiting.")
        sys.exit()
    elif (Choice == "both" or Choice == "transcript") and not GTF:
        print("[Error:] No transcriptome assembly path supplied. Exiting.")
        sys.exit()

    #Start running the analysis
    print("Welcome to the deswomanUTIL bed converter...\n\nStarting the file transformation now!\n")
    if Choice == "neORF":
        generate_final_file_orf(DESwoMANPath, Outname)
        print("Your neORF bedfile was successfully transformed!!")
    elif Choice == "transcript":
        generate_final_file_transcript(GTF, Outname, DESwoMANPath)
        print("Your transcript bedfile was successfully transformed!!")
    elif Choice == "both":
        generate_final_file_orf(DESwoMANPath, Outname)
        generate_final_file_transcript(GTF, Outname, DESwoMANPath)
        print("Your files (neORF + transcript) were successfully transformed!!")
    print("Goodybe :D")
    
if __name__ == "__main__":
    main()