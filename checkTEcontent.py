#!/usr/bin/python3

import argparse
import sys
import os
import subprocess
from Bio import SeqIO
import matplotlib.pyplot as plt

"""
CHECK THE TE CONTENT OF NEORFS USING A TE LIBRARY
Author: Marie
Date: 12.09.2025
Edited: 12.09.2025
"""

def count_items_higherX(llist:list, thresh:float)-> float: 
    """
    Takes the Blast tabular output and calculates the percentage of the neORFs CDS that is covered by a TE. Prints summary stats. 

    Parameters: 
    -llist (list): List of TE coverage percentages
    -thresh (float): Theshold

    Returns:
    - float: Percentage of sequences (from the list) passing that theshold
    """
    count = 0
    for i in llist:
        if i >= thresh:
            count +=1
    return (count/len(llist))*100

def merge_intervals(intervals:list)->list:
    """
    Takes the Blast tabular output and calculates the percentage of the neORFs CDS that is covered by a TE. Prints summary stats. 

    Parameters: 
    -intervals (list): List of intervals that may or may not overlap with each other

    Returns:
    -list: A list with concatenated non overlapping/nested intervals
    """
    if not intervals:
        return []
    
    # Sort intervals by starting point
    intervals.sort(key=lambda x: x[0])
    
    merged = [intervals[0]]
    
    for current in intervals[1:]:
        prev = merged[-1]
        
        # If intervals overlap or are nested
        if current[0] <= prev[1]:
            merged[-1][1] = max(prev[1], current[1])  # extend the interval
        else:
            merged.append(current)
    
    return merged


def get_te_overlap():
    """
    Takes the Blast tabular output and calculates the percentage of the neORFs CDS that is covered by a TE. Prints summary stats. 

    Parameters: 
    -Nothing

    Returns:
    -Nothing
    """
    OverlapDict = {}
    #Extract the overlapping TEs
    with open("Blast_out", "r") as Blast:
        for line in Blast:
            l = line.strip().split("\t")
            if f"{l[0]}#{l[4]}" not in OverlapDict:
                interval = sorted([int(l[2]),int(l[3])])
                OverlapDict[f"{l[0]}#{l[4]}"] = [interval] 
            else:
                interval = sorted([int(l[2]),int(l[3])])
                OverlapDict[f"{l[0]}#{l[4]}"] += [interval] 
    #Now get the percentage of overlap
    Overlaps = []
    for key, value in OverlapDict.items():
        NR_ints = merge_intervals(value)
        TELength = 0
        for interval in NR_ints:
            TELength += int(interval[1]) - int(interval[0]) +1
        Overlaps.append((TELength/int(key.split("#")[1]))*100)
    print(f"Out of all neORFs with a TE hit: {round(count_items_higherX(Overlaps,80),2)} % (nearly) fully overlap with a TE (>= 80%)")
    print(f"Out of all neORFs with a TE hit: {round(100 - count_items_higherX(Overlaps,80),2)} % overlap < than 80 %")

def run_blast_te(PathToFasta:str, PathToTE:str, Coverage:float, Evalue:float, Identity:float):
    """
    Generates the input file for Blast (i.e. all neORFs merged in one file)

    Parameters: 
    -PathToFasta (str): Path to a TE database. 
    -PathToTE (str): Path to the TE database 
    -Coverage (float): Qcovhsp for blast
    -Evalue(float): Evalue threshold for blast
    -Identity (float): Percent identity for blast

    Returns:
    -Nothing
    """
    #Run the blast
    subprocess.call(f'blastn -query {PathToFasta} -subject {PathToTE} -perc_identity {Identity} -qcov_hsp_perc {Coverage} -evalue {Evalue} -outfmt "6 qseqid sseqid qstart qend qlen" -out Blast_out', shell = True) 
    #Safe all neORF nucleotides in a dictionary and fill an empty TE dictionary
    Nuc_dict = SeqIO.to_dict(SeqIO.parse(PathToFasta, "fasta"))
    TE_dict = {}
    for key in Nuc_dict:
        TE_dict[key] = ""
    #Add the TE info from the Blast 
    with open("Blast_out") as Blast:
        for line in Blast:
            l = line.split("\t")
            TE_dict[l[0]] += l[1]
    NoTEcount = 0
    TEcount = 0
    for key, value in TE_dict.items():
        if not value:
            NoTEcount +=1
        else:
            TEcount +=1
    print("-------------------------------")
    print(f"Number of neORFs without TE overlap: {NoTEcount}\nNumber of neORFs with TE overlap: {TEcount}")
    print("-------------------------------")
    get_te_overlap()
    print("-------------------------------")
    os.remove("Blast_out")

def main():
    """
    The main function to run the program

    Parameters: 
    None (need to be parsed to the commandline)

    Returns:
    Nothing
    """
    #Initialize all arguments
    parser = argparse.ArgumentParser(description="Get simple TE stats",epilog="-------------------------")
    parser.add_argument("--NeORF", help="Path to the NeORF fasta file", type=str)
    parser.add_argument("--TE_db", help="Path to the TE database (fasta)", type=str)
    parser.add_argument("--evalue", help="Blast evalue treshold (default = 0.01)", type=float, default = "0.01")
    parser.add_argument("--perc_ident", help="Blast percent identity treshold (default = 80)", type=float, default = "0.01")
    parser.add_argument("--cov", help="Blast coverage treshold (default = 0.01)", type=float, default = "80")
    print("-------------------------\nCheck TE content of neORFs\nV.1.0\nAuthor:Marie Lebherz\n-------------------------\n")

    #Read the input    
    args = parser.parse_args()
    NeORF = args.NeORF
    TE = args.TE_db
    Eval = args.evalue
    Pident = args.perc_ident
    Cov = args.cov

    #Check if all is correct
    if not NeORF:
        print("No DESwoMAN neORF fasta file supplied. Exiting.")
        sys.exit()
    elif not TE:
        print("No TE database path supplied. Exiting.")
        sys.exit()

    #Start running the analysis
    print("Welcome to the TE stats tool")
    run_blast_te(NeORF, TE, Cov, Eval, Pident)
    print("Sucessfully generated simple TE stats.")
    print("Goodybe :)")
        
# Example usage:
#Query ="/global/students/research/m_lebh01/DESwoMAN/FinalTestDESwoMAN/MarieCorrections/genomes_strat1/AK5/denovo_nucl.fa"
#TEDB = "/global/scratch2/m_lebh01/TranscriptModelling/ReferenceNCBI/TE_DatabaseFlybase/FlyBase_TE.fa" 
#run_blast_te(Query,TEDB, 0,10,0)

if __name__ == "__main__":
    main()