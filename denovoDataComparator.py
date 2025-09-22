#!/usr/bin/python3

import argparse
import sys
import os
import subprocess
from Bio import SeqIO

"""
Compare your neORFs with another de novo gene library (using blast) and get some stats of how similar they are
Author: Marie
Date: 18.09.2025
Edited: 18.09.2025
"""

def run_blastn(PathToFasta:str, PathToDB:str, Coverage:float, Evalue:float, Identity:float, Strand:str):
    """
    Generates the input file for Blast (i.e. all neORFs merged in one file)

    Parameters: 
    -PathToFasta (str): Path to a TE database. 
    -PathToTE (str): Path to the TE database 
    -Coverage (float): Qcovhsp for blast
    -Evalue(float): Evalue threshold for blast
    -Identity (float): Percent identity for blast
    -Strand(str): Only blast against the forward strand

    Returns:
    -Nothing
    """
    #Run the blast
    cmd = [
    "blastn",
    "-query", PathToFasta,
    "-subject", PathToDB,
    "-perc_identity", str(Identity),
    "-qcov_hsp_perc", str(Coverage),
    "-evalue", str(Evalue),
]
    # If a strand is specified
    if Strand:
        cmd += Strand.split() 
    cmd += [
        "-outfmt", "6 qseqid sseqid qstart qend qlen",
        "-out", "Blast_out"
    ]
    subprocess.run(cmd, check = True)
    #Safe all neORF nucleotides in a dictionary and fill an empty TE dictionary
    Nuc_dict = SeqIO.to_dict(SeqIO.parse(PathToFasta, "fasta"))
    Hitdict = {} #Store the hits
    for key in Nuc_dict:
        Hitdict[key] = "" #Initiate an (empty) entry for each neORF
    #Add the TE info from the Blast 
    with open("Blast_out") as Blast:
        for line in Blast:
            l = line.split("\t")
            Hitdict[l[0]] += l[1]
    NoHITcount = 0
    HITcount = 0
    for key, value in Hitdict.items():
        if not value:
            NoHITcount +=1
        else:
            HITcount +=1
    print("-------------------------------")
    print(f"Number of neORFs without a hit: {NoHITcount}\nNumber of neORFs with a hit: {HITcount}")
    print(f"{round((HITcount/len(Hitdict))*100,2)}% of the neORFs have a hit in the specified dataset {PathToDB}.")
    print("-------------------------------")
    os.remove("Blast_out")

def run_blastp(PathToFasta:str, PathToDB:str, Coverage:float, Evalue:float, Identity:float):
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
    cmd = [
    "blastp",
    "-query", PathToFasta,
    "-subject", PathToDB,
    "-qcov_hsp_perc", str(Coverage),
    "-evalue", str(Evalue),
    "-outfmt", "6 qseqid sseqid qstart qend qlen pident",
    "-out", "Blast_out"
]
    subprocess.run(cmd, check=True)
    #Safe all neORF nucleotides in a dictionary and fill an empty TE dictionary
    Nuc_dict = SeqIO.to_dict(SeqIO.parse(PathToFasta, "fasta"))
    Hitdict = {} #Store the hits
    for key in Nuc_dict:
        Hitdict[key] = "" #Initiate an (empty) entry for each neORF
    #Add the TE info from the Blast 
    with open("Blast_out") as Blast:
        for line in Blast:
            l = line.split("\t")
            if float(l[-1]) >= Identity:
                Hitdict[l[0]] += l[1]
    NoHITcount = 0
    HITcount = 0
    for key, value in Hitdict.items():
        if not value:
            NoHITcount +=1
        else:
            HITcount +=1
    print("-------------------------------")
    print(f"Number of neORFs without a hit: {NoHITcount}\nNumber of neORFs with a hit: {HITcount}")
    print(f"{round((HITcount/len(Hitdict))*100,2)}% of the neORFs have a hit in the specified dataset {PathToDB}.")
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
    parser = argparse.ArgumentParser(description="Get simple Fasta/Protein stats",epilog="-------------------------")
    parser.add_argument("--NeORF", help="Path to the NeORF fasta file", type=str)
    parser.add_argument("--denovo_db", help="Path to a second de novo gene dataset to compare it to", type=str)
    parser.add_argument("--evalue", help="Blast evalue treshold (default = 0.01)", type=float, default = "0.01")
    parser.add_argument("--perc_ident", help="Blast percent identity treshold (default = 80)", type=float, default = "0.01")
    parser.add_argument("--cov", help="Blast coverage treshold (default = 0.01)", type=float, default = "80")
    parser.add_argument("--type", help="Input (nuc for nucleotides = default, aa for amino acids)", type=str, default = "nuc")
    parser.add_argument("--plus", help="Only check the forward strand when checking for nucleotide hits",action="store_true")

    print("-------------------------\nCompare two de novo gene datasets (or other fasta files)\nV.1.0\nAuthor:Marie Lebherz\n-------------------------\n")

    #Read the input    
    args = parser.parse_args()
    NeORF = args.NeORF
    Hit = args.denovo_db
    Eval = args.evalue
    Pident = args.perc_ident
    Cov = args.cov
    Input = args.type
    if args.plus:
        Strand = "-strand plus"
    else:
        Strand = ""

    #Check if all is correct
    if not NeORF:
        print("No DESwoMAN neORF fasta file supplied. Exiting.")
        sys.exit()
    elif not Hit:
        print("No comparable de novo gene/protein database path supplied. Exiting.")
        sys.exit()

    #Start running the analysis
    print("Welcome to the de novo comparator stats tool :D\n")
    print(f"Your input dataframe is: {NeORF}")
    print(f"You comparison dataset is: {Hit}")
    print(f"Your databases are of type: {Input}")
    print(f"Your blast parameters are: Evalue of {Eval}, Coverage of {Cov}% and Percent Identitiy of {Pident}%!")
    if Input == "nuc":
        if Strand:
            print("You chose to only check for hits on the forward strand.")
        run_blastn(NeORF, Hit, Cov, Eval, Pident, Strand)
    elif Input == "aa":
        run_blastp(NeORF, Hit, Cov, Eval, Pident)
    else:
        print(f"{Input} is not a valid input. Please specify either nuc for nucleotides (dna) or aa for amino acids (protein)")
        sys.exit()
    print("Finished comparing the two datasets.")
    print("Goodybe :D")
        
# Example usage:
#Query ="/global/students/research/m_lebh01/DESwoMAN/FinalTestDESwoMAN/MarieCorrections/genomes_strat1/AK5/denovo_nucl.fa"
#TEDB = "/global/scratch2/m_lebh01/TranscriptModelling/ReferenceNCBI/TE_DatabaseFlybase/FlyBase_TE.fa" 
#run_blast_te(Query,TEDB, 0,10,0)

if __name__ == "__main__":
    main()