#!/usr/bin/python3

import argparse
import sys
import os
import subprocess
from Bio import SeqIO

"""
CHECK THE TE CONTENT OF NEORFS USING A TE LIBRARY
Author: Marie
Date: 12.09.2025
Edited: 12.09.2025
"""

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
    subprocess.call(f'blastn -query {PathToFasta} -subject {PathToTE} -perc_identity {Identity} -qcov_hsp_perc {Coverage} -evalue {Evalue} -outfmt 6 -out Blast_out', shell = True) 
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
    print(f"Number of neORFs without TE overlap: {NoTEcount}\nNumber of neORFs with TE overlap: {TEcount}")

    
Query ="/global/students/research/m_lebh01/DESwoMAN/FinalTestDESwoMAN/MarieCorrections/genomes_strat1/AK5/denovo_nucl.fa"
TEDB = "/global/scratch2/m_lebh01/TranscriptModelling/ReferenceNCBI/TE_DatabaseFlybase/FlyBase_TE.fa" 
run_blast_te(Query,TEDB, 80,0.01,80)