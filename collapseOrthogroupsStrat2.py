#!/usr/bin/python3
"""
COLLAPSE THE STRATEGY2 OUTPUT
Author: Marie
Date: 26.02.2026
"""
import argparse
import sys
from Bio import SeqIO
import random

def get_orthogroups(PathToDESwoMAN:str)->dict:
    """
    Read the orthogrups File and return the groups

    Parameters:
    -PathToDESwoMAN (str): Path to DESwoMAN output folder

    Output:
    -dict: Dictionary with the orthogroups
    """
    Groups = {}
    GroupFile = open(f"{PathToDESwoMAN}Orthogroup_output_step3.txt", "r")
    for line in GroupFile:
        l = line.strip().split(":")
        Groups[l[0]] = l[1].split(",")
    return Groups

def collapse_by_ends(PathToDESwoMAN:str,Out:str)->dict:
    """
    Get the whole genomic locus that an orthogroup is spanning (i.e. earliest start, latest end)
    [Note!] Some of the resulting sequences will not be ORFs anymore!

    Parameters:
    -PathToDESwoMAN (str): Path to DESwoMAN output folder
    -Out(str): Path to the output folder to save your data
    """
    #Make a dict but reverse to above
    Groups = {}
    GroupFile = open(f"{PathToDESwoMAN}Orthogroup_output_step3.txt", "r")
    for line in GroupFile:
        l = line.strip().split(":")
        for neORF in l[1].split(","):
            Groups[neORF] = l[0]
    GroupFile.close()
    #Now parse the infoFile
    InfoFile =  open(f"{PathToDESwoMAN}information_file.txt", "r")
    InfoFile.readline()
    CoordinateDict = {}
    for line in InfoFile:
        l = line.strip().split(",")
        Group = Groups[l[9]]
        if Group in CoordinateDict:
            CoordinateDict[Group][3].append(int(l[-3]))
            CoordinateDict[Group][4].append(int(l[-2])+1)
        else:
            CoordinateDict[Group] = [l[-1], l[1],l[2], [int(l[-3])], [int(l[-2])+1]]

    Out = open(f"{Out}/genomic_position_per_orthogroup.bed", "w")
    for key,value in CoordinateDict.items():
        Out.write(f"{value[1]}\t{str(min(value[3]))}\t{str(max(value[4]))}\t{key}|{value[0]}\t.\t{value[2]}\n")


def collapse_by_longest(Groups:dict)->dict:
    """
    Select the (first) longest ORF per orthogroup

    Parameters:
    -Groups(dict): Orthogroups in a dictionary where key is the group ID and values are the neORF IDs from that group

    Returns:
    -dict: Dictionary where for each group the (first) longest ORF is chosen
    """
    LongestDict = {}
    for group in Groups:
        LongestLength,LongestID = 0,""
        for neORF in Groups[group]:
            Length = int(neORF.split("_")[-1])-int(neORF.split("_")[-2])
            if Length > LongestLength: #if its the same length the first will be chosen
                LongestLength = Length
                LongestID = neORF
        LongestDict[LongestID] = group
    return LongestDict

def select_random(Groups:dict)->dict:
    """
    Select a random ORF per orthogroup

    Parameters:
    -Groups(dict): Orthogroups in a dictionary where key is the group ID and values are the neORF IDs from that group

    Returns:
    -dict: Dictionary where for each group an ORF is randomly chosen
    """
    RandomDict = {}
    for group in Groups:
        RandomDict[random.sample(Groups[group],1)[0]] = group
    return RandomDict

def redoInfoFile(ReprDict:dict, PathToDESwoMAN:str, Out:str):
    """
    Redo the information file

    Parameters:
    -ReprDict (dict): Dictionary with one representative Sequence per orthogroup
    -PathToDESwoMAN (str): Path to the DESwoMAN output
    -Out(str): Path to the output folder to save your data
    """
    InfoFile = open(f"{PathToDESwoMAN}information_file.txt")
    NewInfoFile = open(f"{Out}/information_file_collapsed.txt", "w")
    NewInfoFile.write(InfoFile.readline().replace("neORF_candidate", "orthogroup_ID"))
    for line in InfoFile:
        l = line.strip().split(",")
        if l[9] in ReprDict:
            NewInfoFile.write(line.replace(l[9],f"{ReprDict[l[9]]} representative={l[9]}"))
    NewInfoFile.close()
    return None

def redoFastaFile(ReprDict:dict, PathToDESwoMAN:str,Out:str):
    """
    Redo the nucleotide and protein fasta file

    Parameters:
    -ReprDict (dict): Dictionary with one representative Sequence per orthogroup
    -PathToDESwoMAN (str): Path to the DESwoMAN output
    -Out(str): Path to the output folder to save your data
    """
    Nuc = SeqIO.to_dict(SeqIO.parse(f"{PathToDESwoMAN}denovo_nucl.fa" , "fasta"))
    NewNuc = open(f"{Out}/denovo_nucl_collapsed.fa", "w")
    Prot = SeqIO.to_dict(SeqIO.parse(f"{PathToDESwoMAN}denovo_protein.fa" , "fasta"))
    NewProt = open(f"{Out}/denovo_protein_collapsed.fa", "w")
    for neORF in ReprDict:
        #if neORF not in Nuc:
            #print("LALALA:",neORF)
            #continue
        NewNuc.write(f">{ReprDict[neORF]} representative={neORF}\n{Nuc[neORF].seq}\n")
        NewProt.write(f">{ReprDict[neORF]} representative={neORF}\n{Prot[neORF].seq}\n")
    NewNuc.close()
    NewProt.close()
    return None

def collapseORF(PathToDESwoMAN:str, collapse_choice:str,Out:str):
    """
    Run the whole pipeline

    Parameters:
    -PathToDESwoMAN (str): Path to the DESwoMAN output folder
    -collapse_choice(str): how to collapse the orthogroups (random or longest)
    -Out(str): Path to the output folder to save your data
    """
    #Get the orthogroups
    Groups = get_orthogroups(PathToDESwoMAN)
    #Collapse by chosen parameter
    if collapse_choice == "longest":
        ReprDict = collapse_by_longest(Groups)
    elif collapse_choice == "random":
        ReprDict = select_random(Groups)
    #Get the collapsed files
    redoInfoFile(ReprDict, PathToDESwoMAN, Out)
    redoFastaFile(ReprDict, PathToDESwoMAN, Out)

def main():
    """
    Main function
    """
    #Initialize all arguments
    parser = argparse.ArgumentParser(description="Collapse the strategy 2 output.",epilog="-------------------------")
    parser.add_argument("--deswoman", help="Path to the DESwoMAN folder", type=str)
    parser.add_argument("--choice", help="Strategy to find a representative sequence per orthogroup", type=str, default = "longest")
    parser.add_argument("--outpath", help="Path to the folder where you want to save your data", type=str, default = ".")
    print("-------------------------\nCollapse Strat2 Output\nV.1.0\nAuthor:Marie Lebherz\n-------------------------\n")

    #Read the Parameters    
    args = parser.parse_args()
    DESwoMANPath = args.deswoman
    Choice = args.choice
    Out = args.outpath

    #Check if all is correct
    if not DESwoMANPath:
        print("[Error:] No DESwoMAN path supplied. Exiting.")
        sys.exit()
    elif Choice not in ["longest", "random", "locus"]:
        print(f"[Error:] Invalid collapsing choice. Please chose from: longest, random. You chose: {Choice}")
        sys.exit()

    print(f"Your strategy to collapse the Orthogroups is: [{Choice}]")
    if Choice in ["longest", "random"]:
        collapseORF(DESwoMANPath,Choice,Out)
    else:
        print("\n!Info!\nNote that because of your choice some sequences might not be ORFs. The earliest start and latest end will be used!\nResults will be  reported in a bedfile!")
        collapse_by_ends(DESwoMANPath,Out)
    print("\nAll your files are now collapsed. Finished :) ")
    
if __name__ == "__main__":
    main()