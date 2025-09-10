from Bio import SeqIO
import os
import subprocess
from collections import Counter
import argparse
import sys
from Bio import Phylo
from io import StringIO

"""
GET A SUBSET OF HIGH CONFIDENCE DE NOVO ORFs
This script filters the DESwoMAN output based on the noncoding homologs and a phylogeny
Author: Marie
Date: 16.05.25
Edited: 04.08.25

How to run:
usage: filterDESwoMAN.py [-h] [--te_db TE_DB] [--te_cov TE_COV] [--te_idt TE_IDT] [--te_eval TE_EVAL] [--ortho ORTHO] [--deswoman DESWOMAN] [--rna_check] [--te_check]
                         [--tr_db TR_DB] [--tr_cov TR_COV] [--tr_idt TR_IDT] [--tr_eval TR_EVAL] [--tr_strand TR_STRAND] [--out OUT] [--tree TREE]
                         [--species_file SPECIES_FILE] [--accepted_mutations ACCEPTED_MUTATIONS] [--frameshift_score FRAMESHIFT_SCORE]


CRITERIA FOR FILTERING:

    Less than 80% identity and coverage overlap with a TE (curated Drosophila database)
    At least one noncoding homolog in an outgroup species
    If there is a species with no coding mutation present that species is added to the coding species.
"""

def get_populations(SpeciesList:str)->list:
    """
    Read in a text file with one population/id on each line

    Parameters: 
    - Filepath (str): Path to the file with the sample ids

    Returns:
    -list: list of lines
    """
    Strains = []
    Info = open(SpeciesList, "r")
    for line in Info:
        l = line.strip()
        Strains.append(l)
    return Strains

def find_common_node(newick:str, species_list:list)->str:
    """
    Input a tree and a list of species and return their most recent common ancestor

    Parameters: 
    - newick (str): Path to the newick file
    - species_list (list): List of species names (need to all be in the newick file)

    Returns:
    -list: list of lines
    """
    tree = Phylo.read(newick, "newick")
    # Normalize species names (remove dots if needed, optional depending on your data)
    species_list = [s.replace(" ", "").strip() for s in species_list]
    # Find the MRCA
    mrca = tree.common_ancestor(species_list)
    # Return the name of the MRCA node
    return mrca.name

def get_outgroup_from_clade_name(tree:str, clade_name:str)->list:
    """
    Input a clade name (needs to be included in the tree), it will return all terminal branches who are outgroups

    Parameters: 
    - tree (str): Path to the newick file
    -clade_name (str): The name of the clade of interest. 

    Returns:
    -list: list of lines
    """
    # Step 1: Find the target clade by name
    tree = Phylo.read(tree, "newick")
    target_clade = None
    for clade in tree.find_clades():
        if clade.name == clade_name:
            target_clade = clade
            break
    if not target_clade:
        raise ValueError(f"Clade '{clade_name}' not found in the tree.")
    # Step 2: Get all terminal names in the target clade (ingroup)
    ingroup = {t.name for t in target_clade.get_terminals()}
    # Step 3: Get all terminal names in the whole tree
    all_terminals = {t.name for t in tree.get_terminals()}
    # Step 4: Outgroup = all - ingroup
    outgroup = sorted(all_terminals - ingroup)
    return outgroup

def read_orthogroups_info(OrthoPath:str, PathToDESwoMAN:str, SpeciesList:str)->dict:
    """
    Make a dictionary with all coding homologs

    Parameters: 
    -OrthoPath(str): Path to the "Orthogroups.txt" foöe
    -PathToDESwoMAN(str): Path to the DESwoMAN folder
    - SpeciesList: File with all species (.txts)

    Returns:
    -dict: Dictionary of coding homologs. 
    """
    #For the species instead of looking at the number of species in an orthogroup also add the tree
    Specs = get_populations(SpeciesList)
    if OrthoPath == "NotSupplied":
        CodingDict = {}
        for line in Specs:
            Prots =  SeqIO.to_dict(SeqIO.parse(f"{PathToDESwoMAN}/{line}/denovo_protein.fa", "fasta"))
            for key in Prots:
                CodingDict[key+"_" + line] = [line[:4]]
     
    else:
        Orthogroups = {}
        Info = open(OrthoPath, "r")
        for line in Info:
            l = line.strip()
            OG_ID = l.split(":")[0]
            OGs = l.split(":")[1][1:].split(" ")
            Spec = []
            for i in OGs:
                if i:
                    Spec.append(i.split("_")[-1][:4])

            Species_Checker = set(Spec)
            Orthogroups[OG_ID] = [OGs, Species_Checker]  
        Info.close()
        CodingDict = {}
        for key, value in Orthogroups.items():
            Seqs = value[0]
            for i in Seqs:
                CodingDict[i] = list(value[1])  #Create the final dictionary

    return CodingDict

def get_acceptable_noncoding_homologs(PathToDESwoMAN:str, CodingDict:dict, SpeciesList:str, mutations:list, frame:float)->tuple[dict,dict]:
    """
    Make a dictionary with all coding homologs

    Parameters: 
    -PathToDESwoMAN(str): Path to the DESwoMAN folder
    -CodingDict (dict): Dictionary with the neORFs
    -SpeciesList (str): File with all species (.txts)
    -mutations (list): List with all accepted mutations
    -frame (float): Maximum percentage of the sequence not affected by frameshift to still count as noncoding

    Returns:
    -tuple[dict,dict]: Dictionary of coding homologs (updated), Accepted noncoding homologs (= with mutation)
    """
    x = get_populations(SpeciesList)
    validated_count = 0
    unvalidated_count = 0
    not_present_count = 0
    AcceptedHomologs = {}

    for line in x:
        HomologInfo = open(f"{PathToDESwoMAN}/{line}/Table_output_step3.txt", "r")
        HomologInfo.readline()
        for li in HomologInfo:
            l = li.strip().split(",")
            neORFID, TargetSpecies = f"{l[0]}_{line}", l[1][0:4]
            start, stop, Indels,perc_seq_not_affected_by_frameshift, substitutions, premature_stop, transcription_status = l[4], l[5], l[6], l[7], l[8], l[9], l[11]

            #Really comblicated way deal with each mutation
            startcheck = "P" if "start" in mutations else "placeholder"
            stopcheck = "P" if "stop" in mutations else "placeholder"
            premature_stopcheck = "A" if "premature_stop" in mutations else "placeholder"
            completecheck = "complete" if "complete" in mutations else "placeholder"
            partialcheck = "partial" if "partial" in mutations else "placeholder"
            antisensecheck = "reverse" if "reverse" in mutations else "placeholder"

            if l[2] == "P" and start != "NA": #Only include present homologs. 
                if start != startcheck or stop != stopcheck  or premature_stop != premature_stopcheck or float(perc_seq_not_affected_by_frameshift) <= frame or (transcription_status != completecheck and transcription_status !=partialcheck and transcription_status !=antisensecheck):#transcription_status != "complete": #Change here (!)
                #if start != "P" or stop != "P"  or premature_stop != "A" or int(Indels) > 0 or transcription_status != "complete":
                    #Check if any mutation is suggesting that the homolog is non coding
                    validated_count += 1
                    if neORFID in AcceptedHomologs:
                        AcceptedHomologs[neORFID] += [TargetSpecies]
                    else: 
                        AcceptedHomologs[neORFID] = [TargetSpecies]
                else:
                    CodingDict[neORFID] += [TargetSpecies]
                    CodingDict[neORFID] = list(set(CodingDict[neORFID])) #In case there is no noncoding mutation add that species to the potentially coding ones

                    #If all features suggest a coding status then exclude this target as an accepted homolog
                    unvalidated_count += 1
                    if neORFID not in AcceptedHomologs:
                        AcceptedHomologs[neORFID] = [] #If no entry is present add an empty list
            else:
                #If no hit was present at all then also exclude this target as an accepted homolog
                not_present_count += 1
                if neORFID not in AcceptedHomologs:
                    AcceptedHomologs[neORFID] = [] #If no entry is present add an empty list
        HomologInfo.close()
    print(f"Number of valid homologs:{validated_count}\nNumber of possibly coding homologs:{unvalidated_count}\nNumber of absent homologs:{not_present_count}")
  
    return AcceptedHomologs, CodingDict

def generateBlastInput(PathToDESwoMAN:str, SpeciesList:list)->list:
    """
    Generates the input file for Blast (i.e. all neORFs merged in one file)

    Parameters: 
    -PathToDESwoMAN(str): Path to the DESwoMAN folder
    -SpeciesList: File with all species (.txts)

    Returns:
    -list: List of ORFs
    """
    x = get_populations(SpeciesList)
    ORF_list = []
    BlastFile = open("neORFs_Nuc_all.fa", "w")
    for line in x: 
        Nuc_dict = SeqIO.to_dict(SeqIO.parse(f"{PathToDESwoMAN}/{line}/denovo_nucl.fa", "fasta"))
        for key, value in Nuc_dict.items():     
            BlastFile.write(f">{key}_{line}\n{value.seq}\n")
            ORF_list.append(key + "_" + line)
    BlastFile.close()
    return ORF_list

def run_blast_te(PathToTE:str, Coverage:float, Evalue:float, Identity:float, Strand:str):
    """
    Generates the input file for Blast (i.e. all neORFs merged in one file)

    Parameters: 
    -PathToTE (str): Path to the TE database 
    -Coverage (float): Qcovhsp for blast
    -Evalue(float): Evalue threshold for blast
    -Identity (float): Percent identity for blast
    -Strand (str): Search forward only/reverse/both?

    Returns:
    -Nothing
    """
    subprocess.call(f'blastn -query neORFs_Nuc_all.fa -subject {PathToTE} -perc_identity {Identity} -qcov_hsp_perc {Coverage} -strand {Strand} -evalue {Evalue} -outfmt "6 qseqid sseqid evalue" -out Blast_out', shell = True) 

def get_TE_neORFs(ORF_list:list)->tuple[list,list]:
    """
    Get the hits and no hits from the blast output

    Parameters: 
    -ORF_list: List with all ORFs

    Returns:
    -tuple[list,list]: Two lists (one with hits one with no hits)
    """
    HitList = []
    BlastOut = open("Blast_out","r")
    for line in BlastOut:
        l = line.strip().split("\t")
        HitList.append(l[0])
    BlastOut.close()
    NoHit = set(ORF_list)- set(HitList)
    PercentageExcluded = round((len(set(HitList))/len(set(ORF_list)))*100)
    print(f"Percentage excluded due to Sequence overlap:{PercentageExcluded} %\nNumber kept:{len(NoHit)}\nNumber removed:{len(set(HitList))}")
    return list(NoHit), list(set(HitList))

def run_blast_operation(PathToDESwoMAN:str, PathToTE:str, Coverage:float, Evalue:float, Identity:float, Strand:str, SpeciesList:str)->tuple[list,list]:
    """
    Function for the complete blast workflow. 

    Parameters: 
    -PathToDESwoMAN (str): Path to the DESwoMAN folder
    -PathToTE (str): Path to the TE database 
    -Coverage (float): Qcovhsp for blast
    -Evalue(float): Evalue threshold for blast
    -Identity (float): Percent identity for blast
    -Strand (str): Search forward only/reverse/both?
    -SpeciesList(str): Path to the Species File


    Returns:
    -tuple[list,list]: Two lists (one with hits one with no hits)
    """
    ORF_list =  generateBlastInput(PathToDESwoMAN, SpeciesList)
    run_blast_te(PathToTE, Coverage, Evalue, Identity, Strand)
    NoHit, Hit = get_TE_neORFs(ORF_list)
    os.remove("Blast_out")
    os.remove("neORFs_Nuc_all.fa")
    return NoHit, Hit

def compare2lists(list1:list, list2:list)->bool:
    """
    Compare two lists. If any item from list1 is in list2 it returns True.

    Parameters: 
    -list1 (list)
    -list2(list)

    Returns:
    -bool: True or False
    """
    for item in list1:
        if item in list2:
            return True
    return False

def get_save_neORF(CodingDict:dict, NoncodingDict:dict, TEHitList:list, Tree:str)->list:
    """
    Compare two lists. If any item from list1 is in list2 it returns True.

    Parameters: 
    -CodingDict(dict): Coding neORFs
    -NoncodingDict(dict): Noncoding homologs
    -TEHitList(list): List of neORFs with a TE match
    -Tree(str): Path to the newick folder

    Returns:
    -list: A list of validated homologs
    """

    Validated, ToRemove, ValidSpecies, AllSpecies = [], [], [], []

    for neORF in CodingDict:

        AllSpecies.append(neORF.split("_")[-1][:4])

        if neORF not in TEHitList:
            ToRemove.append(neORF)
            continue

        Coding_Homologs,  NoncodingHomologs = CodingDict[neORF], list(set(NoncodingDict[neORF])) #Initiate the species with a coding/noncoding homolog
        NodeCoding = find_common_node(Tree, Coding_Homologs) #Check the node in the tree
        List_Accepted_Homologs = get_outgroup_from_clade_name(Tree, NodeCoding)#Get the accepted homologs for that node
        index = 0

        #Modify the result so its the species names and NOT the population names
        for elem in List_Accepted_Homologs:
            List_Accepted_Homologs[index] = elem[:4]
            index +=1
        

        if compare2lists(List_Accepted_Homologs, NoncodingHomologs): #check if any of them is present
            Validated.append(neORF)
            ValidSpecies.append(neORF.split("_")[-1][:4])
            #print(neORF, Coding_Homologs, NodeCoding, NoncodingHomologs, compare2lists(List_Accepted_Homologs, NoncodingHomologs))
        else:
            ToRemove.append(neORF)

    PercentageValid = round((len(Validated)/(len(Validated)+ len(ToRemove)))*100, 2)
    print(f"Percentage of valid neORF:{PercentageValid} %\nNumber of safe validated neORFs:{len(Validated)}\nNumber of potentially non de novo neORF:{len(ToRemove)}")
    
    return Validated


def filter_neORFs(Orthopath:str, PathToDESwoMAN:str, PathToTE:str, Coverage:float, Evalue:float, Identity:float, Strand:str, PathToTr:str, CoverageTr:float, EvalueTr:float, IdentityTr:float,Te_check:bool, rna_check:bool, Tree:str, SpeciesList:str, mutations:list, frame:float)->list:
    """
    Run the whole filtering analysis:

    Parameters: 
    -Orthopath (str): Path to the Orthofinder outputfile
    -PathToDESwoMAN (str): Path to the DESwoMAN output 
    -PathToTE (str): Path to the TE database 
    -Coverage (float): Coverage for TE blast 
    -Evalue (float): Evalue for TE blast
    -Identity (float): Percent identity for TE blast
    -Strand (str): Strand for TE blast
    -PathToTr (str): Path to transcript database (e.g. ncrna)
    -CoverageTr(float): Coverage for nuc blast
    -EvalueTr(float): Evalue for nuc blast
    IdentityTr(float): Percent identity for nuc blast
    -Te_check (bool): Filter TE
    -rna_check (bool): Filter nucleotides e.g. RNA
    -Tree (str): Path to  the newick file
    -SpeciesList (str): File with Species info
    -mutations (list): list with accepted mutations
    -frame (float): Maximum percentage of the sequence not affected by frameshift to still count as noncoding

    Returns:
    -list: A list of validated homologs
    """
    print("#####################\nFiltering the neORFs detected with DESwoMAN.\n#####################")
    CodingHomologs = read_orthogroups_info(Orthopath, PathToDESwoMAN, SpeciesList)
    SpeciesWithHomologs, CodingHomologs = get_acceptable_noncoding_homologs(PathToDESwoMAN, CodingHomologs,SpeciesList, mutations, frame)
    if Te_check:
        print("#####################\nStarting the TE filter:")
        NoHit1, Hit1 = run_blast_operation(PathToDESwoMAN, PathToTE, Coverage, Evalue, Identity, "both", SpeciesList)
    else:
        print("#####################\nNo TE search specified... Skipping TE Blast...\n#####################")
        NoHit1, Hit1 = [], []
    if rna_check:
        print("#####################\nStarting the RNA filter:")
        print("Hi")
        NoHit2, Hit2 = run_blast_operation(PathToDESwoMAN, PathToTr, CoverageTr, EvalueTr, IdentityTr, Strand, SpeciesList)
        print("#####################") 
    else:
        print("#####################\nNo RNA search specified... Skipping Nucleotide Blast...\n#####################")
        NoHit2, Hit2 = [], []
    NoHit = list(set(NoHit1 + NoHit2) -  set(Hit1 + Hit2)) 

    if not Te_check and not rna_check:
        x = get_populations(SpeciesList)
        for line in x: 
            Prots =  SeqIO.to_dict(SeqIO.parse(f"{PathToDESwoMAN}/{line}/denovo_protein.fa", "fasta"))
            for key in Prots:
                NoHit.append(key +  "_" + line)
    Validated = get_save_neORF(CodingHomologs, SpeciesWithHomologs, NoHit, Tree)
    return Validated

def create_output(PathToDESwoMAN:str, Valid_neORFs:list, Outpath:str, SpeciesList:str):
    """
    Redo all the DESwoMAN output files but now only with the accepted neORFs

    Parameters: 
    -PathToDESwoMAN (str): Path to the DESwoMAN output 
    - Valid_neORFs (list): List of accepted de novo ORFs
    -Outpath: Name of the output file
    -SpeciesList (str): File with Species info

    Returns:
    -Nothing
    """
    Valid_Transcripts = []
    for i in Valid_neORFs:
        Valid_Transcripts.append(i.split("_")[0]+ "_" + i.split("_")[-1])

    subprocess.call(f"mkdir {Outpath}", shell = True)

    x = get_populations(SpeciesList)
    print("#####################\nStarting to filter the DESwoMAN output\n#####################")
    for line in x:
        
        print(f"Filtering for {line}")

        ProtFile = SeqIO.to_dict(SeqIO.parse(f"{PathToDESwoMAN}/{line}/denovo_protein.fa", "fasta"))
        NuclFile = SeqIO.to_dict(SeqIO.parse(f"{PathToDESwoMAN}/{line}/denovo_nucl.fa", "fasta"))
        TranscriptFile = SeqIO.to_dict(SeqIO.parse(f"{PathToDESwoMAN}/{line}/denovo_transcripts.fa", "fasta"))
        IntronFile = SeqIO.to_dict(SeqIO.parse(f"{PathToDESwoMAN}/{line}/denovo_unspliced_lowered_introns.fa", "fasta"))
        InfoFile = open(f"{PathToDESwoMAN}/{line}/information_file.txt", "r")
        HomologFile= open(f"{PathToDESwoMAN}/{line}/Table_output_step3.txt", "r") 

        #make a new directory for the filtered files
        subprocess.call(f"mkdir {Outpath}/{line}/", shell = True)

        #Redo the protein file
        NewProt = open(f"{Outpath}/{line}/denovo_protein.fa", "w")
        for key, value in ProtFile.items():
            if key + "_" + line in Valid_neORFs:
                NewProt.write(f">{key}\n{value.seq}\n")
        NewProt.close()

        #Redo the nucleotide file
        NewNuc = open(f"{Outpath}/{line}/denovo_nucl.fa", "w")
        for key, value in NuclFile.items():
            if key + "_" + line in Valid_neORFs:
                NewNuc.write(f">{key}\n{value.seq}\n")
        NewNuc.close()

        #Redo the transcript file
        NewTranscript = open(f"{Outpath}/{line}/denovo_transcripts.fa", "w")
        for key, value in TranscriptFile.items():
            if key + "_" + line in Valid_Transcripts:
                NewTranscript.write(f">{key}\n{value.seq}\n")
        NewTranscript.close()

        #Redo the intron file
        NewInt = open(f"{Outpath}/{line}/denovo_unspliced_lowered_introns.fa", "w")
        for key, value in IntronFile.items():
            if key + "_" + line in Valid_neORFs:
                NewInt.write(f">{key}\n{value.seq}\n")
        NewInt.close()

        #Redo the InfoFile:
        NewInfo = open(f"{Outpath}/{line}/information_file.txt", "w")
        for li in InfoFile:
            l = li.split(",")
            if l[5] + "_" + line in Valid_Transcripts:
                NewInfo.write(li)
            elif li[0] == "g":
                NewInfo.write(li)
        NewInfo.close()

        #Redo the HomologFile: 
        NewHomolog = open(f"{Outpath}/{line}/Table_output_step3.txt", "w")
        for li in HomologFile:
            l = li.split(",")
            if l[0] + "_" + line in Valid_neORFs:
                NewHomolog.write(li)
            elif li[0] == "d":
                NewHomolog.write(li)
        NewHomolog.close()


#Give the option to choose which enabling mutations to include and what threshold
        
if __name__ == "__main__":

    #Initialize all arguments
    parser = argparse.ArgumentParser(description="Get high confidence neORFs from the DESwoMAN output",epilog="-------------------------")
    parser.add_argument("--te_db", help="Path to theTE database", type=str)
    parser.add_argument("--te_cov", help="Coverage Blast against TE (Default = 80)", type=float, default= 80)
    parser.add_argument("--te_idt", help="Percent Identity Blast against TE (Default = 80)", type=float, default= 80)
    parser.add_argument("--te_eval", help="E value threshold Blast against TE (Default = 0.001)", type=float, default= 0.001)
    parser.add_argument("--ortho", help="Path to the Orthogroups.txt file from OrthoFinder", type=str)
    parser.add_argument("--deswoman", help="Path to the DESwoMAN Outputs folder (Default = working directory)", type=str, default = "")
    parser.add_argument("--rna_check", help="Blast against a rna database",action="store_true")
    parser.add_argument("--te_check", help="Blast against a TE database",action="store_true")
    parser.add_argument("--tr_db", help="Path to the rna database", type=str)
    parser.add_argument("--tr_cov", help="Coverage Blast against rna(Default = 80)", type=float, default= 80)
    parser.add_argument("--tr_idt", help="Percent Identity Blast against rna(Default = 80)", type=float, default= 80)
    parser.add_argument("--tr_eval", help="E value threshold Blast against rna (Default = 0.001)", type=float, default= 0.001)
    parser.add_argument("--tr_strand", help="Strand to consider in  Blast against rna (Default = plus)", type=str, default="plus")
    parser.add_argument("--out", help="Name of the output folder (Default: FilteredNeORFs)", type=str, default="FilteredNeORFs")
    parser.add_argument("--tree", help="Newick tree file; needs all internal nodes", type=str)
    parser.add_argument("--species_file", help="Text file with all species", type=str)
    #Add mutations as argument
    parser.add_argument("--accepted_mutations", help="List all mutations to be accepted as noncoding", type=str, default="complete,start,stop,premature-stop,frameshift")
    parser.add_argument("--frameshift_score", help="Frameshift_score to be accepted as enough to deem a sequence noncoding.", type=float, default=50)
    

    print("-------------------------\nFilter neORFs\nV.1.0\nAuthor:Marie Lebherz\n-------------------------\n")

    #Read the input    
    args = parser.parse_args()
    TEPath =args.te_db
    DESwoMANPath = args.deswoman
    OrthoPath = args.ortho
    Coverage = args.te_cov
    Evalue = args.te_idt
    Identity = args.te_eval
    CoverageTr = args.tr_cov
    EvalueTr = args.tr_idt
    IdentityTr = args.tr_eval
    PathToTr = args.tr_db
    Strand = args.tr_strand
    DoTranscript = args.rna_check
    DoTE = args.te_check
    Outpath = args.out
    Tree = args.tree
    SpeciesList = args.species_file

    #mutations
    muts = args.accepted_mutations
    mutations = muts.split(",")
    frame = args.frameshift_score
    
    #Check if all is correct
    if not DESwoMANPath:
        print("No argument supplied. Exiting...")
        sys.exit()
    if os.path.isdir(DESwoMANPath) != True:
        print("One of your input directories is not a directory. Exiting...")
        sys.exit()
    if not OrthoPath:
        OrthoPath = "NotSupplied"
        print("You did not supply an Orthogroup File. All neORFs are treated as if they have no detected coding Homologs.")
    #Start running the analysis
    print("Starting the analysis!\n\n")
    Valid_neORFs = filter_neORFs(OrthoPath, DESwoMANPath, TEPath, Coverage, Evalue, Identity, Strand, PathToTr, CoverageTr, EvalueTr, IdentityTr,DoTE, DoTranscript, Tree, SpeciesList, mutations, frame)     
    create_output(DESwoMANPath, Valid_neORFs, Outpath, SpeciesList)
    print("\nFinished!")
    print("Goodybe ;)")

#Mutations available:
#-Start: present or absent
#-Stop: present or absent
#- Frameshift > Score (num)
#- Premature Stop: present or absent
#- Transcription: present, incomplete, antisense or absent