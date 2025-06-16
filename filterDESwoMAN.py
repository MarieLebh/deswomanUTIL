from Bio import SeqIO
import os
import subprocess
from collections import Counter
import argparse
import sys
from Bio import Phylo
from io import StringIO

"""
Create a subset of high confidence neORFs
Author: Marie
Date: 16.05.25

This script filters the DESwoMAN output based on the following phylogeny:

CRITERIA FOR FILTERING:

    Less than 80% identity and coverage overlap with a TE (curated Drosophila database)
    At least one noncoding homolog in an outgroup species
    If there is a species with no coding mutation present that species is added to the coding species.
"""

def get_populations(SpeciesList):
    #Read in a text file with one population in each line
    Strains = []
    Info = open(SpeciesList, "r")
    for line in Info:
        l = line.strip()
        Strains.append(l)
    return Strains

def find_common_node(newick, species_list):
    #Input a tree and a list of species and return their most recent common ancestor
    tree = Phylo.read(newick, "newick")
    # Normalize species names (remove dots if needed, optional depending on your data)
    species_list = [s.replace(" ", "").strip() for s in species_list]
    # Find the MRCA
    mrca = tree.common_ancestor(species_list)
    # Return the name of the MRCA node
    return mrca.name

def get_outgroup_from_clade_name(tree, clade_name):
    #Input a clade name (needs to be included in the tree), it will return all terminal branches who are outgroups
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

def read_orthogroups_info(OrthoPath, PathToDESwoMAN, SpeciesList):
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

def get_acceptable_noncoding_homologs(PathToDESwoMAN, CodingDict, SpeciesList):
    #Filter the DESwoMAN output for mutations that are acceptable
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
            if l[2] == "P" and start != "NA": #Only include present homologs. 
                if start != "P" or stop != "P"  or premature_stop != "A" or float(perc_seq_not_affected_by_frameshift) < 60 or (transcription_status != "complete" and transcription_status != "partial"):#transcription_status != "complete": #Change here (!)
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

def generateBlastInput(PathToDESwoMAN, SpeciesList):
    #Merge all fasta files of neORFs together to create 
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

def run_blast_te(PathToTE, Coverage, Evalue, Identity, Strand):
    #Run blast against a database
    subprocess.call(f'blastn -query neORFs_Nuc_all.fa -subject {PathToTE} -perc_identity {Identity} -qcov_hsp_perc {Coverage} -strand {Strand} -evalue {Evalue} -outfmt "6 qseqid sseqid evalue" -out Blast_out', shell = True) 

def get_TE_neORFs(ORF_list):
    #Read  the blast output and give a list with all that have no hit
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

def run_blast_operation(PathToDESwoMAN, PathToTE, Coverage, Evalue, Identity, Strand, SpeciesList):
    #This function merges all blast related tasks
    ORF_list =  generateBlastInput(PathToDESwoMAN, SpeciesList)
    run_blast_te(PathToTE, Coverage, Evalue, Identity, Strand)
    NoHit, Hit = get_TE_neORFs(ORF_list)
    os.remove("Blast_out")
    os.remove("neORFs_Nuc_all.fa")
    return NoHit, Hit

def compare2lists(list1, list2):
    #Comoare two lists. If any item from list1 is in list2 it returns True.
    for item in list1:
        if item in list2:
            return True
    return False

def get_save_neORF(CodingDict, NoncodingDict, TEHitList, Tree):
    #Filter ORFs based on whether they have a noncoding homolog in synteny and have no hit in a Drosophlia TE database.
    
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
    
    #Count the percentage of neORFs retained in each species
    #DmelPerc = (ValidSpecies.count("Dmel")/AllSpecies.count("Dmel"))*100
    #DsimPerc = (ValidSpecies.count("Dsim")/AllSpecies.count("Dsim"))*100
    #DsuzPerc = (ValidSpecies.count("Dsuz")/AllSpecies.count("Dsuz"))*100
    #DanaPerc = (ValidSpecies.count("Dana")/AllSpecies.count("Dana"))*100
    #DsubPerc = (ValidSpecies.count("Dsub")/AllSpecies.count("Dsub"))*100
    #DimmPerc = (ValidSpecies.count("Dimm")/AllSpecies.count("Dimm"))*100
    
    #print(f"#####################\nPercentages of neORF retained after filtering (per Species):\nDmel:{round(DmelPerc,2)} %\nDsim:{round(DsimPerc,2)} %\nDsuz:{round(DsuzPerc,2)} %\nDana:{round(DanaPerc,2)} %\nDsub:{round(DsubPerc,2)} %\nDimm:{round(DimmPerc,2)} %")

    return Validated


def filter_neORFs(Orthopath, PathToDESwoMAN, PathToTE, Coverage, Evalue, Identity, Strand, PathToTr, CoverageTr, EvalueTr, IdentityTr,Te_check, rna_check, Tree, SpeciesList):
    #This runs the whole filtering analysis
    print("#####################\nFiltering the neORFs detected with DESwoMAN.\n#####################")
    CodingHomologs = read_orthogroups_info(Orthopath, PathToDESwoMAN, SpeciesList)
    SpeciesWithHomologs, CodingHomologs = get_acceptable_noncoding_homologs(PathToDESwoMAN, CodingHomologs,SpeciesList)
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

def create_output(PathToDESwoMAN, Valid_neORFs, Outpath, SpeciesList):
    #Redo all the DESwoMAN output files but now only with the accepted neORFs

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
    Valid_neORFs = filter_neORFs(OrthoPath, DESwoMANPath, TEPath, Coverage, Evalue, Identity, Strand, PathToTr, CoverageTr, EvalueTr, IdentityTr,DoTE, DoTranscript, Tree, SpeciesList)     
    create_output(DESwoMANPath, Valid_neORFs, Outpath, SpeciesList)
    print("\nFinished!")
    print("Goodybe ;)")

 #python3 GetHighConfidenceNeORFs.py --te_db /global/scratch2/m_lebh01/TranscriptModelling/ReferenceNCBI/TE_DatabaseFlybase/FyBase_TE.fa --ortho Orthofinder/OrthoFinder/Results_May19/Orthogroups/Orthogroups.txt --deswoman ./
