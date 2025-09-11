#!/usr/bin/python3
import argparse
import sys
import os
def get_populations(File:str)-> list:
    """
    Read in a text file with one population/id on each line

    Parameters: 
    - File (Path to a text file containing all query (!) species, one per line)

    Returns:
    list of species
    """
    Strains = []
    Info = open(File, "r")
    for line in Info:
        l = line.strip()
        Strains.append(l)
    return Strains


def read_orthogroups_info2(OrthoPath:str)->dict:
    """
    Reads in the orthogroups info from orthofinder (only used to filter the neORF homologs)

    Parameters: 
    - OrthoPath (str): Path to the Orthofinder Orthogroups.txt output file
    -Important (!): The species/line name needs to be included here (e.g. STRG1.1_1_20_305_DmelFI)

    Returns:
    - dict: A dictionary where the neORF id is the key and the corresponding populations in its orthogroup are the value 
    """
    Orthogroups = {}
    Info = open(OrthoPath, "r")
    for line in Info:
        l = line.strip()
        OG_ID = l.split(":")[0]
        OGs = l.split(":")[1][1:].split(" ")
        Pops = []
        for i in OGs:
            if i:
                Pops.append(i.split("_")[-1])
        Orthogroups[OG_ID] = [OGs, list(set(Pops))]  #Fill the dictionary with the populations in an orthogroup
    print("Done")
    Final = {}
    for key, value in Orthogroups.items(): #Now do it per ORF
        Seqs = value[0]
        for i in Seqs:
            Final[i] = value[1]
    return Final



def filter_output(OrthoPath:str, PathDESwoMAN:str, File:str):
    """
    Filters out coding homologs from the DESwoMAN output. Returns 2 files:
    - Table_output_step3_filteredbyOGs.txt (A file where the coding homologs are removed, use this one for all follow up analysis)
    - Table_output_step3_RemovedCoding.txt (A file with the removed coding homologs)

    Parameters: 
    - OrthoPath (str): Path to the Orthofinder Orthogroups.txt output file
    - PathDESwoMAN (str): Path to the DESwoMAN folder 

    Returns:
    - nothing
    """
    x = get_populations(File)
    print("Start reading in the orthogroup info!")
    OrthogroupInfo = read_orthogroups_info2(OrthoPath)
    print("Finished reading in the orthogroup info!\nStart filtering the output files.")
    for line in x:
        print("Start the filtering for :", line)
        lines_to_keep = []
        lines_exclude = []
        path = PathDESwoMAN + line + "/Table_output_step3.txt"
        Outfile = open(path, "r")
        for li in Outfile:
            l = li.strip().split(",")
            strain = l[0] + "_" + line
            target = l[1]
            if l[0] == "denovo":
                continue
            if target in OrthogroupInfo[strain]:
                lines_exclude.append(li)
            else:
                lines_to_keep.append(li)
        Out = open(PathDESwoMAN + line + "/Table_output_step3_InvalidHomologs.txt", "w") #Makes the file to use later
        Out.write("denovo,target_genome,validated_hit,genomic_position_homolog,start,stop,Indels,perc_seq_not_affected_by_frameshift,substitutions,premature_stop,pos_premature_stop,transcription_status,intron_lowered_homolog,homolog_in_ali,denovo_in_ali\n")
        for i in lines_exclude:
            Out.write(i)
        Out.close()
        Out2 = open(PathDESwoMAN + line + "/Table_output_step3_ValidatedNoncodingHomologs.txt" , "w") #Makes a file with the coding homologs (for control etc.)
        Out2.write("denovo,target_genome,validated_hit,genomic_position_homolog,start,stop,Indels,perc_seq_not_affected_by_frameshift,substitutions,premature_stop,pos_premature_stop,transcription_status,intron_lowered_homolog,homolog_in_ali,denovo_in_ali\n")
        for i in lines_to_keep:
            Out2.write(i)
        Out2.close()
        print("Finished")
        print("Number kept (= no coding homolog in orthogroup with the neORF):", len(lines_to_keep))
        print("Number removed (= coding homologs in an orthogroup with the neORF):", len(lines_exclude), "\n")        
        


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Filter coding homologs from the final DESwoMAN output using Orthofinder Orthogroups containing both query and outgroup expressed ORFs.",epilog="-------------------------")
    parser.add_argument("--ortho", help="Path to the Orthogroups.txt file from OrthoFinder (e.g. /home/usr/orthofinder/orthogroups.txt)", type=str)
    parser.add_argument("--deswoman", help="Path to the DESwoMAN Outputs folder (e.g. /home/usr/deswoman/, Default = working directory)", type=str, default = "")
    parser.add_argument("--query_file", help="Path to a text file where each query species is on one line.", type=str)

    print("-------------------------\nRemove noncoding homologs that have an expressed homolog via Orthofinder\nV.1.0\nAuthor:Marie Lebherz\n-------------------------\n")

    args = parser.parse_args()

    DESwoMANPath = args.deswoman
    OrthoPath = args.ortho
    File = args.query_file

    if not DESwoMANPath or not OrthoPath or not DESwoMANPath:
        print("Filepath missing! Exiting...")
        sys.exit()

    print("Starting the analysis!")  
    filter_output(OrthoPath, DESwoMANPath,File)  

    print("Finished!")
    print("Goodybe ;)")
