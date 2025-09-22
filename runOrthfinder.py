import subprocess
import sys
import argparse
from Bio import SeqIO

"""
Prepare and run Orthofinder
Author: Marie
Date: 16.05.25
"""

def get_populations(SpeciesFile:str)->list:
    """
    Read in a text file with one population/id on each line

    Parameters: 
    - SpeciesFile (str): Path to the file with the sample ids

    Returns:
    -list: list of lines
    """
    Strains = []
    Info = open(SpeciesFile, "r")
    for line in Info:
        l = line.strip()
        Strains.append(l)
    return Strains

def save_data(SpeciesFile:str, DESwoMAN:str):
    """
    Read in a text file with one population/id on each line

    Parameters: 
    - SpeciesFile (str): Path to the file with the sample ids
    - DESwoMAN (str): Path to the DESwoMAN folder

    Returns:
    - Nothing (but creates Orthofinder input files)
    """
    x = get_populations(SpeciesFile)
    subprocess.run(["mkdir", "Orthofinder"], check = True)
    for line in x:
        Out_dict = {}
        path = f"{DESwoMAN}{line}/denovo_protein.fa"
        Nuc_dict = SeqIO.to_dict(SeqIO.parse(path, "fasta"))
        for key, value in Nuc_dict.items():
            key = key + "_" + line
            value = value
            Out_dict[key] = value
        Out = open("Orthofinder/" + line + "_DeNovoProt.fa", "w")
        for key, value in Out_dict.items():
            Out.write(">" + key + "\n")
            Out.write(str(value.seq) + "\n")
        Out.close()

def run_Orthofinder(num_threads:int):
    """
    Run Orthofinder on the created folder

    Parameter:
    - num_threads(int): Number of threads

    Returns:
    - Nothing (but runs Orthofinder)
    """
    print("Starting to run Orthofinder:")
    sys.stdout.flush()
    cmd =["orthofinder",
          "-f", "Orthofinder",
          "-t", str(num_threads)]
    subprocess.run(cmd, check = True)
    print("Finished running Orthofinder.")
    sys.stdout.flush()

def main():
    """
    Main function

    Parameter:
    - none

    Returns:
    - Nothing
    """
    parser = argparse.ArgumentParser(description="Run Orthofinder using neORF protein files",epilog="-------------------------")
    parser.add_argument("--threads", help="Number of threads for Orthofinder to use", type=int)
    parser.add_argument("--create_folder", help="Create the folder an exit",action="store_true")
    parser.add_argument("--species_file", help="Text file with all species", type=str)
    parser.add_argument("--deswoman", help="Path to the DESwoMAN Outputs folder (Default = working directory)", type=str, default = "")

    print("-------------------------\nGet NeORF Orthogroups\nV.1.0\nAuthor:Marie Lebherz\n-------------------------\n")
    args = parser.parse_args()
    num_threads = args.threads
    folder = args.create_folder
    desw = args.deswoman
    specFile = args.species_file

    print("Starting the analysis!")

    if folder:
        save_data(desw,specFile)
    else:
        save_data(desw,specFile)
        run_Orthofinder(num_threads)

    print("Finished!")
    print("Goodybe :)")
    
if __name__ == "__main__":
    main()
