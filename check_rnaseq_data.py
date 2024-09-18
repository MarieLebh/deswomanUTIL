#!/usr/bin/python3

"""
Check your RNAseq experiment with RSeQC.

Prequisite: Have RSeQC and all dependencies installed on your workstation. 
For more information on that: https://rseqc.sourceforge.net/
"""
import subprocess

def get_populations():
   #read in population names from a file and output them in a list
    Strains = []
    Info = open("File.txt", "r")
    for line in Info:
        l = line.strip()
        Strains.append(l)
    return Strains
  
def gff2bed(gff, Strain):
    #Transform GeMoMa gff to very basic bedfile
    File = open(gff, "r")
    #Bed = open(gff.split(".")[0] + ".bed", "w")
    Bed = open(Strain + ".bed", "w")
    for line in File:
        if line[0] == "#":
            continue
        else:
            l = line.strip().split("\t")
            Chrom = l[0]
            Source = l[1]
            Type =l[2]
            Start =int(l[3]) -1
            End =l[4] 
            Score =l[5]
            Strand =l[6]
            Phase =l[7]
            Name =l[8].split(";")[0].split("=")[1]

            Bed.write(Chrom + "\t" + str(Start) + "\t" + End + "\t" + Type + "_" + Name + "\t" + Score + "\t" + Strand + "\n")
    Bed.close()

def check_strandedness(Strain):
    #Check the strandedness of one strain with RSeQCs "infer_experiment" function
    Path_Bam = "/global/students/homes/m_lebh01/Droso_TranscriptomeAssemblies/" + Strain + "/S2_Hisat2/hisat2/accepted_hits.sorted.bam"
    Path_GFF3 = "/global/students/homes/m_lebh01/Droso_GenomeAnnotations/Genome_annotations/" + Strain + ".gff"
    print("Starting to check the strandedness of: " + Strain +"\n")
    print("This is the genome annotation file used: ", Path_GFF3)
    print("This is the Bamfile used: ", Path_Bam)
    print("Start running RSeQC! \n")
    gff2bed(Path_GFF3, Strain)
    #Run program
    subprocess.call("infer_experiment.py -r "+ Strain + ".bed  -i " + Path_Bam, shell = True)
    subprocess.call("rm " + Strain + ".bed", shell = True)


def check_strandedness_all():
    #Run strandedness analysis for all populations spefified in a file
    print("Starting to use RSeQC to check the strandedness of the RNAseq data! \n For more information look at the documentation at: 'https://rseqc.sourceforge.net/'.")
    pops = get_populations()
    for pop in pops:
        if "Dsim" in pop or "Dsuz" in pop:
            print("\n\n\nDataset: RNAseq dataset 1 (Dsim/Dsuz)\n")
        else:
            print("\n\n\nDataset: RNAseq dataset 2 (Dana/Dimm/Dmel/Dsub)\n")            
        check_strandedness(pop)
 
