#!/usr/bin/python3
import sys
def gff2bed():
    #Transform gff/gtf files to bedfiles
    #Usage python3 gff2bed.py path_to_gff outfilename
    #Example python3 gff2bed.py GenomeAnnotation.gff GenomeAnnotation
  
    gff = sys.argv[1] #input gff path
    outname = sys.argv[2] #input Outfile name 

    File = open(gff, "r")

    Bed = open(outname + ".bed", "w")
    for line in File:
        if line[0] == "#": #ignore hashtags
            continue
        else: #Get all the informaton from the different columns
            l = line.strip().split("\t") 
            Chrom = l[0]
            Source = l[1]
            Type = l[2]
            Start = int(l[3]) -1 #Adapt the coordinates (gff is 1 based, bed is 0 based)
            End = l[4] 
            Score = l[5]
            Strand = l[6]
            Phase = l[7]
            Name = l[8].split(";")[0].split("=")[1]

            Bed.write(Chrom + "\t" + str(Start) + "\t" + End + "\t" + Type + "_" + Name + "\t" + Score + "\t" + Strand + "\n")
    Bed.close()


if __name__ == "__main__":
  gff2bed()
