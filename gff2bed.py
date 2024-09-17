#!/usr/bin/python3

def gff2bed(gff, Strain):
  #Very basic way to transform a gff file into a bedfile 
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

