#!/usr/bin/python3
"""
Calculate TE overlap between a TE and a sequence bedfile.

At the moment this works for TE annotation from either the FasTE pipeline/Repeat Masker or the Transposon Ultimate pipeline.
"""
from pybedtools import BedTool # type: ignore

def merge_intervals(intervals):
  #Merge overlapping intervals to non-overlapping intervals
    if not intervals:
        return []

    intervals.sort(key=lambda x: x[0])  # Sort intervals by start time

    merged = []
    for interval in intervals:
        if not merged or interval[0] > merged[-1][1]:
            merged.append(interval)
        else:
            merged[-1][1] = max(merged[-1][1], interval[1])

    return merged

def overlap_intervals(intervals1, intervals2):
    overlapping_intervals = []
    intervals1 = [intervals1[0]]
    for i in intervals1:
        for j in intervals2:
            if (int(i[0]) <= int(j[1]) and int(j[0]) <= int(i[1])):
                # Merge overlapping intervals
                start = max(int(i[0]), int(j[0]))
                end = min(int(i[1]), int(j[1]))
                overlapping_intervals.append([start, end])
    return overlapping_intervals

#Get the TE overlap between 2 sequences

def get_TE_overlap(Seq_Path, TE_Path, Outfile, Program):
    #get the TE overlap between a TE annotation bedfile and a sequence file (e.g. genes, transcripts, ...)
  
    #Specify all necessary dictionaries
    IntervalsTE = {}
    IntervalsSeq = {}
    Intervals = {}
    List_all = []
    Control_dict = {}
    Len_dict = {}
    
    TE_type_dict = {}
    
    #Read in Bedfile of the sequence file
    Sequence = BedTool(Seq_Path)
    Sequence = Sequence.sort()
    
    #Read in Bedfile of the TE file
    TE = BedTool(TE_Path)
    TE = TE.sort()
    
    #Now run the bedtools intersection
    Intersection = Sequence.intersect(TE, wao=True) 
    for gene in Intersection:  

        if "CDS" in gene[3]:
             continue

        Len_dict[gene[3]] = int(gene[2]) - int(gene[1]) 
        if gene[3] in Control_dict.keys():
            Control_dict[gene[3]] += [[gene[1], gene[2], gene[5], gene[7], gene[8], gene[11]]]
        else:
             Control_dict[gene[3]] = [[gene[1], gene[2], gene[5], gene[7], gene[8], gene[11]]]

        List_all.append(gene[3]) #append all gene names

        if int(gene[-1]) != 0: 

            #Check the strand
            if gene[5] == "." or gene[5] == "+":
                    strand = "+"
            else:
                    strand = "-"

            #Save the overlapping intervals
            if gene[3] in IntervalsTE and strand == gene[11]:
                IntervalsTE[gene[3]] += [[gene[7], gene[8]]]
                IntervalsSeq[gene[3]] += [[gene[1], gene[2]]]

                if Program == "TransposonUltimate": #Save the TE types in a dictionary depending on the specified program
                    TE_type_dict[gene[3]] += [[gene[9].split("(")[1][:-1], gene[12]]]   
                elif Program == "FasTE":
                    TE_type_dict[gene[3]] += [[gene[12].replace("__", ""), gene[16]]]

            elif gene[3] not in Intervals and strand == gene[11]:
                IntervalsTE[gene[3]] = [[gene[7], gene[8]]]
                IntervalsSeq[gene[3]] = [[gene[1], gene[2]]]

                if Program == "TransposonUltimate": #Save the TE types in a dictionary depending on the specified program
                    TE_type_dict[gene[3]] = [[gene[9].split("(")[1][:-1], gene[12]]]   
                elif Program == "FasTE": 
                    TE_type_dict[gene[3]] = [[gene[12].replace("__", ""), gene[16]]]

    #Extract the intervals where the TE actually overlaps the gene sequence
    for gene, intervals in IntervalsTE.items():
         Intervals[gene] = overlap_intervals(IntervalsSeq[gene], intervals)

    #Now get the actual not overlaping intervals
    Final_dict = {}
    for gene, intervals in Intervals.items():
        new_intervals = merge_intervals(intervals)
        sum_overlap = 0
        Number_TE = len(IntervalsTE[gene])

        for i in new_intervals:
            sum_overlap += int(i[1])-int(i[0])

        Final_dict[gene] = [sum_overlap, sum_overlap/Len_dict[gene], Number_TE]

    #Add the ids without TE overlap back to the dictionary
    for i in List_all:
         if i not in Final_dict.keys():
              Final_dict[i] = [0,0,0]
    
    #Save the results in an output file

    #First TE proportion
    Out1 = open(Outfile + "_" + Program +  "TE_overlap_info.tsv", "w")
    Out1.write("Seq_id\tOverlap_bp\tPercentage_overlap\tTE_count\n")
    for key, value in Final_dict.items():
         Out1.write(key + "\t" + str(value[0]) + "\t" + str(value[1]) + "\t" + str(value[2]) + "\n")
    Out1.close()

    #Then TE type
    Out2 = open(Outfile + "_" + Program + "TE_Types.tsv", "w")
    Out2.write("Seq_id\tTE_Type\tTE_overlap\n")
    for key, value in TE_type_dict.items():
         for TE in value:
              Out2.write(key + "\t" + str(TE[0]) + "\t" + str(TE[1]) + "\n")
    Out2.close()


    return Final_dict, Control_dict, TE_type_dict
