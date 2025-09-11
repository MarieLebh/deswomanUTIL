#!/usr/bin/python3

import argparse
import sys
import os

"""
GENERATE AN OUTPUT GFF WITH ALL DESwoMAN NeORFs

usage: getGFFfromOut.py [-h] [--deswoman DESWOMAN] [--gtf GTF] [--outname OUTNAME]

For each neORF this will extract:
(1) Exon positions (Can simply take from gtf) 
(2) CDS positions
(3) UTRs (3prime/5prime)
(4) Start codon
(5) Stop codon
(6) Transcript position 
(7) NeORF position 
"""

def extract_cds_bed(PathToDESwoMAN:str ,PathToTranscriptome:str)->dict:
    """
    Extract the CDS of neORFs into a bedfile 
 
    Parameters: 
    - PathToDESwoMAN (str): Path to the DESwoMAN info file
    - PathToTranscriptome (str): Path to the transcriptome gtf file

    Returns:
    - dict: dictionary with the CDS (stored as gff lines)
    """
    GTF = open(PathToTranscriptome , "r")
    Exons = {}
    for lin in GTF:
        if lin[0] == "#":
            continue  #SkiptStringTieHeader
        l = lin.strip().split("\t")
        if l[2] == "exon":
            Tr_ID = l[8].split(";")[1].split(" ")[2]
            Tr_ID = Tr_ID[1:len(Tr_ID)-1] 
            if Tr_ID not in Exons:
                Exons[Tr_ID] = [[[l[3], l[4]]], l[6], l[0]]
            else:
                Exons[Tr_ID] [0].append([l[3], l[4]])
    exon_dict = Exons
    GTF.close()
    
    InfoFile =  open(PathToDESwoMAN , "r")
    NeORFs = {}
    for li in InfoFile:        
        l = li.strip().split(",")
        if l[0] != "gene_name":
            NeORFs[l[9]] = [l[5],int(l[9].split("_")[-2]),int(l[9].split("_")[-1])+3] 
    InfoFile.close()
    orf_dict = NeORFs

    FinalDict = {}
    for orfname, (transcript, orf_start, orf_end) in orf_dict.items():
        exon_list, strand, chrom = exon_dict[transcript]
        # Build list of exons in order of transcription (not genome!)
        if strand == "-":
            exons = exon_list[::-1]  # reverse, so we're always moving 5'->3'
        else:
            exons = exon_list
        transcript_pointer = 1  
        for ex_idx, (exon_start, exon_end) in enumerate(exons):
            exon_start, exon_end = int(exon_start), int(exon_end)
            exon_len = exon_end - exon_start + 1

            exon_tr_start = transcript_pointer
            exon_tr_end = transcript_pointer + exon_len - 1

            # Does the ORF overlap this exon?
            overlap_start = max(exon_tr_start, orf_start)
            overlap_end = min(exon_tr_end, orf_end)

            if overlap_start <= overlap_end:  # Overlap exists
                offset_start = overlap_start - exon_tr_start
                offset_end = overlap_end - exon_tr_start

                if strand == "+":
                    cds_genomic_start = exon_start + offset_start
                    cds_genomic_end = exon_start + offset_end
                else:
                    # on minus strand, exons are counted from the high coordinate down
                    cds_genomic_end = exon_end - offset_start
                    cds_genomic_start = exon_end - offset_end

                bed_start = min(cds_genomic_start, cds_genomic_end) - 1  
                bed_end = max(cds_genomic_start, cds_genomic_end)       
                bed_line = [
                    chrom,  
                    "deswomanUTIL",
                    "CDS",             
                    str(bed_start+1),
                    str(bed_end),
                    ".",
                    strand,
                    ".",
                    f'gene_id "{orfname.split("_")[0][:-2]}"; transcript_id "{orfname.split("_")[0]}"; neORF_id "{orfname}"'
                ]

                if orfname.split("_")[0] in FinalDict:
                    FinalDict[orfname.split("_")[0]]+= ["\t".join(bed_line)+"\n"]
                else:
                    FinalDict[orfname.split("_")[0]] = ["\t".join(bed_line)+"\n"]
            transcript_pointer += exon_len
    return FinalDict


def generate_final_file(DESwoMAN:str, GTF:str, Outname:str):
    """
    Generate a final file from the DESwoMAN output
 
    Parameters: 
    - DESwoMAN (str): Path to the DESwoMAN info file
    - GTF (str): Path to the transcriptome gtf file

    Returns:
    - nothing
    """
    DESWOMAN = open(DESwoMAN, "r")
    GTFfile = open(GTF, "r")
    Out = open(Outname + ".gff", "w")
    Out.write(f"# deswomanUTIL -in {DESwoMAN} -gtf {GTF} \n# version 1\n")    
    CDS= extract_cds_bed(DESwoMAN, GTF)

    SaveDict = {}

    DESWOMAN.readline()
    for line in DESWOMAN:
        l = line.split(",")
        SaveDict[l[5]] = [f'{l[1]}\tDESwoMAN\tneORF\t{int(l[-3])+1}\t{int(l[-2])+1}\t.\t{l[2]}\t.\tgene_id "{l[0]}"; transcript_id "{l[5]}"; neORF_id "{l[-4]}";\n']


        if l[2] == "+":
            SaveDict[l[5]] += [f'{l[1]}\tdeswomanUTIL\tstart_codon\t{int(l[-3])+1}\t{int(l[-3])+4}\t.\t{l[2]}\t.\tgene_id "{l[0]}"; transcript_id "{l[5]}"; neORF_id "{l[-4]}";\n']
            SaveDict[l[5]] += [f'{l[1]}\tdeswomanUTIL\tstop_codon\t{int(l[-2])-2}\t{int(l[-2])+1}\t.\t{l[2]}\t.\tgene_id "{l[0]}"; transcript_id "{l[5]}"; neORF_id "{l[-4]}";\n']
            SaveDict[l[5]] += [f'{l[1]}\tdeswomanUTIL\t3_prime_utr\t{int(l[11])+2}\t{int(l[7])}\t.\t{l[2]}\t.\tgene_id "{l[0]}"; transcript_id "{l[5]}"; neORF_id "{l[-4]}";\n']
            SaveDict[l[5]] += [f'{l[1]}\tdeswomanUTIL\t5_prime_utr\t{int(l[6])}\t{int(l[10])}\t.\t{l[2]}\t.\tgene_id "{l[0]}"; transcript_id "{l[5]}"; neORF_id "{l[-4]}";\n']
        else:
            SaveDict[l[5]] += [f'{l[1]}\tdeswomanUTIL\tstop_codon\t{int(l[-3])+1}\t{int(l[-3])+4}\t.\t{l[2]}\t.\tgene_id "{l[0]}"; transcript_id "{l[5]}"; neORF_id "{l[-4]}";\n']
            SaveDict[l[5]] += [f'{l[1]}\tdeswomanUTIL\tstart_codon\t{int(l[-2])-2}\t{int(l[-2])+1}\t.\t{l[2]}\t.\tgene_id "{l[0]}"; transcript_id "{l[5]}"; neORF_id "{l[-4]}";\n']
            SaveDict[l[5]] += [f'{l[1]}\tdeswomanUTIL\t5_prime_utr\t{int(l[11])+2}\t{int(l[7])}\t.\t{l[2]}\t.\tgene_id "{l[0]}"; transcript_id "{l[5]}"; neORF_id "{l[-4]}";\n']
            SaveDict[l[5]] += [f'{l[1]}\tdeswomanUTIL\t3_prime_utr\t{int(l[6])}\t{int(l[10])}\t.\t{l[2]}\t.\tgene_id "{l[0]}"; transcript_id "{l[5]}"; neORF_id "{l[-4]}";\n']
       

    for line in GTFfile:
        if line[0] != "#":
            l = line.split("\t")
            if l[8].split(";")[1].split(" ")[2][1:-1] not in SaveDict:
                continue
            SaveDict[l[8].split(";")[1].split(" ")[2][1:-1]] += [f'{l[0]}\tStringtie\t{l[2]}\t{l[3]}\t{l[4]}\t{l[5]}\t{l[6]}\t.\tgene_id "{l[8].split(";")[0].split(" ")[1][1:-1]}"; transcript_id "{l[8].split(";")[1].split(" ")[2][1:-1]}";\n']

    for key in CDS:
        for value in CDS[key]:
            SaveDict[key] += value

    for key in SaveDict:
        for item in SaveDict[key]:
            Out.write(item)
            
    #Close everything again
    DESWOMAN.close()
    GTFfile.close()
    Out.close()

#Test paths
#DESwoMAN = "/global/students/research/m_lebh01/DESwoMAN/FinalTestDESwoMAN/MarieCorrections/genomes_strat1/DESwoMAN_denovo_output/information_file.txt"
#GTF = "/global/students/research/m_lebh01/DESwoMAN/FinalTestDESwoMAN/MarieCorrections/transcriptomes_strat1/AK5.gtf"
#generate_final_file(DESwoMAN, GTF)


if __name__ == "__main__":
    #Initialize all arguments
    parser = argparse.ArgumentParser(description="Get high confidence neORFs from the DESwoMAN output",epilog="-------------------------")
    parser.add_argument("--deswoman", help="Path to the DESwoMAN info file", type=str)
    parser.add_argument("--gtf", help="Path to the transcriptome gtf file", type=str)
    parser.add_argument("--outname", help="Name of the output file (no file extension)", type=str, default = "DESwoMAN")
    print("-------------------------\nGet a DESwoMAN GFF file\nV.1.0\nAuthor:Marie Lebherz\n-------------------------\n")

    #Read the input    
    args = parser.parse_args()
    GTF =args.gtf
    DESwoMANPath = args.deswoman
    Outname = args.outname

    #Check if all is correct
    if not DESwoMANPath:
        print("No DESwoMAN path supplied. Exiting.")
        sys.exit()
    elif not GTF:
        print("No transcriptome assembly path supplied. Exiting.")
        sys.exit()

    #Start running the analysis
    print("Welcome to deswomanUTIL gff converter...\n")
    generate_final_file(DESwoMANPath, GTF, Outname)
    print("Your file was successfully transformed!!")
    print("Goodybe ;)")

