import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import argparse

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

def count_proportion_retained(SpeciesList:str, Raw:str, Filtered:str)->pd.DataFrame:
    """
    Count the number of ORFs before and after filtering and return it in a df

    Parameters: 
    - Filepath (str): Path to the file with the sample ids
    - Raw (str): Path to the DESwoMAN folder (before filter)
    - Filtered (str): Path to the DESwoMAN folder (after filter)

    Returns:
    -pd.DataFrame: dataframe with the proportion per population
    """
    pops = get_populations(SpeciesList)
    OutDict = {}
    for pop in pops:
        RawNum = len(SeqIO.to_dict(SeqIO.parse(f"{Raw}/{pop}/denovo_protein.fa", "fasta")))
        FilteredNum = len(SeqIO.to_dict(SeqIO.parse(f"{Filtered}/{pop}/denovo_protein.fa", "fasta")))
        PercRetained = (FilteredNum/RawNum) * 100
        OutDict[pop] =float('%.1f' % round(PercRetained, 1)) # Gives you '5.6'
    PercRetained
    df = pd.DataFrame(OutDict.items(), columns=['Line', 'Percent_kept'])
    return df

def make_plot(df:pd.DataFrame):
    """
    Plot the figure

    Parameters: 
    -pd.DataFrame: dataframe with the proportion per population

    """
    plt.figure(dpi = 200)
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)
    ax = sns.barplot(data = df, x = "Line", y = "Percent_kept")
    #The palette I use for my data
    #palette = 5*[ "#440154"] + 5*["#414487"] + 5*["#2a788e"] + 5*["#22a884"]+ 7*["#7ad151"]+10* ["#fde725"] 
    #for container in ax.containers:
    #    ax.bar_label(container, rotation = 70)
    plt.xticks(rotation=90)
    plt.ylabel("% neORFs kept after filtering")
    plt.xlabel(" ")
    plt.tight_layout()
    plt.savefig("Filtering.jpg")

def main():
    """
    Main function
    """
    parser = argparse.ArgumentParser(description="Plot your DESwoMAN filtering",epilog="-------------------------")
    parser.add_argument("--speciesFile", help="Path to the file with all populations/samples (.txt)", type=str)
    parser.add_argument("--raw", help="Path to the raw (!) DESwoMAN output", type=str)
    parser.add_argument("--filtered", help="Path to the filtered(!) DESwoMAN output", type=str)

    print("-------------------------\nPlot filtering output neORFs\nV.1.0\nAuthor:Marie Lebherz\n-------------------------\n")
    args = parser.parse_args()
    Specs = args.speciesFile
    Raw = args.raw
    Filt = args.filtered
    df = count_proportion_retained(Specs,Raw,Filt)
    make_plot(df)
    print("Done")

if __name__ == "__main__":
    main()
