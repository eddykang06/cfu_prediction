# Script to extract Fcnts data from a folder, CFU data from a file, then
# bind everything into 1 dataframe that can be used for training regression models
# 
# Input: 
#      
# Note that the input folder should be the default

# Packages
import pandas as pd
import numpy as np
import os
import glob
import sys

# Store command-line input and output
fcnts_path = sys.argv[0] # folder
cfu_path   = sys.argv[1]
outfig     = sys.argv[2]

def read_fcnts(folder_path):
    """
    Extracts all Fcnt csvs from a given folder as a list of dataframes
    Args:
        folder_path : file path containing Fcnts output from pipeline
    Output: 
        fcnt_df_list : list of Fcnts dataframes (converted from csvs)
    """

    # Extract all file names within folder
    files = os.listdir(folder_path)

    # Filter out the summary files and keep only NDC comparisons
    files = [csv for csv in files if ".summary" not in csv and "NDC0hr" in csv]

    # Attach path to each file
    filenames = ["".join([fcnts_path, "/" , csv]) for csv in files]

    # Load each file as a dataframe
    fcnt_df_list = [pd.read_table(csv, sep = "\t", header = 0, skiprows = 1) for csv in filenames]

    return fcnt_df_list

# Function to calculate TPMs from existing data
def tpm_convert(fcnt_df_list):
    """
    Converts a list of Fcnts dataframes to a single binded TPM dataframe
    Args:
        fcnt_df_list : 
    Output:
        tpm_df_list :

    """
    tpm_df_list = fcnt_df_list

    for df in tpm_df_list:

        # Move gene names to index
        df = df.set_index("Geneid")

        # Remove all columns except for Length and Fcnts
        df = df.loc[:, [col for col in list(df.columns) if col == "Length" or col.startswith("/")]]

        # Convert all entries to integers
        df = df.apply(lambda column: [int(entry) for entry in column])

        # Convert gene length from b -> kb
        df["Length"] = df["Length"].apply(lambda column: column / 1000)

        # Select Fcnts columns by excluding Length column
        fcnts_cols = [col fr col in list(df.columns) if col != "Length"]
        
        # Fcnts / gene length = (Counts per kb)
        df[fcnts_cols] = df[fcnts_cols].apply(lambda column: column / df["Length"])

        # (Counts per kb) * 10^6 / (Total counts/kb) = TPM
        df[fcnts_cols] = df[fcnts_cols].apply(lambda column: column * 10**6 / sum(column))

        # Remove length column 
        df.drop(columns = "Length", inplace = True)

        # Strip sample names to get easy to read samples
        def sample_name_strip(name):
            """
            Convert a sample file name into an easy to read sample name
            Input:
                name : (ex: "/ExpOut/260107_AV242502_RNASeq_miniHT_SpnT4WT_CEF_CIP/Out/Rep/Bams/T4-wt12CEF12CIP1hr-a.bam")
            Output:
                new_name : (ex: 12CEF12CIP1hr-a)
            """
            # Find index of the last / and remove entire prefix (OG file path)
            samplename_start_idx = name.strip().rfind("/") + 1
            new_name = name[samplename_start_idx:]

            # Find index of . (.bam is at end of sample name) and remove filetag
            filetag_start_idx = new_name.rfind(".")
            new_name = new_name[:filetag_start_idx]
    
            # Remove "T4-wt"
            new_name = new_name.replace("T4-wt", "")
    
            return new_name

        # Apply stripper to each column names
        df = df.rename(columns = lambda column: sample_name_strip(column))

    return tpm_df_list

# Function to bind all tpm df list into 1 df, then remove redundant NDC columns***
def bind_tpm_data(tpm_df_list):
    """
    Function to take a list of TPM dataframes, then bind all into 1 dataframe
    Args: 
        tpm_df_list : list of TPM dataframes
    Output:
        all_tpms : dataframe with all TPM values (samples on row, genes on column)
    """

    # Iterate a join by index (Reduce in R, outerjoin)


    all_tpms

    return all_tpms



# Function to calculate TPMs from


# Function to load the CFU data
def read_cfus(file_path):




# Function to bind TPMs
def bind_all_data():
    

# Call main
def main(fcnts_path, cfu_path, outfig):

    stored_fcnts = read_fcnts(fcnts_path)
    stored_tpms  = tpm_convert(stored_fcnts)
    tpm_df       = bind_tpm_data(stored_tpms)
    cfu_df       = read_cfus(cfu_path)
    all_data     = bind_all_data(tpm_df, cfu_df)

    return all_data

if __name__ == "__main__":
    main()