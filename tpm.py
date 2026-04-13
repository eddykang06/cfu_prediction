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
    '''
    Extracts all Fcnt csvs from a given folder as a list of dataframes
    Args:
        folder_path : file path containing Fcnts output from pipeline
    Output: 
        fcnt_df_list : list of Fcnts dataframes (converted from csvs)
    '''

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
    '''
    Converts a list of Fcnts dataframes to a single binded TPM dataframe
    Input:
        fcnt_df_list : 
    Output:
        all_tpm_df :

    '''
    for df in fcnt_df_list:

        # Move gene names to index
        df = df.set_index("Geneid")

        # Remove all columns except for Length and Fcnts
        df = df.loc[:, [col for col in list(df.columns) if col == "Length" or col.startswith("/")]]

        # Convert all entries to integers
        df = df.apply(lambda column: [int(entry) for entry in column])

        # Convert gene length from b -> kb
        df["Length"] = df["Length"].apply(lambda column: column / 1000)

        # Select Fcnts columns and operate





        # Remember to 


        # 
        all_tpm_df = ?
    
    return all_tpm_df

# can use df.rename(lambda x: ) to rename the samples?



# Function to calculate TPMs from


# Function to load the CFU data
def read_cfus(file_path):




# Function to bind TPMs
def bind_data():

# Call main
def main(fcnts_path, cfu_path, outfig):

    stored_fcnts = read_fcnts(fcnts_path)
    stored_tpms  = tpm_convert(stored_fcnts)
    stored_cfus  = read_cfus()
    all_data     = bind_data(stored_fcnts, stored_cfus)
    # 

    return 

if __name__ == "__main__":
    main()