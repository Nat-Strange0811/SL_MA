import os
import sys
import shutil
import duckdb as dd
import pandas as pd

'''
delete.py

The purpose of this script is to delete temporary files created during the MA process for a specified
protein ID. The input to the script is the ID which is accessed through the command line.

'''

def main():
    '''
    Main Function:
    
        Inputs:
            - seq_id: The ID of the protein for which temporary files are to be deleted.
            
        Outputs:
            - Deletes the temporary files associated with the specified protein ID.
    '''
    
    print("Accessing system arguments...")
    # Access the protein ID from the command line arguments
    seq_id = sys.argv[1]
    mode = sys.argv[2]
    
    #Define the directory path where the temporary files are stored for the given protein ID
    save_dir = f'/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Code/TemporaryFiles/{seq_id}/'
    
    # Create an empty list to store the dataframes read from the temporary files
    dataframes = []
    
    print("Reading temporary files and creating lookup table...")
    # Read each temporary file in the specified directory and extract the relevant columns to create a list of dataframes
    for file in os.listdir(save_dir):
        dataframes.append(pd.read_csv(os.path.join(save_dir, file), sep='\t', usecols=['SNPID', 'dbSNP_ID', 'Chr', 'Pos'], low_memory=False))
    
    print("Concatenating dataframes")    
    # Concatenate all the dataframes in the list into a single dataframe, ignoring the index to create a continuous index
    temp_files = pd.concat(dataframes, ignore_index=True)
    
    print("Removing duplicate rows to create lookup table...")
    # Remove duplicate rows from the concatenated dataframe based on the specified columns to create a lookup table for dbSNP_ID, Chr, and Pos
    lookup = temp_files.drop_duplicates(subset=['SNPID'])
    
    # Establish a connection to an in-memory DuckDB database
    conn = dd.connect(database=':memory:')
    
    print("Registering lookup table...")
    # Register the lookup dataframe as a table in the DuckDB database to enable SQL queries on it
    conn.register('lookup', lookup)
    
    print("Registering MetaAnalysis results...")
    # Read the MetaAnalysis result file for the specified protein ID and mode, and register it as a table in the DuckDB database to enable SQL queries on it
    conn.register('ma', pd.read_csv(f'/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Code/MA/Results/{seq_id}/MetaAnalysis_{seq_id}_{mode}_1.tbl.gz', sep='\t'))
    
    print("Executing SQL query to merge MetaAnalysis results with lookup table...")
    result = conn.execute("""
        SELECT 
            ma.MarkerName,
            lookup.Chr,
            lookup.Pos,
            lookup.dbSNP_ID,
            ma.Allele1,
            ma.Allele2,
            ma.Freq1,
            ma.FreqSE,
            ma.MinFreq,
            ma.MaxFreq,
            ma.Effect,
            ma.StdErr,
            ma.Pvalue,
            ma.Direction,
            ma.HetISq,
            ma.HetChiSq,
            ma.HetDf,
            ma.HetPVal,
            ma.TotalSampleSize
        FROM ma
        INNER JOIN lookup ON ma.MarkerName = lookup.SNPID
    """).df()
    
    print("Sorting results...")
    result = result.sort_values(by=['Chr', 'Pos'])
    
    print("Removing old MetaAnalysis result file...")
    os.remove(f'/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Code/MA/Results/{seq_id}/MetaAnalysis_{seq_id}_{mode}_1.tbl.gz')
    result.to_csv(f'/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Code/MA/Results/{seq_id}/MetaAnalysis_{seq_id}_{mode}.tbl.gz', sep='\t', index=False, compression='gzip')
    
    #Remove the directory and all its contents
    shutil.rmtree(save_dir)         

if __name__ == "__main__":
    main()