import sys
import shutil

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
    
    # Access the protein ID from the command line arguments
    seq_id = sys.argv[1]
    
    #Define the directory path where the temporary files are stored for the given protein ID
    save_dir = f'/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Code/TemporaryFiles/{seq_id}/'
    #Remove the directory and all its contents
    shutil.rmtree(save_dir)
            
if __name__ == "__main__":
    main()