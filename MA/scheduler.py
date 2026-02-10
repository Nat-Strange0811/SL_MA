import subprocess
import sys
import os

'''

scheduler.py

This file acts as a python scheduler for the entire MA process, it uses multi-core processing to
schedule the execution of the main.py script for each cohort in parallel, it then uses built in
python function to write the command file for the metal executable and calls the metal executable
to the commandline, finally it calls the delete.py script for cleanup.

This has the benefit of using seperate .sh files as it means fewer 'jobs' need to be submitted to
the cluster avoiding some of the limits on the number of jobs that each user can have. However, 
due to the requirement of lots of memory for the main.py process, it causes long queue times. It
has been kept as an example of how the process can be scheduled in python. 

'''


def main():
    
    '''
    
    Main Function:
    
        Inputs:
            - seq_id: The ID of the protein for which the meta-analysis is to be performed, accessed through command line.
            
        Outputs:
            - Executes the main.py script for each cohort in parallel, then executes the metal executable, and finally calls the delete.py script for cleanup.
    
    '''
    
    # Access the protein ID from the command line arguments
    seq_id = sys.argv[1]

    # Define the list of cohorts to be processed in the meta-analysis
    cohorts = [
        "BWHHS_027",
        "BWHHS_019",
        "BWHHS_controls",
        "CHRIS",
        "COPDGene_cases",
        "COPDGene_controls",
        "EPIC_B1_B2_Other",
        "EPIC_T2D_cohort",
        "EPIC_T2D_cases",
        "Fenland_OMICS_final",
        "Fenland_GWAS_final",
        "Fenland_CoreExome_final",
        "Generation_Scotland",
        "INTERVAL_SL_sumstats",
        "WHII",
        "Decode",
        "HUNT_controls",
        "HUNT_prev",
        "HUNT_incident"
    ]
    
    # Initialize a list to keep track of the subprocesses for each cohort
    procs = []
    
    # Loop through each cohort
    for cohort in cohorts:
        #Initialise the log and error file paths for the current cohort, ensure they exist
        cohort_log = f"/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Code/LogsMA/{seq_id}/Tempfiles/{cohort}_log.txt"
        cohort_err = f"/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Code/LogsMA/{seq_id}/Tempfiles/{cohort}_err.txt"
        os.makedirs(os.path.dirname(cohort_log), exist_ok=True)
        os.makedirs(os.path.dirname(cohort_err), exist_ok=True)
        
        #Use 'subprocess' to call the main.py script for the current cohort.
        #Use open to define two files with which python can write the output and error streams to
        #subprocess.Popen takes a list of command line arguments, as well as optional arguments for output files
        with open(cohort_log, "w") as log_file, open(cohort_err, "w") as err_file:
            p = subprocess.Popen(['python', 'MA/main.py', cohort, seq_id], stdout=log_file, stderr=err_file)
            #Append the process to the list of processes so we can wait for them to finish later
            procs.append(p)
    
    #For each process (at the moment creating temporary files for each cohort) is waited on to finish before proceeding    
    for p in procs:
        p.wait()
    
    #Make a temporary file for the command file for the metal executable
    os.makedirs(os.path.dirname(f"/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Code/TemporaryFiles/{seq_id}/meta_analysis_commands.txt"), exist_ok=True)
    #Open the command file
    with open(f"/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Code/TemporaryFiles/{seq_id}/meta_analysis_commands.txt", "w") as f:
        #Write a list of commands, these are 'standard' and do not change
        f.write("SCHEME STDERR\n")
        f.write("MARKER SNPID\n")
        f.write("ALLELE Effect_Allele Other_Allele\n")
        f.write("EFFECT Beta\n")
        f.write("PVALUE Pval\n")
        f.write("WEIGHT N\n")
        f.write("STDERR SE\n")
        f.write("SEPARATOR TAB\n")
        
        #Loop through all temporary files created and append them to the command file
        for file in os.listdir(f"/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Code/TemporaryFiles/{seq_id}/"):
            if file.endswith(".txt.gz"):
                f.write(f"PROCESS /data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Code/TemporaryFiles/{seq_id}/{file}\n")

        #Write the final commads for the metal executable
        f.write(f"OUTFILE /data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Code/MA/Results/{seq_id}/Meta_Aanalysis_{seq_id}_results .tbl\n")
        f.write("ANALYZE HETEROGENEITY\n")
        f.write("QUIT\n")
    
    #Define the metal log and error file paths and ensure they exist
    metal_log = f"/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Code/LogsMA/{seq_id}/MetalOutput/metal_log.txt"   
    metal_err = f"/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Code/LogsMA/{seq_id}/MetalOutput/metal_err.txt" 
    os.makedirs(os.path.dirname(metal_log), exist_ok=True)
    os.makedirs(os.path.dirname(metal_err), exist_ok=True)
    #Call the metal executable using subprocess, passing the command file as an argument
    with open(metal_log, "w") as log_file, open(metal_err, "w") as err_file:
        p = subprocess.Popen(["/data/PHURI-Langenberg/programs/METAL/random-metal-0.1.0/executables/metal", f"/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Code/TemporaryFiles/{seq_id}/meta_analysis_commands.txt"], stdout=log_file, stderr=err_file)
        p.wait()
    
    #Define the cleanup log and error file paths and ensure they exist
    cleanup_log = f"/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Code/LogsMA/{seq_id}/Cleanup/cleanup_log.txt"
    cleanup_err = f"/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Code/LogsMA/{seq_id}/Cleanup/cleanup_err.txt"
    os.makedirs(os.path.dirname(cleanup_log), exist_ok=True)
    os.makedirs(os.path.dirname(cleanup_err), exist_ok=True)
    #Call the cleanup script using subprocess, passing the protein ID as an argument
    with open(cleanup_log, "w") as log_file, open(cleanup_err, "w") as err_file:
        q = subprocess.Popen(['python', 'MA/delete.py', seq_id], stdout=log_file, stderr=err_file)
        q.wait()
    
if __name__ == "__main__":
    main()