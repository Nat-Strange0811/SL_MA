import polars as pl
import sys
import os
import pandas as pd

'''
Main File -

This file creates temporary files for the meta-analysis run by 'METAL' software. It ensures all files
are in the same format and performs the filtering steps to save time.

It also includes a function to print the unique values contained within some of the columns, this is
to check that the original files are in the expected format and to check for any anomalies. This
function is not used in the main code but can be called if needed.
'''

def genFilePaths(seq_id, cohort_id):
    
    '''
    Function - genFilePaths
    
        Inputs:
            seq_id: An identifier for the protein that is being analysed, used to generate the file path
            cohort_id: An identifier for the cohort that is being analysed, used to generate the file path
            
        Output:
            The file path for the specified cohort and protein, this is used to read in the data
    '''
    
    #Base directory for all the files, this is the same for all cohorts and proteins
    base_dir = '/data/PHURI-Langenberg/people/Mine/SL_MA/'
    #Seq_id needs to be split into two parts to account for the different naming conventions used in the files
    seq1 = seq_id.split('_')[0]
    seq2 = seq_id.split('_')[1]
    
    #Dictionary used to store the file paths for each cohort
    file_paths = {
        "BWHHS_019" : os.path.join(base_dir, f"BWHHS_019/bwhhs019_case_{seq_id}_formatted.txt.gz"),
        "BWHHS_027" : os.path.join(base_dir, f"BWHHS_027/bwhhs027_{seq_id}_formatted.txt.gz"),
        "BWHHS_controls" : os.path.join(base_dir, f"BWHHS_controls/bwhhs019_control_{seq_id}_formatted.txt.gz"),
        "CHRIS" : os.path.join(base_dir, f"CHRIS/{seq1}-{seq2}_3_CHRIS_03102025_MF.tsv.gz"),
        "COPDGene_cases" : os.path.join(base_dir, f"COPDGene/output/cases/COPDcases_{seq_id}_formatted.txt.gz"),
        "COPDGene_controls" : os.path.join(base_dir, f"COPDGene/output/controls/COPDcontrols_{seq_id}_formatted.txt.gz"),
        "Decode" : os.path.join(base_dir, f"Decode/DECODE_{seq_id}_b37.txt.gz"),
        "EPIC_B1_B2_Other" : os.path.join(base_dir, f"EPIC_B1_B2_Other/all_invn_X{seq_id}_MarkerName_fastGWA.gz"),
        "EPIC_T2D_cases" : os.path.join(base_dir, f"EPIC_T2D/cases/cases_invn_X{seq_id}_MarkerName_fastGWA.gz"),
        "EPIC_T2D_cohort" : os.path.join(base_dir, f"EPIC_T2D/cohort/cohort_invn_X{seq_id}_MarkerName_fastGWA.gz"),
        "Fenland_GWAS_final" : os.path.join(base_dir, f"Fenland/Fenland_GWAS_final/output/Fenland_GWAS_res_invn_X{seq_id}_forMA.txt.gz"),
        "Fenland_OMICS_final" : os.path.join(base_dir, f"Fenland/Fenland_OMICS_final/output/Fenland_OMICS_res_invn_X{seq_id}_forMA.txt.gz"),
        "Fenland_CoreExome_final" : os.path.join(base_dir, f"Fenland/Fenland_CoreExome_final/output/Fenland_CoreExome_res_invn_X{seq_id}_forMA.txt.gz"),
        "Generation_Scotland" : os.path.join(base_dir, f"Generation_Scotland/genscot_{seq_id}_formatted.txt.gz"),
        "HUNT_controls" : os.path.join(base_dir, f"HUNT/controls/hunt_controls_{seq_id}_formatted.txt.gz"),
        "HUNT_incident" : os.path.join(base_dir, f"HUNT/incident/hunt_incident_{seq_id}_formatted.txt.gz"),
        "HUNT_prev" : os.path.join(base_dir, f"HUNT/prev/hunt_prev_{seq_id}_formatted.txt.gz"),
        "INTERVAL_SL_sumstats" : os.path.join(base_dir, f"INTERVAL_SL_sumstats/seq.{seq1}.{seq2}/{seq1}-{seq2}_3_INTERVAL_20250107_SCGP.tsv.gz"),
        "WHII" : os.path.join(base_dir, f"WHII/CLEANED.seq.{seq1}.{seq2}.fastGWA.gz")
    }
    
    return file_paths[cohort_id]

def printUniqueValues(data, ID):
    '''
    Function - printUniqueValues
    
        Inputs:
            data: A Polars DataFrame containing the data for a specific cohort
            ID: An identifier for the cohort, used to label the output
            
        Output:
            Prints the unique values for the chromosome, reference allele and alternate allele 
            columns, as well as a check for the validity of the SNPID format.
    '''
    
    #Initial print calls
    print()
    print(f"Unique values for {ID}:")
    
    #Regex pattern to check the SNPID column format
    pattern = r"^chr(?:\d+|X|Y|MT):\d+_[A-Z]+_[A-Z]+$"
    
    #Add a column for each unique component of the SNPID
    data = data.with_columns(
        pl.col("SNPID")
        .str.split(":", inclusive=False)
        .list.get(0)
        .str.replace("chr", "")
        .alias("chr"),

        pl.col("SNPID")
        .str.split(":", inclusive=False)
        .list.get(1)
        .str.split("_")
        .list.get(0)
        .alias("pos"),
        
        pl.col("SNPID")
        .str.split(":", inclusive=False)
        .list.get(1)
        .str.split("_")
        .list.get(1)
        .alias("ref"),
        
        pl.col("SNPID")
        .str.split(":", inclusive=False)
        .list.get(1)
        .str.split("_")
        .list.get(2)
        .alias("alt")
    )
    
    #Print all unique values in each of the columns above as well as checking the SNPID format
    print("Chromosomes:", data.select(pl.col("chr").unique().sort()).to_series().to_list())
    print("Ref Alleles:", data.select(pl.col("ref").unique().sort()).to_series().to_list())
    print("Alt Alleles:", data.select(pl.col("alt").unique().sort()).to_series().to_list())
    print("Valid SNPID format check:", data.select(pl.col("valid_snp").all()).to_series().item())
    
def main():
    
    '''
    Main function - main
    
        Inputs:
            None, but utilises command line arguments for the cohort ID and the seq_id
            
        Output:
            Creates a temporary file for the specified cohort and seq_id
    '''
    
    #Dictionary to store the column names for each cohort, this is used to read in the data and select the relevant columns
    ValuesDict = {
            "BWHHS_027"                 :   ("PVAL", "BETA", "MarkerName", "SE", "EAF_QTL", "INFO", "EFFECT_ALLELE", "OTHER_ALLELE"),
            "BWHHS_019"                 :   ("PVAL", "BETA", "MarkerName", "SE", "EAF_QTL", "INFO", "EFFECT_ALLELE", "OTHER_ALLELE"),
            "BWHHS_controls"            :   ("PVAL", "BETA", "MarkerName", "SE", "EAF_QTL", "INFO", "EFFECT_ALLELE", "OTHER_ALLELE"),
            "CHRIS"                     :   ("LOG10P", "BETA", "Harmonized_SNPID", "SE", "EAF", "INFO", "EFFECT_ALLELE", "NON_EFFECT_ALLELE"),
            "COPDGene_cases"            :   ("PVAL", "BETA", "MarkerName", "SE", "EAF", "INFO", "EFFECT_ALLELE", "NON_EFFECT_ALLELE"),
            "COPDGene_controls"         :   ("PVAL", "BETA", "MarkerName", "SE", "EAF", "INFO", "EFFECT_ALLELE", "NON_EFFECT_ALLELE"),
            "EPIC_B1_B2_Other"          :   ("P", "BETA", "MarkerName", "SE", "AF1", "INFO", "A1", "A2"),
            "EPIC_T2D_cohort"           :   ("P", "BETA", "MarkerName", "SE", "AF1", "INFO", "A1", "A2"),
            "EPIC_T2D_cases"            :   ("P", "BETA", "MarkerName", "SE", "AF1", "INFO", "A1", "A2"),
            "Fenland_OMICS_final"       :   ("res_invn_X-log10p", "res_invn_X_beta", "MarkerName", "res_invn_X_se", "af", "info", "a_0", "a_1"),
            "Fenland_GWAS_final"        :   ("res_invn_X-log10p", "res_invn_X_beta", "MarkerName", "res_invn_X_se", "af", "info", "a_0", "a_1"),
            "Fenland_CoreExome_final"   :   ("res_invn_X-log10p", "res_invn_X_beta", "MarkerName", "res_invn_X_se", "af", "info", "a_0", "a_1"),
            "Generation_Scotland"       :   ("PVAL", "BETA", "MarkerName", "SE", "EAF", "INFO", "EFFECT_ALLELE", "NON_EFFECT_ALLELE"),
            "INTERVAL_SL_sumstats"      :   ("MLOG10P", "BETA", "SNPID", "SE", "EAF", None, "EA", "NEA"),
            "WHII"                      :   ("PVAL", "BETA", "cpaid", "SE", "EAF", "INFO", "EFFECT_ALLELE", "OTHER_ALLELE"),
            "Decode"                    :   ("Pval", "Beta", "MarkerName", "SE", "ImpMAF", None, "effectAllele", "otherAllele"),
            "HUNT_controls"             :   ("P", "BETA", "MarkerName", "SE", "AF1", None, "A1", "A2"),
            "HUNT_prev"                 :   ("P", "BETA", "MarkerName", "SE", "AF1", None, "A1", "A2"),
            "HUNT_incident"             :   ("P", "BETA", "MarkerName", "SE", "AF1", None, "A1", "A2")
        }
    
    #Read system variables
    label = sys.argv[1]
    seq_id = sys.argv[2]
    
    #Generate the file path for the specified cohort and seq_id
    file = genFilePaths(seq_id, label)
    #Specify where to save the temporary file
    save_dir = f'/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Code/TemporaryFiles/{seq_id}/'
    
    #Make the save directory incase it doesn't exist
    os.makedirs(save_dir, exist_ok=True)
    
    #If we can't find the cohort file, skip the cohort, some cohorts don't have data for all proteins
    if not os.path.exists(file):
        print(f"File not found: {file}. Skipping {label}.")
        return
    
    #The EPIC B1_B2 files are space separated and therefore need to be treated uniquely
    if "B1_B2" in label:
        data = pl.scan_csv(file, separator=" ", null_values=["NA"])
    #Due to errors in the BWHHS files we need to manually define the column names and skip the first row
    elif "BWHHS" in label: 
        data = pl.scan_csv(file,
            separator="\t", 
            null_values=["NA", ""],
            has_header=False,
            new_columns=[
                "EFFECT_ALLELE",
                "OTHER_ALLELE",
                "EAF_QTL",
                "BETA",
                "SE",
                "PVAL",
                "INFO",
                "IMPUTED",
                "MarkerName",
                "_extra",
                "N"
            ]
            ).slice(1, None)
    #Otherwise just read normally
    else:
        data = pl.scan_csv(file, separator="\t", null_values=["NA"], truncate_ragged_lines=True)
    
    #Specify which column contains the Marker Name information, hence forth referred to as 'snp_col'
    snp_col = ValuesDict[label][2]
    
    #snp_col generally comes in two forms chr:pos_ref_alt or chr:pos:ref_alt
    
    #Ensure the snp_col is treated as a string
    snp = pl.col(snp_col).cast(pl.Utf8)
    #Split the snp_col by :, which seperates it into constituent parts. Chr and Pos_A1_A2
    parts = snp.str.split(":", inclusive=False)
    #We can ensure we get the part with the alleles by selecting the last element of the split.
    pos_alleles = parts.list.get(-1, null_on_oob = True).str.split("_")
    #The chromosome is the first part of the split, we also need to remove the 'chr' prefix and replace 'X' with '23' to ensure consistency across cohorts
    chrom = parts.list.get(0, null_on_oob = True).str.replace("chr", "").replace("X", "23")
    #The position is the first part of the second split
    pos = parts.list.get(1, null_on_oob = True).str.split("_").list.get(0)
    #we can then split the allele section to ensure we get the reference and alternate alleles
    ref_allele = pos_alleles.list.get(-2, null_on_oob = True)
    alt_allele = pos_alleles.list.get(-1, null_on_oob = True)
    
    #Define a pattern to ensure the constructed SNPID is correct
    pattern = r"^chr(?:\d+|X|Y|MT):\d+_[A-Z]+_[A-Z]+$"
    
    #Construct 'data' in the form we want, ensuring we only select relevant columns and construct the SNPID
    data = data.select(
        [
            (
                pl.lit("chr") + chrom + 
                pl.lit(":") + pos + 
                pl.lit("_") + ref_allele +
                pl.lit("_") + alt_allele
            ).alias("SNPID"),
            
            pl.col(ValuesDict[label][6]).alias("Effect_Allele"),
            pl.col(ValuesDict[label][7]).alias("Other_Allele"),
            pl.col(ValuesDict[label][1]).alias("Beta"),
            pl.col(ValuesDict[label][0]).alias("Pval"),
            pl.col(ValuesDict[label][3]).alias("SE"),
            pl.col(ValuesDict[label][4]).alias("EAF").cast(pl.Float64),
            pl.col(ValuesDict[label][5]).alias("INFO").cast(pl.Float64) if ValuesDict[label][5] is not None else pl.lit(1).alias("INFO"),
            pl.col("N").alias("N"),
        ]
    )
    
    #The WHII dataset contains some errors, where the rows are space seperated instead of tab seperated
    if label == "WHII":
        
        #Identify rows where the SNPID does not match the expected pattern
        data = data.with_columns(
            pl.col("SNPID").str.contains(pattern).alias("valid_snp")
        )
        
        #Separate the data into valid and invalid SNPIDs
        errors = data.filter(~pl.col("valid_snp"))
        data = data.filter(pl.col("valid_snp"))
        
        #For the rows with invalid SNPIDs, we need to split the SNPID column by space to extract the relevant information
        column = pl.col("SNPID")
        splitColumns = column.str.split(" ", inclusive=False)
        snpID = splitColumns.list.get(0, null_on_oob = True)
        eff = splitColumns.list.get(2, null_on_oob = True)
        other = splitColumns.list.get(3, null_on_oob = True)
        eaf  = splitColumns.list.get(4, null_on_oob = True)
        beta = splitColumns.list.get(5, null_on_oob = True)
        se = splitColumns.list.get(6, null_on_oob = True)
        pval = splitColumns.list.get(7, null_on_oob = True)
        n = splitColumns.list.get(8, null_on_oob = True)
        info = splitColumns.list.get(10, null_on_oob = True)
        
        fixed = errors.with_columns(
            snpID.alias("SNPID"),
            eff.alias("Effect_Allele"),
            other.alias("Other_Allele"),
            beta.cast(pl.Float64).alias("Beta"),
            pval.cast(pl.Float64).alias("Pval"),
            se.cast(pl.Float64).alias("SE"),
            eaf.cast(pl.Float64).alias("EAF"),
            info.cast(pl.Float64).alias("INFO"),
            n.cast(pl.Int64).alias("N")
        ).select([
            "SNPID",
            "Effect_Allele",
            "Other_Allele",
            "Beta",
            "Pval",
            "SE",
            "EAF",
            "INFO",
            "N"
        ])

        #In order to merge we need to have the same shapes
        fixed = fixed.with_columns(
            pl.col("SNPID").str.contains(pattern).alias("valid_snp")
        )
        
        #Sometimes the INFO column contains non-numeric values, which causes errors when we try to cast it to float
        #In this case we just discard the fixed data, pending updates.
        try:
            fixed.collect()
            data = pl.concat([data, fixed])
        except Exception as e:
            print(f"Error collecting fixed data: {e}")
            print("Skipping fixed data.")
    
    #Apply QC filters
    data = data.filter((pl.col("INFO") >= 0.4) & (pl.col("EAF") >= 0.001) & (pl.col("EAF") <= 0.999))
    
    #Select important columns for METAL
    data = data.select([
        "SNPID",
        "Effect_Allele",
        "Other_Allele",
        "Beta",
        "Pval",
        "SE",
        "N"
    ])
    
    #Collect the data into memory
    data = data.collect()
    
    #Save as a tab separated file with gzip compression
    data.write_csv(
        save_dir + f"{label}.txt.gz",
        separator="\t",
        include_header=True,
    )
        
if __name__ == "__main__":
    main()
        