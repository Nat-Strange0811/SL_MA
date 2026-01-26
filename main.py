import sys
import os
import random
import duckdb as dd
import pandas as pd
from pathlib import Path
import sys
import polars as pl
import re

'''
Script for initial QC analysis for SomaLogic Meta Analaysis project. Two main components:

    Class - Analysis:
        
        Class that is instantiated by the main method, holds all major details for each cohort and executes the random file selection and subsequent analysis. Purpose is to generate and write a csv file containing data for each analysed file within the total cohort.
        
    Method - Main:
    
        Checks that usage is appropriate and initiated the analysis class, serves as a launching point for the entire analysis.
'''

class Analysis():
    '''
    Analysis Class:
    
        Methods:
            init - Defines sturcture and launches/saves analysis
            calculateQCMetrics - Runs analysis protocol, generates a pandas dataframe containing all results
    '''
    
    
    def __init__(self, file_path, ID):
        '''
        init Method:
        
            Inputs:
                Self - The object
                file_path - File path that holds all the .gz files on the HPC
                ID - Cohort ID, defines particular behaviour and quirks
                
            Output:
                Saved .csv file for the cohort that is being run.
        '''
        
        #Tracking
        print("Analysis object initalised\n\n\n")
        
        #Define initial variables
        self.file_path = file_path
        self.ID = ID
        self.outFile = "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results"
        
        '''
        self.columndDict = {
            "BWHHS_027"         : ("EFFECT_ALLELE", "OTHER_ALLELE", "EAF_QTL", "BETA", "SE", "PVAL", "INFO", "IMPUTED", "MarkerName", "N"),
            "BWHHS_019"         : ("EFFECT_ALLELE", "OTHER_ALLELE", "EAF_QTL", "BETA", "SE", "PVAL", "INFO", "IMPUTED", "MarkerName", "N"),
            "CHRIS"             : ("CHR", "POS", "SNPID", "NON_EFFECT_ALLELE", "EFFECT_ALLELE", "EAF", "BETA", "SE", "LOG10P", "INFO", "N", "Harmonized_SNPID"),
            "COPD_cases"        : ("EFFECT_ALLELE", "NON_EFFECT_ALLELE", "EAF", "BETA", "SE", "PVAL", "INFO", "N", "MACH_R2", "MarkerName"),
            "COPD_controls"     : ("EFFECT_ALLELE", "NON_EFFECT_ALLELE", "EAF", "BETA", "SE", "PVAL", "INFO", "N", "MACH_R2", "MarkerName"),
            "EPIC"              : ("CHR", "SNP", "POS", "A1", "A2", "N", "AF1", "BETA", "SE", "P", "INFO", "MarkerName"),
            "EPIC_T2D"          : ("CHR", "SNP", "POS", "A1", "A2", "N", "AF1", "BETA", "SE", "P", "INFO", "MarkerName"),
            "EPIC_Cases"        : ("CHR", "SNP", "POS", "A1", "A2", "N", "AF1", "BETA", "SE", "P", "INFO", "MarkerName"),
            "Fenland_Omics"     : ("chr", "rsid", "pos", "a_0", "a_1", "af", "info", "res_invn_X_beta", "res_invn_X_se", "res_invn_X_t", "res_invn_X-log10p", "N", "MarkerName"),
            "Fenland_GWAS"      : ("chr", "rsid", "pos", "a_0", "a_1", "af", "info", "res_invn_X_beta", "res_invn_X_se", "res_invn_X_t", "res_invn_X-log10p", "N", "MarkerName"),
            "Fenland_CoreExome" : ("chr", "rsid", "pos", "a_0", "a_1", "af", "info", "res_invn_X_beta", "res_invn_X_se", "res_invn_X_t", "res_invn_X-log10p", "N", "MarkerName"),
            "GenScot"           : ("EFFECT_ALLELE", "NON_EFFECT_ALLELE", "EAF", "BETA", "SE", "PVAL", "INFO", "N", "MarkerName"),
            "INTERVAL"          : ("CHR", "POS", "SNPID", "EA", "NEA", "EAF", "N", "BETA", "SE", "MLOG10P", "CHISQ"),
            "WHII"              : ("cpaid", "SNP", "STRAND", "EFFECT_ALLELE", "OTHER_ALLELE", "EAF", "BETA", "SE", "PVAL", "N", "MAC", "INFO"),
            "VIKING"            : ("BETA", "SE", "PVAL", "SNPID", "Harmonised_SNPID", "chr", "POS", "EFFECT_ALLELE", "NON_EFFECT_ALLELE", "EAF", "INFO", "N", "IMPUTED", "HWE_P", "CALL_RATE"),
            "Decode"            : ("Chrom_b37", "Pos_b37", "Chrom_b38", "Pos_b38", "Name", "rsids", "effectAllele", "otherAllele", "Beta", "Pval", "minus_log10_pval", "SE", "N", "ImpMAF", "MarkerName"),
            "HUNT_controls"     : ("A1", "A2", "N", "AF1", "BETA", "SE", "P", "SNP", "Freq_Prev", "N_CHR_Prev", "Freq_Control", "N_CHR_Control", "Freq_Incident", "N_CHR_Incident", "r2", "MarkerName"),
            "HUNT_prev"         : ("A1", "A2", "N", "AF1", "BETA", "SE", "P", "SNP", "Freq_Prev", "N_CHR_Prev", "Freq_Control", "N_CHR_Control", "Freq_Incident", "N_CHR_Incident", "r2", "MarkerName"),
            "HUNT_incident"     : ("A1", "A2", "N", "AF1", "BETA", "SE", "P", "SNP", "Freq_Prev", "N_CHR_Prev", "Freq_Control", "N_CHR_Control", "Freq_Incident", "N_CHR_Incident", "r2", "MarkerName")
        }
        '''
        
        self.ValuesDict = {
            "BWHHS_027"         :   ("SE", "EAF_QTL", "Extracted", "EFFECT_ALLELE", "OTHER_ALLELE", "MarkerName", "INFO"),
            "BWHHS_019"         :   ("SE", "EAF_QTL", "Extracted", "EFFECT_ALLELE", "OTHER_ALLELE", "MarkerName", "INFO"),
            "BWHHS_Controls"    :   ("SE", "EAF_QTL", "Extracted", "EFFECT_ALLELE", "OTHER_ALLELE", "MarkerName", "INFO"),
            "CHRIS"             :   ("SE", "EAF", "Constructed", "CHR", "POS", "EFFECT_ALLELE", "NON_EFFECT_ALLELE", "INFO"),
            "COPD_cases"        :   ("SE", "EAF", "Extracted", "EFFECT_ALLELE", "NON_EFFECT_ALLELE", "MarkerName", "INFO"),
            "COPD_controls"     :   ("SE", "EAF", "Extracted", "EFFECT_ALLELE", "NON_EFFECT_ALLELE", "MarkerName", "INFO"),
            "EPIC"              :   ("SE", "AF1", "Constructed", "CHR", "POS", "A1", "A2", "INFO"),
            "EPIC_T2D"          :   ("SE", "AF1", "Constructed", "CHR", "POS", "A1", "A2", "INFO"),
            "EPIC_Cases"        :   ("SE", "AF1", "Constructed", "CHR", "POS", "A1", "A2", "INFO"),
            "Fenland_Omics"     :   ("res_invn_X_se", "af", "Constructed", "chr", "pos", "a_0", "a_1", "info"),
            "Fenland_GWAS"      :   ("res_invn_X_se", "af", "Constructed", "chr", "pos", "a_0", "a_1", "info"),
            "Fenland_CoreExome" :   ("res_invn_X_se", "af", "Constructed", "chr", "pos", "a_0", "a_1", "info"),
            "GenScot"           :   ("SE", "EAF", "Extracted", "EFFECT_ALLELE", "NON_EFFECT_ALLELE", "MarkerName", "INFO"),
            "INTERVAL"          :   ("SE", "EAF", "Constructed", "CHR", "POS", "EA", "NEA", None),
            "WHII"              :   ("SE", "EAF", "Extracted", "EFFECT_ALLELE", "OTHER_ALLELE", "cpaid", "INFO"),
            "VIKING"            :   ("SE", "EAF", "Extracted", "EFFECT_ALLELE", "NON_EFFECT_ALLELE", "Harmonised_SNPID", "INFO"),
            "Decode"            :   ("SE", "ImpMAF", "Constructed", "Chrom_b37", "Pos_b37", "effectAllele", "otherAllele", None),
            "HUNT_controls"     :   ("SE", "AF1", "Extracted", "A1", "A2", "MarkerName", None),
            "HUNT_prev"         :   ("SE", "AF1", "Extracted", "A1", "A2", "MarkerName", None),
            "HUNT_incident"     :   ("SE", "AF1", "Extracted", "A1", "A2", "MarkerName", None)
        }
         
        #Tracking
        print("Loading dbSNP\n\n\n")
        
        #Set up duckdb connection
        self.conn = dd.connect(database=':memory:')
        
        #Read in dbSNP
        snp = self.conn.read_csv(
            '/data/PHURI-Langenberg/people/SL_MA/02_dbSNP/00-All_filtered.tsv.gz',
            delimiter='\t',
            header=True,
            comment='#',
            skiprows=56,
            columns={
                'chr': 'VARCHAR',
                'pos': 'INTEGER',
                'rsID': 'VARCHAR',
                'ref': 'VARCHAR',
                'alt': 'VARCHAR',

            }
        )
        
        #Register snp table
        self.conn.register("snp", snp)
        
        #Tracking
        print("Loading files\n\n\n")
        
        #Interval folder structure is slightly unique, so to load all appropriate files we need a nuanced method. Utilises os.walk to check all sub directories in the original folder.
        if ID == "INTERVAL":
            self.files = []
            for root, _, files in os.walk(self.file_path):
                for f in files:
                    if f.endswith(".gz"):
                        self.files.append(os.path.join(root, f))   
        #All other IDs are relatively simple to handle     
        else:
            self.files = [os.path.join(self.file_path, f) for f in os.listdir(self.file_path) if f.endswith(".gz")]
        
        #Define n to avoid errors incase of using dummy data, generate a random subset of files
        n = min(100, len(self.files))
        self.random_files = random.sample(self.files, n)
        
        #Instantiate results
        self.results = pd.DataFrame(columns= ["FileID", "SE", "MAC", "MAF", "Rare Variant", "INFO", "rsID"])
        self.length = 0
        
        #Tracking
        print("Running results\n\n\n")
        
        #Run analysis
        self.calculateQCMetrics()
        
        #Save results
        self.results.to_csv(os.path.join(self.outFile, f"{self.ID}.csv"), index=False)    
    
    def ExtractID(self, file):
        '''
        ExtractID Method:
        
            Inputs:
                self - Allows access to object level variables
                file - File path for the specific file being analysed
                
            Outputs:
                Extracted ID from the file name
        '''
        
        regex_dict = {
            "Decode"            : r'^DECODE_(?P<seq_id>[^_]+(?:_[^_]+)*)_b37*$',
            "CHRIS"             : r'^(?P<seq_id>[^_]+)_.*$',
            "INTERVAL"          : r'^(?P<seq_id>[^_]+)_.*$',
            "BWHHS_027"         : r'bwhhs027_(?P<seq_id>[^_]+_[^_]+)_formatted\.txt\.gz$',
            "BWHHS_019"         : r'^bwhhs019_case_(?P<seq_id>[^_]+_[^_]+)_formatted\.txt\.gz$',
            "BWHHS_Controls"    : r'^bwhhs019_control_(?P<seq_id>[^_]+_[^_]+)_formatted\.txt\.gz$',
            "COPD_cases"        : r'^COPDcases_(?P<seq_id>[^_]+_[^_]+)_formatted\.txt\.gz$',
            "COPD_controls"     : r'^COPDcontrols_(?P<seq_id>[^_]+_[^_]+)_formatted\.txt\.gz$',
            "EPIC"              : r'^all_invn_X(?P<seq_id>[^_]+_[^_]+)_MarkerName_fastGWA\.gz$',
            "GenScot"           : r'^genscot_(?P<seq_id>[^_]+_[^_]+)_formatted\.txt\.gz$',
            "HUNT_controls"     : r'^hunt_controls_(?P<seq_id>[^_]+_[^_]+)_formatted\.txt\.gz$',
            "HUNT_prev"         : r'^hunt_prev_(?P<seq_id>[^_]+_[^_]+)_formatted\.txt\.gz$',
            "HUNT_incident"     : r'^hunt_incident_(?P<seq_id>[^_]+_[^_]+)_formatted\.txt\.gz$',
            "EPIC_T2D"          : r'^cohort_invn_X(?P<seq_id>[^_]+_[^_]+)_MarkerName_fastGWA\.gz$',
            "EPIC_Cases"        : r'^cases_invn_X(?P<seq_id>[^_]+_[^_]+)_MarkerName_fastGWA\.gz$',
            "WHII"              : r'^CLEANED\.seq\.(?P<seq_id>[^.]+\.[^.]+)\.fastGWA\.gz$',
            "VIKING"            : r'^(?P<seq_id>[^_]+)_.*$',
            "Fenland_Omics"     : r'^Fenland_OMICS_res_invn_X(?P<seq_id>[^_]+_[^_]+)_.*$',
            "Fenland_GWAS"      : r'^Fenland_GWAS_res_invn_X(?P<seq_id>[^_]+_[^_]+)_.*$',
            "Fenland_CoreExome" : r'^Fenland_CoreExome_res_invn_X(?P<seq_id>[^_]+_[^_]+)_.*$'
        }    
        
        if "test" in file:
            return "test_file"
        
        ID = re.search(regex_dict[self.ID], file).group("seq_id")
        
        return ID
    
    def calculateQCMetrics(self):
        '''
        calculateQCMetric Method:
        
            Inputs:
                self - Allows access to object level variables
                
            Outputs:
                Constructs the self.results dataframe

        '''
        
        #Establish initial mean values for the cohort
        mean_se_pct = 0
        mean_maf_pct = 0
        mean_mac_pct= 0
        
        if self.ValuesDict[self.ID][-1]:
            mean_info_pct = 0
        else:
            mean_info_pct = None
        
        mean_rarevar_pct = 0    
        mean_rsID_pct = 0
        
        #Loop over all files selected
        for file in self.random_files:
            filename = Path(file).name
            
            #Tracking
            print(f"Running for file {self.length + 1}")
            
            #Establish initial counts for number of lines that match given condition
            se_count = 0
            maf_count = 0
            mac_count = 0
            rarevar_count = 0
            n_rows = 0
            
            if self.ValuesDict[self.ID][-1]:
                info_count = 0
            else:
                info_count = None
                
            rsID_count = 0
            
            
            #Load the data file with polars
            if self.ID == "EPIC":
                data = pl.scan_csv(file, separator=" ", null_values=["NA"])
            elif "BWHHS" in self.ID:
                #Due to errors in the BWHHS files we need to manually define the column names and skip the first row
                data = pl.scan_csv("/data/PHURI-Langenberg/people/Mine/SL_MA/BWHHS_019/bwhhs019_case_2190_55_formatted.txt.gz", 
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
            else:
                data = pl.scan_csv(file, separator="\t", null_values=["NA"])
            
            #Count total number of rows    
            n_rows = data.select(pl.len()).collect().item()
                
            #Viking cohort needs to be merged
            if self.ID == "VIKING":
                snps = pl.scan_csv(os.path.join(self.file_path, "snpID.txt.tmp"), separator="\t", null_values=["NA"])
                temp = data
                data = pl.concat([snps, temp], how="horizontal")
                
            
            #Calculate QC metrics
            #Check if info is available
            if self.ValuesDict[self.ID][-1]:
                #Instantiate a new df with the required columns, with True/False values for each condition
                res = data.select([
                    (
                        pl.col(self.ValuesDict[self.ID][0]).cast(pl.Float64) > 10
                    ).alias("se"),
                    (
                        pl.min_horizontal(
                            pl.col(self.ValuesDict[self.ID][1]).cast(pl.Float64),
                            1 - pl.col(self.ValuesDict[self.ID][1]).cast(pl.Float64),
                        ) > 0.001
                    ).alias("maf"),
                    (
                        2
                        * pl.col("N").cast(pl.Float64)
                        * pl.min_horizontal(
                            pl.col(self.ValuesDict[self.ID][1]).cast(pl.Float64),
                            1 - pl.col(self.ValuesDict[self.ID][1]).cast(pl.Float64),
                        )
                        > 3
                    ).alias("mac"),
                    (
                        pl.col(self.ValuesDict[self.ID][-1]).cast(pl.Float64) < 0.8
                    ).sum().alias("info")
                ]).collect()
                
                #Update count for how many rows match the Info condition
                info_count += res["info"][0]
                
            else:
                #Same as above but without Info
                res = data.select([
                    (
                        pl.col(self.ValuesDict[self.ID][0]).cast(pl.Float64) > 10
                    ).alias("se"),
                    (
                        pl.min_horizontal(
                            pl.col(self.ValuesDict[self.ID][1]).cast(pl.Float64),
                            1 - pl.col(self.ValuesDict[self.ID][1]).cast(pl.Float64),
                        ) > 0.001
                    ).alias("maf"),
                    (
                        2
                        * pl.col("N").cast(pl.Float64)
                        * pl.min_horizontal(
                            pl.col(self.ValuesDict[self.ID][1]).cast(pl.Float64),
                            1 - pl.col(self.ValuesDict[self.ID][1]).cast(pl.Float64),
                        )
                        > 3
                    ).alias("mac"),
                ]).collect()
            
            #Select and sum the True/False values to get counts, adds a new column "rarevar_count" for the number of 'false' in se caused by rare variants
            res = res.select([
                pl.sum("se").alias("se_count"),
                pl.sum("maf").alias("maf_count"),
                pl.sum("mac").alias("mac_count"),
                ((pl.col("se") == True) & ((pl.col("maf") == False) | (pl.col("mac") == False))).sum().alias("rarevar_count")
            ])
            
            #Update counts
            se_count += res["se_count"][0]
            maf_count +=  res["maf_count"][0]
            mac_count +=  res["mac_count"][0]
            rarevar_count += res["rarevar_count"][0]
            
            
            #Check which mode for SNP details extraction is required
            if self.ValuesDict[self.ID][2] == "Constructed":
                #Extract chr, pos, ref, alt from the relevant columns
                queries = data.select([
                    pl.col(self.ValuesDict[self.ID][3]).cast(pl.Utf8).str.replace("chr", "").alias("chr"),
                    pl.col(self.ValuesDict[self.ID][4]).cast(pl.Int32).alias("pos"),
                    pl.col(self.ValuesDict[self.ID][5]).alias("ref"),
                    pl.col(self.ValuesDict[self.ID][6]).alias("alt"),
                ]).collect()
            else:
                snp_col = self.ValuesDict[self.ID][5]

                queries = data.select([
                    pl.col(snp_col)
                    .str.split(":", inclusive=False)
                    .list.get(0)
                    .str.replace("chr", "")
                    .alias("chr"),

                    pl.col(snp_col)
                    .str.split(":", inclusive=False)
                    .list.get(1)
                    .str.split("_")
                    .list.get(0)
                    .cast(pl.Int32)
                    .alias("pos"),

                    pl.col(self.ValuesDict[self.ID][3]).alias("ref"),
                    pl.col(self.ValuesDict[self.ID][4]).alias("alt"),
                ]).collect()
        
            self.conn.register('data', queries)
            
            rsID_count = self.conn.execute("""
                SELECT COUNT(*) AS match_count
                FROM data
                JOIN snp
                ON data.chr = TRIM(snp.chr)
                AND data.pos = snp.pos
                AND (
                    (data.ref = ANY(string_split(snp.ref, ','))
                    AND data.alt = ANY(string_split(snp.alt, ',')))
                OR (data.ref = ANY(string_split(snp.alt, ','))
                    AND data.alt = ANY(string_split(snp.ref, ',')))
                )
            """).fetchone()[0]

            #Calculate percentages
            se_pct = (n_rows - se_count) / n_rows * 100
            maf_pct = maf_count / n_rows * 100
            mac_pct = mac_count / n_rows * 100
            rarevar_pct = rarevar_count / se_count * 100 if se_count != 0 else 100
            info_pct = (n_rows - info_count) / n_rows * 100 if info_count is not None else None
            rsID_pct = rsID_count / n_rows * 100
            
            #Append to results
            self.results.loc[self.length] = [self.ExtractID(filename), 
                                            se_pct, 
                                            mac_pct, 
                                            maf_pct, 
                                            rarevar_pct,
                                            info_pct, 
                                            rsID_pct
                                            ]
            self.length += 1
            
            #Update cohort averages
            mean_se_pct += se_pct
            mean_maf_pct += maf_pct
            mean_mac_pct += mac_pct
            mean_rarevar_pct += rarevar_pct
            if self.ValuesDict[self.ID][-1]:
                mean_info_pct += info_pct
            mean_rsID_pct += rsID_pct
        
        #Calculate cohort averages    
        mean_se_pct /= self.length
        mean_maf_pct /= self.length
        mean_mac_pct /= self.length
        mean_rarevar_pct /= self.length
        if self.ValuesDict[self.ID][-1]:
            mean_info_pct /= self.length
        mean_rsID_pct /= self.length
            
        self.results.loc[self.length] = ["Average", 
                                         mean_se_pct, 
                                         mean_mac_pct, 
                                         mean_maf_pct, 
                                         mean_rarevar_pct,
                                         mean_info_pct, 
                                         mean_rsID_pct
                                         ]    
            
        
    
def main():
    
    if len(sys.argv) < 3:
        print("Usage: python3 script.py <file_path> <ID>")
        sys.exit(1)
        
    analysis = Analysis(sys.argv[1], sys.argv[2])
    


if __name__ == "__main__":
    main()