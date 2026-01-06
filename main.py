import sys
import os
import random
import pandas as pd
from pathlib import Path
import sys
import polars as pl
import hashlib

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
            "Decode"            : ("Chrom", "Pos", "Name", "rsids", "effectAllele", "otherAllele", "Beta", "Pval", "minus_log10_pval", "SE", "N", "ImpMAF"),
            "HUNT_controls"     : ("A1", "A2", "N", "AF1", "BETA", "SE", "P", "SNP", "Freq_Prev", "N_CHR_Prev", "Freq_Control", "N_CHR_Control", "Freq_Incident", "N_CHR_Incident", "r2", "MarkerName"),
            "HUNT_prev"         : ("A1", "A2", "N", "AF1", "BETA", "SE", "P", "SNP", "Freq_Prev", "N_CHR_Prev", "Freq_Control", "N_CHR_Control", "Freq_Incident", "N_CHR_Incident", "r2", "MarkerName"),
            "HUNT_incident"     : ("A1", "A2", "N", "AF1", "BETA", "SE", "P", "SNP", "Freq_Prev", "N_CHR_Prev", "Freq_Control", "N_CHR_Control", "Freq_Incident", "N_CHR_Incident", "r2", "MarkerName")
        }
        '''
        
        self.ValuesDict = {
            "BWHHS_027"         :   ("SE", "EAF_QTL", "Extracted", "EFFECT_ALLELE", "OTHER_ALLELE", "MarkerName", "INFO"),
            "BWHHS_019"         :   ("SE", "EAF_QTL", "Extracted", "EFFECT_ALLELE", "OTHER_ALLELE", "MarkerName", "INFO"),
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
            "Decode"            :   ("SE", "ImpMAF", "Constructed", "Chrom", "Pos", "effectAllele", "otherAllele", None),
            "HUNT_controls"     :   ("SE", "AF1", "Extracted", "A1", "A2", "MarkerName", None),
            "HUNT_prev"         :   ("SE", "AF1", "Extracted", "A1", "A2", "MarkerName", None),
            "HUNT_incident"     :   ("SE", "AF1", "Extracted", "A1", "A2", "MarkerName", None)
        }
         
        #Tracking
        print("Loading dbSNP\n\n\n")
    
        self.snp = set()
        
        self.snp_df = pl.scan_csv(
                        "/data/PHURI-Langenberg/people/SL_MA/02_dbSNP/00-All_filtered.tsv.gz",
                        separator="\t",
                        skip_rows = 56,
                    ).select(["chr", "pos", "ref", "alt"]
                    )
                    
        offset = 0
        batch_size = 1_000_000
        '''
        while True:
            batch = self.snp_df.slice(offset, batch_size).collect()
            if batch.is_empty():
                break
        
            for chr_, pos, ref, alt in batch.iter_rows():
                self.snp.add(self.key(chr_, pos, ref, alt))
                self.snp.add(self.key(chr_, pos, alt, ref))
                
            offset += batch_size
        '''
        del self.snp_df
        
        
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
        self.results = pd.DataFrame(columns= ["FileID", "SE", "MAC", "MAF", "INFO", "rsID"])
        self.length = 0
        
        #Tracking
        print("Running results\n\n\n")
        
        #Run analysis
        self.calculateQCMetrics()
        
        #Save results
        self.results.to_csv(os.path.join(self.outFile, f"{self.ID}.csv"), index=False)
        
    
    def key(self, chro, pos, ref, alt):
        return
        h = hashlib.blake2b(digest_size=16)
        h.update(str(chro).encode())
        h.update(b'\0')
        h.update(str(pos).encode())
        h.update(b'\0')
        h.update(ref.encode())
        h.update(b'\0')
        h.update(alt.encode())
        return h.digest()
        
    
    def calculateQCMetrics(self):
        '''
        calculateQCMetric Method:
        
            Inputs:
                self - Allows access to object level variables
                
            Outputs:
                Constructs the self.results dataframe

        '''
        
        mean_se_pct = 0
        mean_maf_pct = 0
        mean_mac_pct= 0
        
        if self.ValuesDict[self.ID][-1]:
            mean_info_pct = 0
        else:
            mean_info_pct = None
            
        mean_rsID_pct = 0
        
        #Loop over all files selected
        for file in self.random_files:
            filename = Path(file).name
            
            #Tracking
            print(f"Running for file {self.length + 1}")
            
            se_count = 0
            maf_count = 0
            mac_count = 0
            n_rows = 0
            
            if self.ValuesDict[self.ID][-1]:
                info_count = 0
            else:
                info_count = None
                
            rsID_count = 0
            
            
            #Initialise the data structure and where the QC results will be held
            data = pl.scan_csv(file, separator="\t", null_values=["NA"])
                
            rows = data.select(pl.len()).collect().item()
            n_rows += rows
                
            #Viking cohort needs to be merged
            if self.ID == "VIKING":
                snps = pl.scan_csv(os.path.join(self.file_path, "snpID.txt.tmp"), separator="\t", null_values=["NA"])
                temp = data
                data = pl.concat([snps, temp], how="horizontal")
                
            
            #Calculate QC metrics
            #Check if info is available
            if self.ValuesDict[self.ID][-1]:
                res = data.select([
                    (
                        pl.col(self.ValuesDict[self.ID][0]).cast(pl.Float64) < 10
                    ).sum().alias("se"),
                    (
                        pl.min_horizontal(
                            pl.col(self.ValuesDict[self.ID][1]).cast(pl.Float64),
                            1 - pl.col(self.ValuesDict[self.ID][1]).cast(pl.Float64),
                        ) > 0.001
                    ).sum().alias("maf"),
                    (
                        2
                        * pl.col("N").cast(pl.Float64)
                        * pl.min_horizontal(
                            pl.col(self.ValuesDict[self.ID][1]).cast(pl.Float64),
                            1 - pl.col(self.ValuesDict[self.ID][1]).cast(pl.Float64),
                        )
                        > 3
                    ).sum().alias("mac"),
                    (
                        pl.col(self.ValuesDict[self.ID][-1]).cast(pl.Float64) > 0.8
                    ).sum().alias("info")
                ]).collect()
                
                se_count += res["se"][0]
                maf_count += res["maf"][0]
                mac_count += res["mac"][0]
                info_count += res["info"][0]
                
            else:
                res = data.select([
                    (
                        pl.col(self.ValuesDict[self.ID][0]).cast(pl.Float64) < 10
                    ).sum().alias("se"),
                    (
                        pl.min_horizontal(
                            pl.col(self.ValuesDict[self.ID][1]).cast(pl.Float64),
                            1 - pl.col(self.ValuesDict[self.ID][1]).cast(pl.Float64),
                        ) > 0.001
                    ).sum().alias("maf"),
                    (
                        2
                        * pl.col("N").cast(pl.Float64)
                        * pl.min_horizontal(
                            pl.col(self.ValuesDict[self.ID][1]).cast(pl.Float64),
                            1 - pl.col(self.ValuesDict[self.ID][1]).cast(pl.Float64),
                        )
                        > 3
                    ).sum().alias("mac"),
                ]).collect()
                
                se_count += res["se"][0]
                maf_count += res["maf"][0]
                mac_count += res["mac"][0]
            
            #Check which mode for SNP details extraction is required
            if self.ValuesDict[self.ID][2] == "Constructed":
                queries = data.select([
                    pl.col(self.ValuesDict[self.ID][3]).alias("chr"),
                    pl.col(self.ValuesDict[self.ID][4]).cast(pl.Int64).alias("pos"),
                    pl.col(self.ValuesDict[self.ID][5]).alias("ref"),
                    pl.col(self.ValuesDict[self.ID][6]).alias("alt"),
                ])
            else:
                snp_col = self.ValuesDict[self.ID][5]

                queries = data.select([
                    pl.col(snp_col)
                    .str.split(":", inclusive=False)
                    .list.get(0)
                    .str.replace("chr", "")
                    .cast(pl.Int64)
                    .alias("chr"),

                    pl.col(snp_col)
                    .str.split(":", inclusive=False)
                    .list.get(1)
                    .str.split("_")
                    .list.get(0)
                    .cast(pl.Int64)
                    .alias("pos"),

                    pl.col(self.ValuesDict[self.ID][3]).alias("ref"),
                    pl.col(self.ValuesDict[self.ID][4]).alias("alt"),
                ])
            
            for batch in queries.collect(engine="streaming").iter_slices(n_rows=1_000_000):
                for chr_, pos, ref, alt in batch.iter_rows():
                    if self.key(chr_, pos, ref, alt) in self.snp:
                        rsID_count += 1

            
            se_pct = se_count / n_rows * 100
            maf_pct = maf_count / n_rows * 100
            mac_pct = mac_count / n_rows * 100
            info_pct = info_count / n_rows * 100 if info_count is not None else None
            rsID_pct = rsID_count / n_rows * 100
            
            self.results.loc[self.length] = [filename, 
                                            se_pct, 
                                            mac_pct, 
                                            maf_pct, 
                                            info_pct, 
                                            rsID_pct
                                            ]
            self.length += 1
            
            mean_se_pct += se_pct
            mean_maf_pct += maf_pct
            mean_mac_pct += mac_pct
            if self.ValuesDict[self.ID][-1]:
                mean_info_pct += info_pct
            mean_rsID_pct += rsID_pct
            
        mean_se_pct /= self.length
        mean_maf_pct /= self.length
        mean_mac_pct /= self.length
        if self.ValuesDict[self.ID][-1]:
            mean_info_pct /= self.length
        mean_rsID_pct /= self.length
            
        self.results.loc[self.length] = ["Average", 
                                         mean_se_pct, 
                                         mean_mac_pct, 
                                         mean_maf_pct, 
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