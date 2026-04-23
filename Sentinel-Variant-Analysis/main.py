import sys
import os
import pandas as pd
import glob
import math
from compare_variants import compare_variants



def main():
    one_cohort_only = [10379_19, 10724_45, 10837_131, 10928_7, 11890_2, 20437_9, 22134_1, 22532_59, 23372_112, 23692_19, 24245_2,24682_35,6581_50,7680_205,7691_11,9725_46]
    
    #Inputs will be the path to the folder containing the MA results
    id = sys.argv[1]
    mode = sys.argv[2]
    
    print(f"Processing {id} in {mode} mode...\n")
    
    ma_file = f"/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Code/MA/Results/{id}/MetaAnalysis_{id}_{mode}.tbl.gz"
    mine_ma_file = f"/data/PHURI-Langenberg/people/Mine/SL_MA/01_MA/02_MA/output_rsid_annotated/SeqId_{id}/SL_MA_SeqId_{id}_rsid_annotated.txt.gz"

    if not(os.path.exists(mine_ma_file)):
        print(f"MA results not found for {id}. Skipping.")
        return
    
    
    #Load the MA results into a pandas dataframe
    print(f"Loading MA results for {id}, path: {mine_ma_file}...")
    ma_df = pd.read_csv(mine_ma_file, sep="\t", header=None, comment="#",
                        names=["MarkerName", "CHR", "POS", "dbSNP_ID", "ALLELE1", "ALLELE2", "FREQ1", "FREQSE", "MINFREQ", "MAXFREQ", "EFFECT", "STDERR", "PVALUE", "DIRECTION", "HETISQ", "HETCHISQ", "HETDF", "HETPVAL", "N"]                        
                        )
    
    ma_df['CHR'] = ma_df['CHR'].replace('X', 23).astype(int)
    ma_df['POS'] = ma_df['POS'].astype(int)
    
    ma_df = ma_df.sort_values(by=['CHR', 'POS'])
    
    maxSampleSize = ma_df['N'].max()
    
    
    '''
    Column 1    MarkerName  - chr(chromosome):(position)_(ref)_(alt)
    Column 2    CHR         - chromosome
    Column 3    POS         - position
    Column 4    dbSNP_ID    - dbSNP ID
    Column 5    ALLELE1     - allele 1
    Column 6    ALLELE2     - allele 2
    Column 7    FREQ1       - frequency of allele 1
    Column 8    FREQSE      - standard error of frequency of allele 1
    Column 9    MINFREQ     - minimum frequency of allele 1
    Column 10   MAXFREQ     - maximum frequency of allele 1
    Column 11   EFFECT      - beta
    Column 12   STDERR      - standard error of beta
    Column 13   PVALUE      - p-value
    Column 14   DIRECTION   - direction of effect in each cohort (e.g. +++ means positive effect in all three cohorts)
    Column 15   HETISQ      - I-squared statistic for heterogeneity
    Column 16   HETCHISQ    - chi-squared statistic for heterogeneity
    Column 17   HETDF       - degrees of freedom for heterogeneity
    Column 18   HETPVAL     - p-value for heterogeneity
    Column 19   N           - sample size
    '''
    
    #Identify which variants fall in the same region (same chromosome and within 500kb of each other)
    region = 0
    regions = []
    span = {} 
    prev = None
    
    #Filters
    ma_df = ma_df[ma_df['PVALUE'] <= 5.58e-12]
    #ma_df = ma_df[ma_df['N'] > 0.5 * maxSampleSize]
    #if id not in one_cohort_only:
        #ma_df = ma_df[ma_df['HETDF'] > 0]
    
    if ma_df.empty:
        print(f"No significant variants found for {id} in {mode} mode. Skipping.")
        sentinels = pd.DataFrame(columns=["MarkerName", "CHR", "POS", "dbSNP_ID", "ALLELE1", "ALLELE2", "FREQ1", "FREQSE", "MINFREQ", "MAXFREQ", "EFFECT", "STDERR", "PVALUE", "DIRECTION", "HETISQ", "HETCHISQ", "HETDF", "HETPVAL", "N", "ZSCORE", "REGION", "SPAN", "START", "END"])
    else:
        if ma_df.iloc[0]['CHR'] == 6 and ma_df.iloc[0]['POS'] >= 25500000 and ma_df.iloc[0]['POS'] <= 34000000:
            minimum = min(ma_df.iloc[0]['POS'] - 500000, 25500000)
            maximum = max(ma_df.iloc[0]['POS'] + 500000, 34000000)
        else:
            minimum = max(ma_df.iloc[0]['POS'] - 500000, 0)
            maximum = ma_df.iloc[0]['POS'] + 500000
        
        for _, row in ma_df.iterrows():
            this_minimum = min(row['POS'] - 500000, 25500000) if (row['CHR'] == 6 and row['POS'] >= 25500000 and row['POS'] <= 34000000) else max(row['POS'] - 500000, 0)
            
            if prev is None:
                pass
            elif row['CHR'] != prev['CHR'] or this_minimum > maximum:
                span[region] = (prev['CHR'], minimum, maximum)
                region += 1
                minimum = max(row['POS'] - 500000, 0)
                if row['CHR'] == 6 and row['POS'] >= 25500000 and row['POS'] <= 34000000:
                    minimum = min(row['POS'] - 500000, 25500000)
                    maximum = max(row['POS'] + 500000, 34000000)
                else:
                    maximum = row['POS'] + 500000
            else:
                if row['CHR'] == 6 and row['POS'] >= 25500000 and row['POS'] <= 34000000:
                    minimum = min(row['POS'] - 500000, 25500000, minimum)
                    maximum = max(row['POS'] + 500000, 34000000)
                else:
                    maximum = row['POS'] + 500000
            
            regions.append(region)
            prev = row
        
        span[region] = (prev['CHR'], minimum, maximum)
        print(span)
        
        ma_df['ZSCORE'] = ma_df['EFFECT'] / ma_df['STDERR']
        ma_df['REGION'] = regions
        ma_df['SPAN'] = ma_df['REGION'].map(span)
        ma_df['START'] = ma_df['SPAN'].apply(lambda x: x[1])
        ma_df['END'] = ma_df['SPAN'].apply(lambda x: x[2])
        sentinels = ma_df.loc[ma_df.groupby('REGION')['ZSCORE'].apply(lambda x: x.abs().idxmax())]
    
    os.makedirs("Sentinel-Variant-Analysis/Results", exist_ok=True)
    os.makedirs(f"Sentinel-Variant-Analysis/Results/{id}", exist_ok=True)
    sentinels.to_csv(f"Sentinel-Variant-Analysis/Results/{id}/{mode}_sentinels.csv", index=False)
    
    compare_variants(f"Sentinel-Variant-Analysis/Results/{id}/{mode}_sentinels.csv", f"/data/PHURI-Langenberg/people/Mine/SL_MA/01_MA/04_regional_sentinel_variants/regional_sentinels/{id}_regional_sentinels.txt", id, mode)
    
if __name__ == "__main__":
    main()