import os
import pandas as pd
import sys

def compare_variants(folder1, folder2, id, mode):
    
    os.makedirs("Sentinel-Variant-Analysis/Results/Comparison", exist_ok=True)
    
    comparison_file = f"Sentinel-Variant-Analysis/Results/Comparison/comparison_{mode}.csv"
    
    df_nat = pd.read_csv(folder1, sep=",", comment="#")
    
    df_mine = pd.read_csv(folder2, sep="\t", header=None, comment="#",
                         names=["ID", "MarkerName", "CHR", "POS", "dbSNP_ID", "ALLELE1", "ALLELE2", "FREQ1", "FREQSE", "MINFREQ", "MAXFREQ", "EFFECT", "STDERR", "PVALUE", "DIRECTION", "HETISQ", "HETCHISQ", "HETDF", "HETPVAL", "N", "START", "END"]
                         )
    
    nat_sentinels = set(df_nat['MarkerName'])
    mine_sentinels = set(df_mine['MarkerName'])
    shared_sentinels = nat_sentinels.intersection(mine_sentinels)
    
    merged = pd.merge(df_nat, df_mine, on="MarkerName", suffixes=("_nat", "_mine"))
    
    
    span_match = merged[(merged['START_nat'] == merged['START_mine']) & (merged['END_nat'] == merged['END_mine'])]
    
    results = [id, len(nat_sentinels), len(nat_sentinels - mine_sentinels), len(mine_sentinels), len(mine_sentinels - nat_sentinels), len(shared_sentinels), (len(shared_sentinels)/len(nat_sentinels) * 100) if len(nat_sentinels) > 0 else 100.0, len(span_match), (len(span_match)/len(shared_sentinels) * 100) if len(shared_sentinels) > 0 else 100.0]
        
    if not os.path.exists(comparison_file):
        with open(comparison_file, "w") as f:
            f.write("ID,Nat_Sentinels,Nat_Unique,Mine_Sentinels,Mine_Unique,Shared_Sentinels,Correlation,Span_Match,Span_Correlation\n")
            
    with open(comparison_file, "a") as f:
        f.write(",".join(map(str, results)) + "\n")