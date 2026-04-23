import polars as pl
import sys
import numpy as np
import os




def main(id):
    outfile = f"MA_Comparison/Results/Pearson_Correlation.csv"
    
    if not os.path.exists(os.path.dirname(outfile)):
        os.makedirs(os.path.dirname(outfile))
    
    try:
        with open(outfile, "r") as f:
            pass
    except FileNotFoundError:
        with open(outfile, "w") as f:
            f.write("ID,Correlation\n")
    
    try:
        nat_data = pl.read_csv(f"MA/Results/{id}/MetaAnalysis_{id}_1.tbl.gz", separator="\t").select(
            [
            pl.col("MarkerName").alias("MarkerName_Nat"), 
             pl.col("Effect").alias("Effect_Nat"),
             
             ]
        )
        mine_data = pl.read_csv(
            f"/data/PHURI-Langenberg/people/Mine/SL_MA/01_MA/02_MA/output_new_comp/SeqId_{id}/SL_MA_SeqId_{id}.tbl.gz", separator="\t").select(
            [
            pl.col("MarkerName").alias("MarkerName_Mine"),
             pl.col("Effect").alias("Effect_Mine")
             ]
        )
    except Exception as e:
        print(f"Error reading files for ID {id}: {e}")
        return
    
    merged_data = nat_data.join(mine_data, left_on="MarkerName_Nat", right_on="MarkerName_Mine", how="full")
    
    errors = merged_data.filter(pl.col("Effect_Mine") != pl.col("Effect_Nat")) # Check for mismatches in Effect values
    
    if errors.height > 0:
        print(f"Errors for ID {id}:")
        print(errors)
        
    result = merged_data.select(
        pl.corr("Effect_Nat", "Effect_Mine").alias("Correlation")
    )['Correlation'][0]
    
    with open(outfile, "a") as f:
        f.write(f"{id},{result}\n")

if __name__ == "__main__":
    id = sys.argv[1]
    
    main(id)