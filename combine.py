import os
import pandas as pd

def combine():
    file_path = "/data/PHURI-Langenberg/people/Mine/SL_MA/VIKING/"
    save_path = "/data/PHURI-Langenberg/people/Mine/SL_MA/VIKING/combined/"
    os.makedirs(save_path, exist_ok=True)
    
    txt_file = os.path.join(file_path, "snpID.txt.tmp")
    files = [os.path.join(file_path, f) for f in os.listdir(file_path) if f.endswith(".gz")]
    file_names = [os.path.basename(f) for f in files]
    
    snps = pd.read_csv(txt_file, sep = "\t", quotechar="'")
    
    for i, file in enumerate(files):
        name = file_names[i]
        data = pd.read_csv(file, sep="\t", compression='gzip')
        
        output = pd.concat([snps, data], axis=1)
        
        output.to_csv(os.path.join(save_path, name), sep="\t", index=False, compression="gzip")
        
        return
    

combine()