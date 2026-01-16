import polars as pl
import duckdb as dd

def add_key(df):
    return df.with_columns(
        pl.concat_str(
            ["chr", "pos", "ref", "alt"], separator="|"
        )
        .hash()
        .alias("key")
    ).select(["key"])

def test():
    '''
    snp_df = pl.scan_csv(
                        "/data/PHURI-Langenberg/people/SL_MA/02_dbSNP/00-All_filtered.tsv.gz",
                        separator="\t",
                        skip_rows = 56,
                        schema_overrides={"chr": pl.Utf8}
                    ).select(["chr", "pos", "ref", "alt"]
                    )
    '''
    
    data = pl.scan_csv("/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/test.tsv.gz", separator="\t", null_values=["NA"])
    
    dictionary = {
        "CHRIS" : ("SE", "EAF", "Constructed", "CHR", "POS", "EFFECT_ALLELE", "NON_EFFECT_ALLELE", "INFO")
    }
    
    queries = data.select([
        pl.col(dictionary["CHRIS"][3]).cast(pl.Utf8).str.replace("chr", "").str.strip_chars(" \t\n\r").alias("chr"),
        pl.col(dictionary["CHRIS"][4]).cast(pl.Int32).alias("pos"),
        pl.col(dictionary["CHRIS"][5]).cast(pl.Utf8).str.strip_chars(" \t\n\r").alias("ref"),
        pl.col(dictionary["CHRIS"][6]).cast(pl.Utf8).str.strip_chars(" \t\n\r").alias("alt"),
        pl.col("SNPID").alias("rsID")
    ]).collect()
    
    con = dd.connect(database=':memory:')
    
    con.register('data', queries)
    
    snp = con.read_csv(
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
    
    print(con.table("data").limit(10).df())
    print(con.table("data").dtypes)
    print()
    
    print(snp.limit(10).df())
    print(snp.dtypes)
    
    
    # Perform join
    match_count = con.execute("""
    SELECT COUNT(*) AS match_count
    FROM 
        data
    JOIN 
        snp
    ON data.chr = TRIM(snp.chr)
    AND data.pos = snp.pos
    AND (data.ref = ANY(string_split(snp.ref, ',')) OR data.ref = ANY(string_split(snp.alt, ',')))
    AND (data.alt = ANY(string_split(snp.alt, ',')) OR data.alt = ANY(string_split(snp.ref, ',')))
    """).fetchone()[0]

    count = con.execute("""
                        SELECT COUNT(*) AS total_count
                        FROM data
                        JOIN snp
                        ON data.rsID = snp.rsID
                        """).fetchone()[0]
    print(f"Total SNPs with matching rsID: {count}")
    print(f"Number of matching SNPs: {match_count}")
    
    
test()