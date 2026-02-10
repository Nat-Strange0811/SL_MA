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
    
    print(data[10:15].collect())
    
    dictionary = {
        "BWHHS_019" : ("SE", "EAF_QTL", "Extracted", "EFFECT_ALLELE", "OTHER_ALLELE", "MarkerName", "INFO"),
    }
    
    
    data = data.with_columns(
        pl.min_horizontal(
                            pl.col(dictionary["BWHHS_019"][1]).cast(pl.Float64),
                            1 - pl.col(dictionary["BWHHS_019"][1]).cast(pl.Float64),
                        ).alias("MAF")).with_columns((2 * pl.col("N").cast(pl.Float64) * pl.col("MAF")).alias("MAC")
    )
                        
    test = data.select([
                    "MAF",
                    "MAC",
                    "N",
                    pl.col(dictionary["BWHHS_019"][1]).cast(pl.Float64).alias("EAF_QTL")
                ]).collect()
    
    res = data.select([
                    (
                        pl.col(dictionary["BWHHS_019"][0]).cast(pl.Float64) < 10
                    ).alias("se"),
                    (
                        pl.min_horizontal(
                            pl.col(dictionary["BWHHS_019"][1]).cast(pl.Float64),
                            1 - pl.col(dictionary["BWHHS_019"][1]).cast(pl.Float64),
                        ) > 0.001
                    ).alias("maf"),
                    (
                        2
                        * pl.col("N").cast(pl.Float64)
                        * pl.min_horizontal(
                            pl.col(dictionary["BWHHS_019"][1]).cast(pl.Float64),
                            1 - pl.col(dictionary["BWHHS_019"][1]).cast(pl.Float64),
                        )
                        > 3
                    ).alias("mac"),
                    (
                        pl.col(dictionary["BWHHS_019"][-1]).cast(pl.Float64) > 0.8
                    ).sum().alias("info")
                ]).collect()
    
    
    res = res.select([
        pl.col("se").sum().alias("se_count"),
        pl.col("maf").sum().alias("maf_count"),
        pl.col("mac").sum().alias("mac_count"),
    ])
    
    
    
    '''
    queries = data.select([
        pl.col(dictionary["BWHHS_019"][3]).cast(pl.Utf8).str.replace("chr", "").str.strip_chars(" \t\n\r").alias("chr"),
        pl.col(dictionary["BWHHS_019"][4]).cast(pl.Int32).alias("pos"),
        pl.col(dictionary["BWHHS_019"][5]).cast(pl.Utf8).str.strip_chars(" \t\n\r").alias("ref"),
        pl.col(dictionary["BWHHS_019"][6]).cast(pl.Utf8).str.strip_chars(" \t\n\r").alias("alt"),
    ]).collect()
    '''
    con = dd.connect(database=':memory:')
    
    #con.register('data', queries)
    con.register('raw_data', test)
    
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
    
    #print(con.table("data").limit(10).df())
    #print(con.table("data").dtypes)
    print()
    
    print(snp.limit(10).df())
    print(snp.dtypes)
    print()
    
    print(con.table("raw_data").limit(10).df())
    
    '''
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
    '''
    
test()