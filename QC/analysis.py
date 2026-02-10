import matplotlib.pyplot as plt
import polars as pl
import numpy as np
import seaborn as sns
import pandas as pd
import colorcet as cc


def plot_histogram(data, name, info = True):
    
    axes_labels = {
        "SE": "Percentage of Rows with \nStandard Error < 10",
        "INFO": "Percentage of Rows with \nInfo Score > 0.8",
        "Rare Variant": "Percentage of Rows\n with SE >= 10 and either Minor Allele Count <= 3 or \nMinor Allele Frequency <= 0.1%",
        "MAC": "Percentage of Rows with \nMinor Allele Count > 3",
        "MAF": "Percentage of Rows with \nMinor Allele Frequency > 0.1%",
        "rsID": "Percentage of Rows with \nan Existing rsID in build 38"
    }
    
    fig, axes = plt.subplots(3, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    for idx, col in enumerate(data.columns):
        arr =  data[col].to_numpy()
        arr = arr.astype(np.float64)
        arr = arr[~np.isnan(arr)] 
        
        axes[idx].hist(
            arr,
            bins=15,
            alpha=0.7
        )
        axes[idx].set_title(col, fontsize=14, fontweight='bold')
        axes[idx].set_xlabel(axes_labels[col], fontsize=12)
        axes[idx].set_ylabel("Frequency", fontsize=12)
        axes[idx].ticklabel_format(style='plain', axis='y')
        axes[idx].ticklabel_format(style='plain', axis='x')
    
    if not info:
        axes[-1].axis("off")
        
    fig.suptitle(f"Distribution of QC Metrics for {name}", fontsize=18, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/Plots/{name}_QC_Metrics_Distribution.png")
    plt.close()

def processData(group, averages_rows, info = True, combine = False):
    columns = ["SE", "MAC", "MAF", "Rare Variant", "rsID"]
    
    if info:
        columns.append("INFO")
    
    if combine:
        data = pl.concat(
            [df.slice(0, -1) for df in group[1:]],
            how="vertical"
        )[columns]
    else:
        data = group[1].slice(0, -1)[columns]
    
    data = data.with_columns([pl.col(col).round(2) for col in data.columns])
    
    averages_rows.append(
        {
            "Cohort": group[0],
            "SE": data["SE"].mean(),
            "MAC": data["MAC"].mean(),
            "MAF": data["MAF"].mean(),
            "Rare Variant": data["Rare Variant"].mean(),
            "rsID": data["rsID"].mean(),
            "INFO": data["INFO"].mean() if info else np.nan
        }
    )
    
    return data, averages_rows
    

def run():
    
    bwhhs_019 = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/BWHHS_019.csv",
        separator=",",
        schema_overrides={
            "FileID": pl.Utf8
        }
    )
    
    bwhhs_027 = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/BWHHS_027.csv",
        separator=",",
        schema_overrides={
            "FileID": pl.Utf8
        }
    )
    
    bwhhs_controls = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/BWHHS_Controls.csv",
        separator=",",
        schema_overrides={
            "FileID": pl.Utf8
        }
    )
    
    chris = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/CHRIS.csv",
        separator=",",
        schema_overrides={
            "FileID": pl.Utf8
        }
    )
    
    COPD_cases = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/COPD_cases.csv",
        separator=",",
        schema_overrides={
            "FileID": pl.Utf8
        }
    )
    
    COPD_controls = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/COPD_controls.csv",
        separator=",",
        schema_overrides={
            "FileID": pl.Utf8
        }
    )
    
    decode = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/Decode.csv",
        separator=",",
        schema_overrides={
            "FileID": pl.Utf8
        }
    )
    
    EPIC_cases = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/EPIC_Cases.csv",
        separator=",",
        schema_overrides={
            "FileID": pl.Utf8
        }
    )
    
    Epic = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/EPIC.csv",
        separator=",",
        schema_overrides={
            "FileID": pl.Utf8
        }
    )
    
    EPIC_t2d = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/EPIC_T2D.csv",
        separator=",",
        schema_overrides={
            "FileID": pl.Utf8
        }
    )
    
    Fenland_core = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/Fenland_CoreExome.csv",
        separator=",",
        schema_overrides={
            "FileID": pl.Utf8
        }
    )
    
    Fenland_GWAS = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/Fenland_GWAS.csv",
        separator=",",
        schema_overrides={
            "FileID": pl.Utf8
        }
    )   
    
    Fenland_Omics = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/Fenland_Omics.csv",
        separator=",",
        schema_overrides={
            "FileID": pl.Utf8
        }
    )
    
    GenScot = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/GenScot.csv",
        separator=",",
        schema_overrides={
            "FileID": pl.Utf8
        }
    )
    
    Hunt_Controls = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/HUNT_controls.csv",
        separator=",",
        schema_overrides={
            "FileID": pl.Utf8
        }
    )
    
    Hunt_incident = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/HUNT_incident.csv",
        separator=",",
        schema_overrides={
            "FileID": pl.Utf8
        }
    )
    
    Hunt_prev = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/HUNT_prev.csv",
        separator=",",
        schema_overrides={
            "FileID": pl.Utf8
        }
    )
    
    interval = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/INTERVAL.csv",
        separator=",",
        schema_overrides={
            "FileID": pl.Utf8
        }
    )
    
    viking = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/VIKING.csv",
        separator=",",
        schema_overrides={
            "FileID": pl.Utf8
        }
    )
    
    whii  = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/WHII.csv",
        separator=",",
        schema_overrides={
            "FileID": pl.Utf8
        }
    )
    
    all_groups = [
        ["BWHHS_019",bwhhs_019],
        ["BWHHS_027",bwhhs_027],
        ["BWHHS_Controls",bwhhs_controls],
        ["COPD_cases",COPD_cases],
        ["COPD_controls",COPD_controls],
        ["Epic",Epic],
        ["EPIC_t2d",EPIC_t2d],
        ["EPIC_cases",EPIC_cases],
        ["Fenland_core",Fenland_core],
        ["Fenland_GWAS",Fenland_GWAS],
        ["Fenland_Omics",Fenland_Omics],
        ["Hunt_Controls",Hunt_Controls],
        ["Hunt_incident",Hunt_incident],
        ["Hunt_prev",Hunt_prev],
        ["CHRIS",chris],
        ["Decode",decode],
        ["GenScot",GenScot],
        ["INTERVAL",interval],
        ["VIKING",viking],
        ["WHII",whii]
    ]
    
    groups = [
        ["BWHHS Group", bwhhs_019, bwhhs_027, bwhhs_controls],
        ["COPD Group", COPD_cases, COPD_controls],
        ["EPIC Group", Epic, EPIC_t2d, EPIC_cases],
        ["Fenland Group", Fenland_core, Fenland_GWAS, Fenland_Omics],
        ["HUNT Group", Hunt_Controls, Hunt_incident, Hunt_prev]       
    ]
    
    averages_rows = []
    
    for group in all_groups:
        info = False if group[1]["INFO"].is_null().any() else True
        data, averages_rows = processData(group, averages_rows, info)
        plot_histogram(data, group[0], info)

        
    for group in groups:
        info = False if group[1]["INFO"].is_null().any() else True
        data, averages_rows = processData(group, averages_rows, info, combine=True)
        plot_histogram(data, group[0], info)
        
    averages = pl.DataFrame(averages_rows)

    # Boxplot of cohort averages per metric, points colored by cohort
    averages_long = averages.melt(
        id_vars="Cohort",
        value_vars=["SE", "MAC", "MAF", "INFO", "rsID"],
        variable_name="Metric",
        value_name="Percentage",
    ).to_pandas().round(3)
    
    palette = cc.glasbey[:25]
    
    plt.figure(figsize=(8, 5))
    sns.boxplot(data=averages_long, x="Metric", y="Percentage", color="white", width=0.4)
    sns.stripplot(
        data=averages_long,
        x="Metric",
        y="Percentage",
        hue="Cohort",
        palette=palette,
        dodge=True,
        jitter=True,
        alpha=0.8,
        linewidth=0.5,
        edgecolor="black",
    )
    plt.title("Distribution of Cohort Averages per QC Metric", fontsize=18, fontweight='bold')
    plt.legend(
    title="Cohort",
    bbox_to_anchor=(0.5, -0.3),  # position below the axes
    loc="upper center",
    ncol=6,                       # split into 6 columns
    fontsize=8
    )
    plt.tight_layout()
    plt.savefig(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/Plots/Cohort_Averages_Boxplot.png",
        dpi=300,
    )
    plt.close()
    
    for col in averages.columns:
        if col != "Cohort":
            plt.figure(figsize=(12, 6))
            averages_sorted = averages.sort_values("Cohort") if isinstance(averages, pd.DataFrame) else averages.sort("Cohort")
            averages_pd = averages_sorted.to_pandas() if not isinstance(averages_sorted, pd.DataFrame) else averages_sorted
            
            plt.bar(averages_pd["Cohort"], averages_pd[col], color=palette[:len(averages_pd)])
            plt.xlabel("Cohort", fontsize=12)
            plt.ylabel(col, fontsize=12)
            plt.title(f"{col} by Cohort", fontsize=14, fontweight='bold')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(
                f"/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/Plots/{col}_by_Cohort.png",
                dpi=300,
            )
            plt.close()
        
run()