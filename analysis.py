import matplotlib.pyplot as plt
import polars as pl
import numpy as np
import seaborn as sns
import pandas as pd

def run():
    
    axes_labels = {
        "SE": "Percentage of Rows with \nStandard Error < 10",
        "INFO": "Percentage of Rows with \nInfo Score > 0.8",
        "MAC": "Percentage of Rows with \nMinor Allele Count > 3",
        "MAF": "Percentage of Rows with \nMinor Allele Frequency > 0.1%",
        "rsID": "Percentage of Rows with \nan Existing rsID in build 38"
    }
    
    bwhhs_019 = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/BWHHS_019.csv",
        separator=","
    )
    
    bwhhs_027 = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/BWHHS_027.csv",
        separator=","
    )
    
    bwhhs_controls = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/BWHHS_Controls.csv",
        separator=","
    )
    
    chris = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/CHRIS.csv",
        separator=","
    )
    
    COPD_cases = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/COPD_cases.csv",
        separator=","
    )
    
    COPD_controls = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/COPD_controls.csv",
        separator=","
    )
    
    decode = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/Decode.csv",
        separator=","
    )
    
    EPIC_cases = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/EPIC_Cases.csv",
        separator=","
    )
    
    Epic = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/EPIC.csv",
        separator=","
    )
    
    EPIC_t2d = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/EPIC_T2D.csv",
        separator=","
    )
    
    Fenland_core = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/Fenland_CoreExome.csv",
        separator=","
    )
    
    Fenland_GWAS = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/Fenland_GWAS.csv",
        separator=","
    )   
    
    Fenland_Omics = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/Fenland_Omics.csv",
        separator=","
    )
    
    GenScot = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/GenScot.csv",
        separator=","
    )
    
    Hunt_Controls = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/HUNT_controls.csv",
        separator=","
    )
    
    Hunt_incident = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/HUNT_incident.csv",
        separator=","
    )
    
    Hunt_prev = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/HUNT_prev.csv",
        separator=","
    )
    
    interval = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/INTERVAL.csv",
        separator=","
    )
    
    viking = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/VIKING.csv",
        separator=","
    )
    
    whii  = pl.read_csv(
        "/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/WHII.csv",
        separator=","
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
    
    individual_cohorts = [
        chris,
        decode,
        GenScot,
        interval,
        viking,
        whii
    ]
    
    averages_rows = []
    
    for group in all_groups:
        info = True
        if group[1]["INFO"].is_null().any():
            data = group[1].slice(0, -1)[["SE", "MAC", "MAF", "rsID"]]
            
            data = data.with_columns([
                pl.col("SE").round(2),
                pl.col("MAC").round(2),
                pl.col("MAF").round(2),
                pl.col("rsID").round(2)
            ])
            
            x = 2
            y = 2
            averages_rows.append(
            {
                "Cohort": group[0],
                "SE": data["SE"].mean(),
                "INFO": np.nan,
                "MAC": data["MAC"].mean(),
                "MAF": data["MAF"].mean(),
                "rsID": data["rsID"].mean()
            }
            )
            
            info = False
        else:
            data = group[1].slice(0, -1)[["SE", "MAC", "MAF", "INFO", "rsID"]]
            
            data = data.with_columns([
                pl.col("SE").round(2),
                pl.col("MAC").round(2),
                pl.col("MAF").round(2),
                pl.col("INFO").round(2),
                pl.col("rsID").round(2)
            ])
            
            x = 2
            y = 3
            averages_rows.append(
            {
                "Cohort": group[0],
                "SE": data["SE"].mean(),
                "INFO": data["INFO"].mean(),
                "MAC": data["MAC"].mean(),
                "MAF": data["MAF"].mean(),
                "rsID": data["rsID"].mean()
            }
        )
        
        fig, axes = plt.subplots(x, y, figsize=(12, 10))
        axes = axes.flatten()
        
        for idx, col in enumerate(data.columns):
            arr =  data[col].to_numpy()
            arr = arr.astype(np.float64)         # force float
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
        
        if info:
            axes[-1].axis("off")
            
        fig.suptitle(f"Distribution of QC Metrics for {group[0]}", fontsize=18, fontweight='bold')
        plt.tight_layout()
        plt.savefig(f"/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/Plots/{group[0]}_QC_Metrics_Distribution.png")
        plt.close()

        
    for group in groups:
        info = True
        if group[1]["INFO"].is_null().any():
            x = 2
            y = 2
            combined_data = pl.concat(group[1:], how="vertical").slice(0, -1)[["SE", "MAC", "MAF", "rsID"]]
            combined_data = combined_data.with_columns([
                pl.col("SE").round(2),
                pl.col("MAC").round(2),
                pl.col("MAF").round(2),
                pl.col("rsID").round(2)
            ])
            averages_rows.append(
            {
                "Cohort": group[0],
                "SE": combined_data["SE"].mean(),
                "INFO": np.nan,
                "MAC": combined_data["MAC"].mean(),
                "MAF": combined_data["MAF"].mean(),
                "rsID": combined_data["rsID"].mean()
            }
            )
            info = False
        else:
            x=2
            y=3
            combined_data = pl.concat(group[1:], how="vertical").slice(0, -1)[["SE", "INFO", "MAC", "MAF", "rsID"]]
            
            combined_data = combined_data.with_columns([
                pl.col("SE").round(2),
                pl.col("MAC").round(2),
                pl.col("MAF").round(2),
                pl.col("INFO").round(2),
                pl.col("rsID").round(2)
            ])
            
            averages_rows.append(
            {
                "Cohort": group[0],
                "SE": combined_data["SE"].mean(),
                "INFO": combined_data["INFO"].mean(),
                "MAC": combined_data["MAC"].mean(),
                "MAF": combined_data["MAF"].mean(),
                "rsID": combined_data["rsID"].mean()
            }
            )
        fig, axes = plt.subplots(x, y, figsize=(12, 10))
        axes = axes.flatten()
        
        for idx, col in enumerate(combined_data.columns):
            arr = combined_data[col].to_numpy()
            arr = arr.astype(np.float64, copy=False)
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
        
        if info:
            axes[-1].axis("off") 
        fig.suptitle(f"Distribution of QC Metrics for {group[0]}", fontsize=18, fontweight='bold')
        plt.tight_layout()
        plt.savefig(f"/data/PHURI-Langenberg/people/SL_MA/01_QC/Nat/Results/Plots/{group[0]}_QC_Metrics_Distribution.png")
        plt.close()
        
    averages = pl.DataFrame(averages_rows)

    # Boxplot of cohort averages per metric, points colored by cohort
    averages_long = averages.melt(
        id_vars="Cohort",
        value_vars=["SE", "MAC", "MAF", "INFO", "rsID"],
        variable_name="Metric",
        value_name="Percentage",
    ).to_pandas().round(3)

    
    palette = sns.color_palette("tab20")
    
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
        
run()