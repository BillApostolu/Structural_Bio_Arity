import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# --- Paths ---
PHASE4_DIR = os.path.join("Data", "phase4_results")
PLOTS_DIR = os.path.join(PHASE4_DIR, "plots")
os.makedirs(PLOTS_DIR, exist_ok=True)

# File names
files = {
    "Sweep_W": os.path.join(PHASE4_DIR, "Sweep_W.csv"),
    "Sweep_k": os.path.join(PHASE4_DIR, "Sweep_k.csv"),
    "Sweep_r": os.path.join(PHASE4_DIR, "Sweep_r.csv"),
    "Global_Grid": os.path.join(PHASE4_DIR, "Global_Grid_Search_Stats.csv")
}

# ==========================================
# 1. EXCEL EXPORT (Human Readable)
# ==========================================
print("Generating formatted Excel workbook...")
excel_path = os.path.join(PHASE4_DIR, "Master_Results_Formatted.xlsx")

with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
    for name, path in files.items():
        if os.path.exists(path):
            df = pd.read_csv(path)
            
            # Reorder columns slightly to put Hallucinations up front for easy viewing
            if 'Total_Hallucinations' in df.columns:
                cols = df.columns.tolist()
                cols.insert(4, cols.pop(cols.index('Total_Hallucinations')))
                df = df[cols]
                
            df.to_excel(writer, sheet_name=name, index=False)
            
print(f"Saved Excel file to: {excel_path}")

# ==========================================
# 2. DATA VISUALIZATION
# ==========================================
# Set styling for publication-quality thesis graphs
sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.size': 12, 'figure.figsize': (10, 6)})

def plot_sweep(df, x_col, y_cols, labels, title, filename, ylabel):
    """Helper function to plot parameter sweeps."""
    plt.figure()
    for y_col, label in zip(y_cols, labels):
        if y_col in df.columns:
            sns.lineplot(data=df, x=x_col, y=y_col, marker='o', linewidth=2.5, label=label)
            
    plt.title(title, fontsize=16, fontweight='bold')
    plt.xlabel(f"Hyperparameter ({x_col})", fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, filename), dpi=300)
    plt.close()

# --- Plot W Sweep (Sliding Window / Macro-Domains) ---
if os.path.exists(files["Sweep_W"]):
    print("Plotting W Sweep metrics...")
    df_w = pd.read_csv(files["Sweep_W"])
    
    # 1. DTW Error Drop as W increases
    plot_sweep(df_w, x_col="W", 
               y_cols=["DTW_IDP_Med", "DTW_Glob_Med"], 
               labels=["IDPs", "Globular Proteins"],
               title="Impact of Low-Pass Structural Filter (W) on DTW Alignment Error", 
               filename="W_Sweep_DTW_Error.png",
               ylabel="Median DTW Error")
    
    # 2. Entropy changes across physical scales
    plot_sweep(df_w, x_col="W", 
               y_cols=["Entropy_Hp_AF2_Med", "Entropy_Hp_PDB_Med"], 
               labels=["AlphaFold 2", "PDB (Experimental)"],
               title="Topological Complexity (Entropy) Across Biophysical Scales (W)", 
               filename="W_Sweep_Entropy.png",
               ylabel="Persistence Entropy")

# --- Plot k Sweep (Long-Range Contacts) ---
# --- Plot k Sweep (Long-Range Contacts) ---
if os.path.exists(files["Sweep_k"]):
    print("Plotting k Sweep metrics (grouped by units)...")
    df_k = pd.read_csv(files["Sweep_k"])
    
    # Graph 1: Alignment Errors (DTW and EMD)
    plot_sweep(df_k, x_col="k", 
               y_cols=["DTW_IDP_Med", "DTW_Glob_Med", "EMD_IDP_Med", "EMD_Glob_Med"], 
               labels=["DTW (IDPs)", "DTW (Globular)", "EMD (IDPs)", "EMD (Globular)"],
               title="Impact of Sequence Gap (k) on Global Alignment Errors", 
               filename="k_Sweep_1_Alignment_Errors.png",
               ylabel="Alignment Error / Distance Score")
    
    # Graph 2: The IDP Rescue (Signature Distance for IDP vs Globular)
    plot_sweep(df_k, x_col="k", 
               y_cols=["D_sig_2D_IDP_Med", "D_sig_2D_Glob_Med"], 
               labels=["IDPs (Disordered)", "Globular Folds (Rigid)"],
               title="The IDP Rescue: Topological Convergence via Sequence Gap (k)", 
               filename="k_Sweep_2_IDP_Rescue.png",
               ylabel="Median 2D Signature Distance (AF2 vs PDB)")

    # Graph 3: Topological Fragmentation (Counts)
    plot_sweep(df_k, x_col="k", 
               y_cols=["Max_CC_AF2_Med", "Max_CC_PDB_Med", "PLM_AF2_Med", "PLM_PDB_Med"], 
               labels=["Max CC (AF2)", "Max CC (PDB)", "Local Maxima PLM (AF2)", "Local Maxima PLM (PDB)"],
               title="Topological Fragmentation: AF2 vs Experimental Ground Truth", 
               filename="k_Sweep_3_Fragmentation.png",
               ylabel="Count (Features)")

    # Graph 4: Mean Persistence
    plot_sweep(df_k, x_col="k", 
               y_cols=["Mean_Persistence_AF2_Med", "Mean_Persistence_PDB_Med"], 
               labels=["Mean Persistence (AF2)", "Mean Persistence (PDB)"],
               title="Topological Motif Lifespan Across k-Sweep", 
               filename="k_Sweep_4_Persistence.png",
               ylabel="Mean Persistence (Arity Range)")

    # Graph 5: Packing Density (Percentiles)
    plot_sweep(df_k, x_col="k", 
               y_cols=["p_25_AF2_Med", "p_25_PDB_Med", "p_50_AF2_Med", "p_50_PDB_Med", "p_75_AF2_Med", "p_75_PDB_Med"], 
               labels=["AF2 p25", "PDB p25", "AF2 p50", "PDB p50", "AF2 p75", "PDB p75"],
               title="Packing Density Deficit in AlphaFold Predictions", 
               filename="k_Sweep_5_Density_Percentiles.png",
               ylabel="Arity (Neighbor Count)")
    
# --- Plot r Sweep (Radius Density) ---
if os.path.exists(files["Sweep_r"]):
    print("Plotting r Sweep metrics...")
    df_r = pd.read_csv(files["Sweep_r"])
    
    # Bar plot since r only has two discrete values (10, 15)
    plt.figure()
    df_r_melt = df_r.melt(id_vars=['r'], value_vars=['DTW_IDP_Med', 'DTW_Glob_Med'], 
                          var_name='Cluster', value_name='DTW_Error')
    df_r_melt['Cluster'] = df_r_melt['Cluster'].replace({'DTW_IDP_Med': 'IDPs', 'DTW_Glob_Med': 'Globular'})
    
    sns.barplot(data=df_r_melt, x='r', y='DTW_Error', hue='Cluster', palette='muted')
    plt.title("Effect of Radius (r) on Structural Alignment Error", fontsize=16, fontweight='bold')
    plt.xlabel("Spatial Radius (r in Å)", fontsize=14)
    plt.ylabel("Median DTW Error", fontsize=14)
    plt.legend(title="Protein Type")
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, "r_Sweep_DTW_Error.png"), dpi=300)
    plt.close()

print(f"\nAll plots successfully generated and saved to: {PLOTS_DIR}")