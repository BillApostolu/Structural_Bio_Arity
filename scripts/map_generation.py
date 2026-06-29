import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
PHASE2_OUT_DIR = os.path.join("Data", "phase2_results")
MAPS_OUT_DIR = os.path.join(PHASE2_OUT_DIR, "maps")

def generate_maps():
    csv_files = glob.glob(os.path.join(PHASE2_OUT_DIR, "Signatures_*.csv"))

    for file in csv_files:
        filename = os.path.basename(file)
        # Extract metadata from filename
        parts = filename.replace("Signatures_", "").replace(".csv", "").split("_W")
            
        param_info = parts[0] # like r_5_k_2_weighted
        w_val = parts[1] # like 25
        
        df = pd.read_csv(file)

        df_af2 = df[df['Model_Type'] == 'AF2'].dropna(subset=['p_25','p_50','p_75'])
        df_pdb = df[df['Model_Type'] == 'PDB'].dropna(subset=['p_25','p_50','p_75'])

        fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True, sharey=True)
        GLOBAL_BINS = 30
        GLOBAL_RANGE = [[0, 30], [0, 30]] 
        
        # PDB Density Red
        h_pdb = axes[0].hist2d(df_pdb['p_25'], df_pdb['p_75'], bins=GLOBAL_BINS, range=GLOBAL_RANGE, cmap='Reds', cmin=1)
        fig.colorbar(h_pdb[3], ax=axes[0], label="Number of Proteins")
        axes[0].set_title(f"Experimental PDB (W={w_val})")
        axes[0].set_xlabel("25th Percentile")
        axes[0].set_ylabel("75th Percentile")
        axes[0].grid(True, linestyle='--', alpha=0.3)
        axes[0].set_xlim(0, 30)
        axes[0].set_ylim(0, 30)
        axes[0].plot([0, 30], [0, 30], 'k--', alpha=0.5)
        
        # AF2 Density Blue
        h_af2 = axes[1].hist2d(df_af2['p_25'], df_af2['p_75'], bins=GLOBAL_BINS, range=GLOBAL_RANGE, cmap='Blues', cmin=1)
        fig.colorbar(h_af2[3], ax=axes[1], label="Number of Proteins")
        axes[1].set_title(f"AlphaFold 2 (W={w_val})")
        axes[1].set_xlabel("25th Percentile")
        axes[1].grid(True, linestyle='--', alpha=0.3)
        axes[1].set_xlim(0, 30)
        axes[1].set_ylim(0, 30)
        axes[1].plot([0, 30], [0, 30], 'k--', alpha=0.5) 
        
        # Scatter
        axes[2].scatter(df_pdb['p_25'], df_pdb['p_75'], c='red', alpha=0.3, s=15, label='PDB', edgecolors='none')
        axes[2].scatter(df_af2['p_25'], df_af2['p_75'], c='blue', alpha=0.3, s=15, label='AF2', edgecolors='none')
        axes[2].set_title("Direct Scatter Overlay")
        axes[2].set_xlabel("25th Percentile")
        axes[2].legend()
        axes[2].grid(True, linestyle='--', alpha=0.3)
        axes[2].set_xlim(0, 30)
        axes[2].set_ylim(0, 30)
        axes[2].plot([0, 30], [0, 30], 'k--', alpha=0.5)
        
        plt.suptitle(f"2D Global Arity Maps | {param_info}", fontsize=14)
        plt.tight_layout()
        plt.savefig(os.path.join(MAPS_OUT_DIR, f"Map_2D_{param_info}_W{w_val}.png"), dpi=300)
        plt.clf()
        plt.close('all')

        #3D
        fig = plt.figure(figsize=(18, 7))
        bins_3d = 20

        xyz_pdb = np.vstack([df_pdb['p_25'], df_pdb['p_50'], df_pdb['p_75']]).T
        hist_pdb, edges_pdb = np.histogramdd(xyz_pdb, bins=bins_3d)
        idx_pdb = [np.clip(np.digitize(xyz_pdb[:, i], edges_pdb[i]) - 1, 0, bins_3d - 1) for i in range(3)]
        counts_pdb = hist_pdb[idx_pdb[0], idx_pdb[1], idx_pdb[2]]
        
        ax1 = fig.add_subplot(121, projection='3d')
        sc1 = ax1.scatter(df_pdb['p_25'], df_pdb['p_50'], df_pdb['p_75'], c=counts_pdb, cmap='Reds', alpha=0.8, s=20, edgecolors='none')
        fig.colorbar(sc1, ax=ax1, shrink=0.6, pad=0.1, label="Number of Proteins")
        ax1.set_title(f"PDB 3D Density Map (W={w_val})")
        ax1.set_xlabel('25th (p25)')
        ax1.set_ylabel('50th (p50)')
        ax1.set_zlabel('75th (p75)')

        xyz_af2 = np.vstack([df_af2['p_25'], df_af2['p_50'], df_af2['p_75']]).T
        hist_af2, edges_af2 = np.histogramdd(xyz_af2, bins=bins_3d)      
        idx_af2 = [np.clip(np.digitize(xyz_af2[:, i], edges_af2[i]) - 1, 0, bins_3d - 1) for i in range(3)]
        counts_af2 = hist_af2[idx_af2[0], idx_af2[1], idx_af2[2]]

        ax2 = fig.add_subplot(122, projection='3d')
        sc2 = ax2.scatter(df_af2['p_25'], df_af2['p_50'], df_af2['p_75'],c=counts_af2, cmap='Blues', alpha=0.8, s=20, edgecolors='none', marker='X')
        fig.colorbar(sc2, ax=ax2, shrink=0.6, pad=0.1, label="Number of Proteins")
        ax2.set_title(f"AlphaFold 2 3D Density Map (W={w_val})")
        ax2.set_xlabel('25th (p25)')
        ax2.set_ylabel('50th (p50)')
        ax2.set_zlabel('75th (p75)')

        plt.suptitle(f"Global 3D Arity Binned Maps | {param_info}", fontsize=15)
        plt.tight_layout()
        plt.savefig(os.path.join(MAPS_OUT_DIR, f"Map_3D_{param_info}_W{w_val}.png"), dpi=300)
        plt.clf()
        plt.close('all')

if __name__ == "__main__":
    generate_maps()