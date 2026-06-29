import os
import glob
import numpy as np
import pandas as pd
from scipy.stats import wasserstein_distance
from scipy.signal import find_peaks
from numba import jit

CORE_DATA_DIR = os.path.join("Data", "processed_1d_signals_new", "core_structural_data")
PHASE2_SMOOTHED_DIR = os.path.join("Data", "phase2_results", "smoothed_signals")
PHASE3_OUT_DIR = os.path.join("Data", "phase3_results")
os.makedirs(PHASE3_OUT_DIR, exist_ok=True)

DTW_THRESH = 1.0  # after seeing the mean DTW of a whole code run i decided to have it equal to 1
PAE_THRESH = 15.0
PLDDT_THRESH = 70.0 # minimum mean confidence for the domain to be flagged as a confident error
WINDOWS = [1, 5, 11, 15, 25, 31]

@jit(nopython=True) # this justs makes the script faster for DTW runs
def compute_dtw_cost_matrix(arity_pdb, arity_af2, plddt_af2):
    n,m = len(arity_pdb), len(arity_af2)
    cost_matrix = np.zeros((n,m))
    
    for i in range(n):
        for j in range(m):
            weight = plddt_af2[j] / 100.0
            struct_diff = (arity_pdb[i] - arity_af2[j])**2
            cost_matrix[i, j] = weight * struct_diff
    return cost_matrix

@jit(nopython=True)
def compute_dtw_alignment(cost_matrix):
    n,m = cost_matrix.shape
    dtw_matrix = np.full((n + 1,m + 1), np.inf)
    dtw_matrix[0, 0] = 0
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            cost = cost_matrix[i - 1, j - 1]
            dtw_matrix[i, j] = cost + min(dtw_matrix[i - 1, j],dtw_matrix[i, j - 1],dtw_matrix[i - 1, j - 1])
    return dtw_matrix[n,m] / max(n,m)

def get_mean_cross_domain_pae(pae_matrix, min_seq_sep=50):# we want the pae matrix value between residues that are at least 50 indexes apart from each other
    n = pae_matrix.shape[0]
    y,x = np.ogrid[:n, :n]
    mask = abs(x - y) > min_seq_sep
    if np.any(mask):
        return np.mean(pae_matrix[mask])
    return np.nan

def calculate_quaternary_fold_angle(coords, smoothed_arity):
    n = len(coords)
    if n < 50:
        return np.nan 
    #hinge detection
    # it inverts the smoothed arity signal to turn troughs flexible linkers into peaks
    inverted_signal = -smoothed_arity
    
    # here the code finds local minima which are peaks in the inverted signal
    # prominence=2.0 matches phase 2 topological persistence threshold T_nu = 2.0
    # in this way small atomic fluctuations are ignore
    troughs, properties = find_peaks(inverted_signal, prominence=2.0)
    
    if len(troughs) == 0: # if no major structural trough is found it is just a single rigid domain
        return np.nan
        
    # but if multiple troughs exist the primary hinge is selected
    deepest_trough_idx = np.argmax(properties["prominences"])
    linker_idx = troughs[deepest_trough_idx]
    
    # quaternary fold
    hinge_point = coords[linker_idx]
    com_A = np.mean(coords[:linker_idx], axis=0)
    com_B = np.mean(coords[linker_idx+1:], axis=0)
    v_A = com_A - hinge_point
    v_B = com_B - hinge_point
    norm_A = np.linalg.norm(v_A)
    norm_B = np.linalg.norm(v_B)
    if norm_A == 0 or norm_B == 0:
        return np.nan
    
    cosine_angle = np.dot(v_A, v_B) / (norm_A * norm_B)
    cosine_angle = np.clip(cosine_angle, -1.0, 1.0)
    return np.degrees(np.arccos(cosine_angle))

print("Phase 3 execution")

param_folders = glob.glob(os.path.join(PHASE2_SMOOTHED_DIR, "r_*_k_*"))

for folder in param_folders:
    folder_name = os.path.basename(folder)

    for weight_type in ["unweighted"]:
        for w in WINDOWS:
            directory = os.path.join(folder, weight_type, f"W_{w}")                
            print(f"Process DTW alignment for {folder_name} | {weight_type} | W={w}")
            files = glob.glob(os.path.join(directory,"*_AF2.npy"))
            results = []
            
            for af2_file in files:
                filename = os.path.basename(af2_file)
                uniprot_id = filename.replace('_AF2.npy', '')
                pdb_file = os.path.join(directory, f"{uniprot_id}_PDB.npy")
                    
                # smoothed 1D arities of phase 2
                arr_af2 = np.load(af2_file)
                arr_pdb = np.load(pdb_file)
                
                # 3D Data already downloaded
                coords_af2 = np.load(os.path.join(CORE_DATA_DIR, f"{uniprot_id}_AF2_coords.npy"))
                coords_pdb = np.load(os.path.join(CORE_DATA_DIR, f"{uniprot_id}_PDB_coords.npy"))
                plddt_af2 = np.load(os.path.join(CORE_DATA_DIR, f"{uniprot_id}_AF2_pLDDT.npy"))
                pae_matrix = np.load(os.path.join(CORE_DATA_DIR, f"{uniprot_id}_AF2_PAE.npy"))
                    
                # Synchronize masking by removing missing sequence blocks from both arrays for valid DTW comparison
                mask = ~np.isnan(arr_pdb) & ~np.isnan(arr_af2)
                valid_af2 = arr_af2[mask]
                valid_pdb = arr_pdb[mask]
                valid_coords_af2 = coords_af2[mask]
                valid_coords_pdb = coords_pdb[mask]
                valid_plddt_af2 = plddt_af2[mask]
                
                if len(valid_af2) == 0:
                    continue
 
                # 1D EMD
                emd_score = wasserstein_distance(valid_pdb, valid_af2)
                
                # DTW Alignment
                cost_matrix = compute_dtw_cost_matrix(valid_pdb, valid_af2, valid_plddt_af2.astype(np.float64))
                dtw_error = compute_dtw_alignment(cost_matrix)
                
                # hallucination flagging
                cross_pae = get_mean_cross_domain_pae(pae_matrix)
                mean_plddt = np.nanmean(valid_plddt_af2)
                is_hallucination = False
                if not np.isnan(cross_pae) and not np.isnan(mean_plddt):
                    if (dtw_error <= DTW_THRESH) and (cross_pae > PAE_THRESH) and (mean_plddt >= PLDDT_THRESH):
                        is_hallucination = True
                        
                # quaternary validation only if we have hallucinations
                angle_af2, angle_pdb, angle_deviation = np.nan, np.nan, np.nan
                if is_hallucination:
                    angle_af2 = calculate_quaternary_fold_angle(valid_coords_af2, valid_af2)
                    angle_pdb = calculate_quaternary_fold_angle(valid_coords_pdb, valid_pdb)
                    if not np.isnan(angle_af2) and not np.isnan(angle_pdb):
                        angle_deviation = np.abs(angle_af2 - angle_pdb) 

                results.append({
                    "Protein_ID": uniprot_id,
                    "Window_Size": w,
                    "EMD_Score": emd_score,
                    "DTW_Error": dtw_error,
                    "Mean_pLDDT": mean_plddt,
                    "Cross_Domain_PAE": cross_pae,
                    "Hallucination_Flag": is_hallucination,
                    "AF2_Fold_Angle": angle_af2,
                    "PDB_Fold_Angle": angle_pdb,
                    "Angle_Deviation": angle_deviation})
            
            # save the results per folder/weight/window
            if results:
                df = pd.DataFrame(results)
                out_csv = os.path.join(PHASE3_OUT_DIR, f"Phase3_Validation_{folder_name}_{weight_type}_W{w}.csv")
                df.to_csv(out_csv, index=False)
print("\nPhase 3 complete")