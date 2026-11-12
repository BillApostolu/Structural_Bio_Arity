import os
import glob
import numpy as np
import pandas as pd
PROCESSED_BASE_DIR = os.path.join("Data", "processed_1d_signals_new")
PHASE2_OUT_DIR = os.path.join("Data", "phase2_results")
SMOOTHED_ARRAYS_DIR = os.path.join(PHASE2_OUT_DIR, "smoothed_signals")

# percentiles for the arity signature
PERCENTILES = [10, 25, 50, 75, 90]
# Persistence threshold for Persistent Local Maxima PLM
T_NU = 2.0 
# Window sizes for hierarchical structural smoothing
WINDOWS = [1, 5, 11, 15, 25, 31] #need to have odd numbers

def calculate_arity_signature(arr):
    """Computes the Arity Signature via the inverse CDF function which are the percentiles and ignores NaNs"""
    if np.isnan(arr).all():
        return {f"p_{p}": np.nan for p in PERCENTILES}
    
    #np.nanpercentile computes the inverse CDF and all percentiles at once while safely ignoring any missing data points
    vals = np.nanpercentile(arr, PERCENTILES)
    return {f"p_{p}": val for p, val in zip(PERCENTILES, vals)}

def smooth_1d_signal(arr, window_size):
    """Applies the Sliding Window filter.
    The average is taken over the resolved residues inside each window, and residues that are NaN
    in the input (unresolved in the PDB structure) stay NaN, so smoothing never fills gaps with invented values."""
    # returns the raw arity series when W = 1
    if window_size < 2:
        return arr
    series = pd.Series(arr) #convert the raw 1D numpy array into a pandas Series

    # .rolling(window=window_size): creates the sliding frame.
    # center=True: ensures the calculated average is placed exactly at the center index of the window,
    # as this keeps the smoothed domains aligned with the actual sequence of the protein.
    # min_periods=1: averages over however many valid residues the window contains (NaNs are skipped),
    # so residues next to a gap or at the chain ends still get a value.
    smoothed = series.rolling(window=window_size, center=True, min_periods=1).mean()

    # Convert back to a fast NumPy array for the rest of the pipeline
    # (copy=True gives a writable array; pandas >= 3 can otherwise return a read-only view)
    smoothed = smoothed.to_numpy(dtype=float, copy=True)

    # min_periods=1 also gives a value to NaN positions within W//2 of a resolved residue,
    # which would extend the PDB profile into unresolved regions, so put the original gaps back
    smoothed[np.isnan(arr)] = np.nan
    return smoothed

def compute_arity_filtration(arr):
    # strip out all missing data because persistent homology requires a continuous sequence
    valid_idx = np.where(~np.isnan(arr))[0]
    if len(valid_idx) == 0:
        return np.nan, np.nan, np.nan, np.nan, np.nan
        
    vals = arr[valid_idx]
    n = len(vals)
    
    # sorts the *indices* of the array by their arity values and this dictates the filtration order, as it processes the lowest arity values first
    order = np.argsort(vals) 
    
    # union find data structures
    # parent: Tracks which component a residue belongs to.
    # birth: Tracks the arity value at which a component was born.
    parent = {i: i for i in range(n)}
    birth = {i: vals[i] for i in range(n)}
    persistences = []
    
    # the find function in union find, which finds the oldest member of the component an index belongs to
    def find(i):
        if parent[i] == i: return i
        parent[i] = find(parent[i])
        return parent[i]
        
    active = np.zeros(n, dtype=bool)
    ncc_history = []
    current_ncc = 0
    
    # filtration main loop
    for idx in order:
        val = vals[idx] # Current arity level
        active[idx] = True # if this residue is now under this level
        current_ncc += 1 # assumes it creates a new isolated island initially
        
        # checks immediate neighbors left and right on the 1D chain
        neighbors = []
        if idx > 0 and active[idx-1]: neighbors.append(find(idx-1))
        if idx < n-1 and active[idx+1]: neighbors.append(find(idx+1))
        
        unique_neighbors = list(set(neighbors))
        
        # merging logic:
        # Case 1: only 1 adjacent neighbor is active and the new residue simply attaches to that existing island as no components merge
        if len(unique_neighbors) == 1:
            parent[idx] = unique_neighbors[0]
            current_ncc -= 1 # Correct the assumption we and do not add a new island
            
        # Case 2: found a merge as both left and right neighbors are active
        elif len(unique_neighbors) == 2:
            c1, c2 = unique_neighbors
            
            # Figure out which island is born at a higher arity
            if birth[c1] > birth[c2]:
                younger, older = c1, c2
            else:
                younger, older = c2, c1
                
            # Calculates Persistence which means how much did the arity level rise between the birth of the younger residue and this merging event
            p = val - birth[younger]
            if p > 0:
                persistences.append(p)
                
            # Execute the union
            parent[younger] = older
            parent[idx] = older
            current_ncc -= 2 # Two residues became one and we initially added 1, so -2
            
        ncc_history.append(current_ncc)
        
    # metrics extraction
    max_ncc = max(ncc_history) if ncc_history else 0
    m_prime = len(persistences) # number of merges
    f_cp = m_prime / n if n > 0 else 0.0 # fraction of merges
    mean_p = np.mean(persistences) if m_prime > 0 else 0.0
    plm = sum(1 for p in persistences if p > T_NU) # PLM
    
    # Calculate entropy of the persistence diagram
    if m_prime > 1:
        total_p = sum(persistences)
        probs = [p / total_p for p in persistences if p > 0]
        entropy = -sum(pr * np.log(pr) for pr in probs) / np.log(m_prime)
    else:
        entropy = 0.0

    return entropy, mean_p, f_cp, max_ncc, plm

# main loop
param_folders = glob.glob(os.path.join(PROCESSED_BASE_DIR, "r_*_k_*"))
for folder in param_folders:
    folder_name = os.path.basename(folder)
    for weight_type in ["unweighted", "weighted"]:
        data_dir = os.path.join(folder, weight_type)
        if not os.path.exists(data_dir):
            continue

        files = glob.glob(os.path.join(data_dir, "*.npy"))
        if not files:
            continue
            
        print(f"\nProcessing phase 2: {folder_name} | {weight_type} | {len(files)} files")
        
        # Iterate over each biophysical window size
        for w in WINDOWS:
            records = []
            # Setup specific directory to save smoothed 1D arrays for Phase 3
            current_smooth_dir = os.path.join(SMOOTHED_ARRAYS_DIR, folder_name, weight_type, f"W_{w}")
            os.makedirs(current_smooth_dir, exist_ok=True)
            
            for file in files:
                filename = os.path.basename(file)
                protein_id, model_type = filename.replace('.npy', '').rsplit('_', 1)
                
                arr = np.load(file)
                
                # Apply smoothing for current window w
                arr_smoothed = smooth_1d_signal(arr, window_size=w)
                
                # Save the smoothed array for Phase 3 (DTW Alignment)
                np.save(os.path.join(current_smooth_dir, filename), arr_smoothed)
                
                # extract arity signature percentiles
                sig = calculate_arity_signature(arr_smoothed)
                
                # filtration metrics
                entropy, mean_p, f_cp, max_ncc, plm = compute_arity_filtration(arr_smoothed)
                
                records.append({
                    "Protein_ID": protein_id,
                    "Model_Type": model_type,
                    "Window_Size": w,
                    "p_10": sig['p_10'],
                    "p_25": sig['p_25'],
                    "p_50": sig['p_50'],
                    "p_75": sig['p_75'],
                    "p_90": sig['p_90'],
                    "Entropy_Hp": entropy,
                    "Mean_Persistence": mean_p,
                    "Fraction_CP": f_cp,
                    "Max_CC": max_ncc,
                    "PLM": plm })

            # Save structural metadata to CSV for this specific window size
            df = pd.DataFrame(records)
            csv_name = f"Signatures_{folder_name}_{weight_type}_W{w}.csv"
            df.to_csv(os.path.join(PHASE2_OUT_DIR, csv_name), index=False)
print("\nPhase 2 complete")