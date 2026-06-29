import os
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import confusion_matrix
import warnings
warnings.filterwarnings("ignore")

PHASE2_DIR = os.path.join("Data", "phase2_results")
PHASE3_DIR = os.path.join("Data", "phase3_results")
PHASE4_OUT_DIR = os.path.join("Data", "phase4_results")

def cliffs_delta(x, y):
    x, y = np.asarray(x).flatten(), np.asarray(y).flatten()
    x, y = x[~np.isnan(x)], y[~np.isnan(y)]
    if len(x) == 0 or len(y) == 0: return np.nan
    mat = np.sign(np.subtract.outer(x, y))
    return np.sum(mat) / (len(x) * len(y))

def load_merged_data(r, k, w, weight="unweighted"):
    folder_name = f"r_{r}_k_{k}"
    sig_file = os.path.join(PHASE2_DIR, f"Signatures_{folder_name}_{weight}_W{w}.csv")
    val_file = os.path.join(PHASE3_DIR, f"Phase3_Validation_{folder_name}_{weight}_W{w}.csv")
        
    df_sig_raw = pd.read_csv(sig_file)
    df_sig = df_sig_raw.pivot(index="Protein_ID", columns="Model_Type", 
                              values=['p_10', 'p_25', 'p_50', 'p_75', 'p_90', 'Entropy_Hp', 'Mean_Persistence', 'Fraction_CP', 'Max_CC', 'PLM'])
    df_sig.columns = [f"{col[0]}_{col[1]}" for col in df_sig.columns]
    
    df_val = pd.read_csv(val_file).set_index("Protein_ID")
    df = df_sig.join(df_val, how='inner')
    
    df['D_sig_2D'] = np.sqrt((df['p_25_PDB'] - df['p_25_AF2'])**2 + (df['p_75_PDB'] - df['p_75_AF2'])**2)
    df['D_sig_3D'] = np.sqrt((df['p_25_PDB'] - df['p_25_AF2'])**2 + 
                             (df['p_50_PDB'] - df['p_50_AF2'])**2 + 
                             (df['p_75_PDB'] - df['p_75_AF2'])**2)
    return df

def generate_baseline_reference_modes():
    # Establish baseline biological families using true experimental PDB data (r=10, k=0, W=1)
    df = load_merged_data(10, 0, 1)
    features = df[['p_25_PDB', 'p_75_PDB']].dropna()
    
    agg_cluster = AgglomerativeClustering(n_clusters=4, linkage='ward')
    raw_agg_labels = agg_cluster.fit_predict(features)
    
    cluster_densities = []
    for c in range(4):
        c_data = features[raw_agg_labels == c]
        cluster_densities.append((c, c_data.mean().sum()))
        
    cluster_densities.sort(key=lambda x: x[1])
    agg_mapping = {old_c: new_c for new_c, (old_c, _) in enumerate(cluster_densities)}
    
    agglom_map = {prot_id: agg_mapping[old_c] for prot_id, old_c in zip(features.index, raw_agg_labels)}
        
    print(f"Assigned {len(agglom_map)} proteins to baseline control groups.")
    return agglom_map

def evaluate_parameter_set(r, k, w, agglom_map):
    df = load_merged_data(r, k, w)
    if df is None or len(df) < 10: return None

    results = {'r': r, 'k': k, 'W': w, 'Total_Proteins': len(df)}
    results['D_sig_2D_Mean'] = df['D_sig_2D'].mean()
    results['D_sig_3D_Mean'] = df['D_sig_3D'].mean()

    # Map to the fixed baseline 4 modes
    df['Fixed_Mode'] = df.index.map(agglom_map)
    for c in range(4):
        c_data = df[df['Fixed_Mode'] == c]
        results[f'DTW_Mode{c}_Med'] = c_data['DTW_Error'].median() if len(c_data) > 0 else np.nan
        results[f'EMD_Mode{c}_Med'] = c_data['EMD_Score'].median() if len(c_data) > 0 else np.nan
        results[f'Mode{c}_Size'] = len(c_data)
    
    # topology between AF2 vs PDB
    topo_metrics = ['Entropy_Hp', 'Mean_Persistence', 'Fraction_CP', 'Max_CC', 'PLM', 
                    'p_10', 'p_25', 'p_50', 'p_75', 'p_90']
    for m in topo_metrics:
        x, y = df[f"{m}_AF2"].dropna(), df[f"{m}_PDB"].dropna()
        common = x.index.intersection(y.index)
        if len(common) > 10:
            _, p_wil = wilcoxon(x[common], y[common], zero_method='zsplit')
            results[f"{m}_AF2_Med"] = x[common].median()
            results[f"{m}_PDB_Med"] = y[common].median()
            results[f"{m}_pval"] = p_wil
            results[f"{m}_Cliffs"] = cliffs_delta(x[common], y[common]) 
        else:
            results[f"{m}_AF2_Med"] = results[f"{m}_PDB_Med"] = results[f"{m}_pval"] = results[f"{m}_Cliffs"] = np.nan

    #hallucinations flagging
    hallucinations = df[df['Hallucination_Flag'] == True]
    results['Total_Hallucinations'] = len(hallucinations)
    if len(hallucinations) > 0:
        results['Mean_Angle_Dev'] = hallucinations['Angle_Deviation'].mean()
        results['Max_Angle_Dev'] = hallucinations['Angle_Deviation'].max()
    else:
        results['Mean_Angle_Dev'] = results['Max_Angle_Dev'] = np.nan
        
    return results

def run_ablation_study(w_val, baseline_k, test_k, r_val=10):
    print(f"Ablation study (W={w_val}, r={r_val}, k: {baseline_k} vs {test_k})")
    df_base = load_merged_data(r_val, baseline_k, w_val)
    df_test = load_merged_data(r_val, test_k, w_val)
    
    if df_base is None or df_test is None: return
    df_merged = df_base.join(df_test, lsuffix='_Base', rsuffix='_Test')
    
    for metric in ['DTW_Error', 'EMD_Score', 'D_sig_2D', 'D_sig_3D']:
        base_data = df_merged[f'{metric}_Base'].dropna()
        test_data = df_merged[f'{metric}_Test'].dropna()
        common = base_data.index.intersection(test_data.index)
        base_data, test_data = base_data[common], test_data[common]
        
        if len(common) > 10:
            stat, p_val = wilcoxon(base_data, test_data, zero_method='zsplit')
            effect = cliffs_delta(test_data, base_data) 
            print(f"  Metric: {metric} | p-val: {p_val:.2e} | Effect: {effect:.3f} (Neg = long range is better)")

def run_confusion_matrix_analysis(r_val=10, k_values=[0, 1, 2, 5], w_values=[1, 5, 11, 15, 25, 31]):
    print("\nStarting PDB vs AF2 independent confusion matrices")
    best_acc = 0.0
    best_params = ""

    for w_val in w_values:
        for k in k_values:
            df = load_merged_data(r_val, k, w_val)
            df_valid = df[['p_25_PDB', 'p_75_PDB', 'p_25_AF2', 'p_75_AF2']].dropna() # ensures 1:1 comparison
            
            pdb_features = df_valid[['p_25_PDB', 'p_75_PDB']]
            af2_features = df_valid[['p_25_AF2', 'p_75_AF2']]
            
            #clusters PDB
            agg_pdb = AgglomerativeClustering(n_clusters=4, linkage='ward')
            pdb_labels = agg_pdb.fit_predict(pdb_features)
            
            pdb_densities = []
            for c in range(4):
                pdb_densities.append((c, pdb_features[pdb_labels == c].mean().sum()))
            pdb_densities.sort(key=lambda x: x[1])
            pdb_map = {old_c: new_c for new_c, (old_c, _) in enumerate(pdb_densities)}
            final_pdb_labels = [pdb_map[l] for l in pdb_labels]
            
            #clusters AF2
            agg_af2 = AgglomerativeClustering(n_clusters=4, linkage='ward')
            af2_labels = agg_af2.fit_predict(af2_features)
            
            af2_densities = []
            for c in range(4):
                af2_densities.append((c, af2_features[af2_labels == c].mean().sum()))
            af2_densities.sort(key=lambda x: x[1])
            af2_map = {old_c: new_c for new_c, (old_c, _) in enumerate(af2_densities)}
            final_af2_labels = [af2_map[l] for l in af2_labels]
            
            cm = confusion_matrix(final_pdb_labels, final_af2_labels, labels=[0, 1, 2, 3])
            agreement_pct = np.trace(cm) / np.sum(cm) * 100
            
            if agreement_pct > best_acc:
                best_acc = agreement_pct
                best_params = f"r={r_val}, k={k}, W={w_val}"
            
            print(f"\nConfusion Matrix (r={r_val}, k={k}, W={w_val})")
            print(f"Overall Accuracy: {agreement_pct:.2f}%")
            
            df_cm = pd.DataFrame(cm, index=["True PDB 0 (IDP)", "True PDB 1", "True PDB 2", "True PDB 3 (Rigid)"],
                                 columns=["Pred AF2 0", "Pred AF2 1", "Pred AF2 2", "Pred AF2 3"])
            print(df_cm.to_string())
            
            cm_out_path = os.path.join(PHASE4_OUT_DIR, f"Confusion_Matrix_r_{r_val}_k_{k}_W_{w_val}.csv")
            df_cm.to_csv(cm_out_path)

    print(f"\nBest Confusion matrix Accuracy: {best_acc:.2f}% at {best_params}")

def print_summary(res, name):
    if res is None: return
    print(f"\n{name} (r={res['r']}, k={res['k']}, W={res['W']})")
    
    if 'DTW_Mode0_Med' in res and not np.isnan(res['DTW_Mode0_Med']):
        print(f"DTW/EMD by Mode | Mode 0: {res['DTW_Mode0_Med']:.2f}/{res['EMD_Mode0_Med']:.2f} "
              f"| Mode 1: {res['DTW_Mode1_Med']:.2f}/{res['EMD_Mode1_Med']:.2f} "
              f"| Mode 2: {res['DTW_Mode2_Med']:.2f}/{res['EMD_Mode2_Med']:.2f} "
              f"| Mode 3: {res['DTW_Mode3_Med']:.2f}/{res['EMD_Mode3_Med']:.2f}")
    print(f"Flagged Proteins: {res.get('Total_Triad_Hallucinations', 0)}")
    print(f"Mean Angle Deviation: {res.get('Mean_Angle_Dev', np.nan):.1f}°")

if __name__ == "__main__":
    print("Starting Phase 4: Statistical Benchmarking & Global Grid Search")

    agglom_map = generate_baseline_reference_modes()

    # W=1
    res_baseline_w1 = evaluate_parameter_set(r=10, k=0, w=1, agglom_map=agglom_map)
    print_summary(res_baseline_w1, "BASELINE (W=1)")

    res_k5_w1 = evaluate_parameter_set(r=10, k=5, w=1, agglom_map=agglom_map)
    print_summary(res_k5_w1, "LONG-RANGE (W=1)")

    # W=5
    res_baseline_w5 = evaluate_parameter_set(r=10, k=0, w=5, agglom_map=agglom_map)
    print_summary(res_baseline_w5, "BASELINE SMOOTHED (W=5)")

    res_k5_w5 = evaluate_parameter_set(r=10, k=5, w=5, agglom_map=agglom_map)
    print_summary(res_k5_w5, "LONG-RANGE SMOOTHED (W=5)")

    print("\n 1D sweeps")
    w_results = [evaluate_parameter_set(10, 0, w, agglom_map) for w in [1, 5, 11, 15, 25, 31]]
    pd.DataFrame([r for r in w_results if r]).to_csv(os.path.join(PHASE4_OUT_DIR, "Sweep_W.csv"), index=False)

    k_results = [evaluate_parameter_set(10, k, 1, agglom_map) for k in [0, 1, 2, 5]]
    pd.DataFrame([r for r in k_results if r]).to_csv(os.path.join(PHASE4_OUT_DIR, "Sweep_k.csv"), index=False)

    r_results = [evaluate_parameter_set(r, 0, 1, agglom_map) for r in [10, 15]]
    pd.DataFrame([r for r in r_results if r]).to_csv(os.path.join(PHASE4_OUT_DIR, "Sweep_r.csv"), index=False)

    print("\nGlobal grid search")
    global_results = []
    for r in [10, 15]:
        for k in [0, 1, 2, 5]:
            for w in [1, 5, 11, 15, 25, 31]:
                res = evaluate_parameter_set(r, k, w, agglom_map)
                if res: global_results.append(res)

    if global_results:
        df_global = pd.DataFrame(global_results)
        df_global = df_global.sort_values(by=['r', 'k', 'W'])
        out_path = os.path.join(PHASE4_OUT_DIR, "Global_Grid_Search_Stats.csv")
        df_global.to_csv(out_path, index=False)
        best_idx = df_global['D_sig_2D_Mean'].idxmin()
        best_row = df_global.loc[best_idx]
        print("\nOptimal Hyperparameters (Based on lowest mean D_sig_2D between AF2 and PDB):")
        print(f"Radius (r): {best_row['r']} Å")
        print(f"Sequence Gap (k): {best_row['k']} residues")
        print(f"Window (W): {best_row['W']}")
        print(f"Achieved D_sig_2D_Mean: {best_row['D_sig_2D_Mean']:.3f}\n")
    print("Ablation sweep")
    for w_val in [1, 5, 15, 25]:
        for test_k in [1, 2, 5]:
            run_ablation_study(w_val=w_val, baseline_k=0, test_k=test_k, r_val=10)
    run_confusion_matrix_analysis(r_val=10, k_values=[0, 1, 2, 5], w_values=[1, 5, 11, 15, 25, 31])
    print("\nPhase 4 Analysis Complete.")