import os
import json
import requests
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB.Polypeptide import protein_letters_3to1

RAW_DATA_DIR = os.path.join("Data", "raw")
AF2_DIR = os.path.join(RAW_DATA_DIR, "AlphaFold2")
PDB_DIR = os.path.join(RAW_DATA_DIR, "PDB")
PAE_DIR = os.path.join(RAW_DATA_DIR, "PAE_JSON") # NEW: PAE Directory
TARGET_LIST_FILE = os.path.join("Data", "target_proteins.txt")
PROCESSED_BASE_DIR = os.path.join("Data", "processed_1d_signals_new")

R_VALUES = [5, 8, 10, 12, 15, 20]
K_VALUES = [0, 1, 2, 3, 4, 5, 7, 10]  

# Minimum coverage against the PDB chain length
MIN_ALIGNMENT_COVERAGE = 0.40

def download_pae(uniprot_id):
    """Dynamically fetches the exact PAE JSON URL from the AlphaFold API."""
    pae_path = os.path.join(PAE_DIR, f"AF-{uniprot_id}-pae.json")
    
    # Ask the AF API for the correct URL for this specific protein
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        api_response = requests.get(api_url, timeout=10)
        if api_response.status_code == 200:
            data = api_response.json()
            if data and len(data) > 0:
                # dynamically extract the exact link that handles versions v3, v4, and fragments
                pae_url = data[0].get("paeDocUrl")
                
                if pae_url:
                    pae_response = requests.get(pae_url, timeout=15)
                    if pae_response.status_code == 200:
                        with open(pae_path, 'w') as f:
                            f.write(pae_response.text)
                        return pae_path
    except requests.RequestException:
        pass
    
    return None

def safe_three_to_one(resname):
    resname = resname.strip().upper()
    modifications_map = {'MSE': 'M', 'CSO': 'C', 'SEC': 'U', 'PYL': 'O'}
    if resname in modifications_map:
        return modifications_map[resname]
    return protein_letters_3to1.get(resname, 'X')

def parse_all_chains(filepath):
    """Extracts sequences, CA coordinates, AND B-factors (pLDDT) for all chains."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", filepath)
    model = structure[0]
    
    valid_modifications = {'MSE', 'CSO', 'SEC', 'PYL'}
    chains_data = {}
    
    for chain in model.get_chains():
        seq = []
        coords = []
        b_factors = [] # store pLDDT and B factors
        
        for residue in chain:
            is_standard = residue.id[0] == ' '
            is_valid_mod = residue.get_resname().strip().upper() in valid_modifications
            
            if not (is_standard or is_valid_mod):
                continue
                
            if 'CA' in residue:
                ca_atom = residue['CA']
                seq.append(safe_three_to_one(residue.get_resname()))
                coords.append(ca_atom.get_coord())
                b_factors.append(ca_atom.get_bfactor()) # Extract B-factor
                
        if len(seq) > 0:
            chains_data[chain.id] = ("".join(seq), np.array(coords), np.array(b_factors))
            
    return chains_data

def synchronize_best_chain(af2_sequence, af2_coordinates, af2_plddt, pdb_chains_data):
    """Aligns AF2 against ALL available PDB chains and synchronizes coords + pLDDT."""
    aligner = PairwiseAligner()
    aligner.mode = 'local' 
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10 
    aligner.extend_gap_score = -1 
    
    best_score = -float('inf')
    best_alignment = None
    best_pdb_coords = None
    best_pdb_seq = None
    best_chain_id = None
    
    # chain selection
    for chain_id, (pdb_seq, pdb_coords, pdb_bfactor) in pdb_chains_data.items():
        alignments = aligner.align(af2_sequence, pdb_seq)
        if alignments:
            score = alignments[0].score
            if score > best_score:
                best_score = score
                best_alignment = alignments[0]
                best_pdb_coords = pdb_coords
                best_pdb_seq = pdb_seq
                best_chain_id = chain_id
                
    if best_alignment is None:
        raise ValueError("SW failed to find any sequence overlap across all chains")
        
    matched_residues = sum([end - start for start, end in best_alignment.aligned[1]])
    coverage = matched_residues / len(best_pdb_seq)
    
    if coverage < MIN_ALIGNMENT_COVERAGE:
        raise ValueError(f"Did not find good alignment as best chain {best_chain_id} coverage ({coverage:.1%}) is below {MIN_ALIGNMENT_COVERAGE*100}%.")
    
    pdb_coordinates_aligned = np.full((len(af2_sequence), 3), np.nan)
    
    af2_indices = best_alignment.aligned[0]
    pdb_indices = best_alignment.aligned[1]
    
    for (af2_start, af2_end), (pdb_start, pdb_end) in zip(af2_indices, pdb_indices):
        pdb_coordinates_aligned[af2_start:af2_end] = best_pdb_coords[pdb_start:pdb_end]
        
    return af2_coordinates, pdb_coordinates_aligned, af2_plddt, best_chain_id

def compute_1d_arities(coordis, r_values, k_values):
    n_residues = len(coordis)
    results = {}

    if np.isnan(coordis).all():
        for r in r_values:
            for k in k_values:
                results[(r,k)] = {'weighted': np.full(n_residues, np.nan), 'unweighted': np.full(n_residues, np.nan)}
        return results

    difference = coordis[:, np.newaxis,:] - coordis[np.newaxis, :, :]
    D = np.linalg.norm(difference, axis=-1)
    S = np.abs(np.arange(n_residues)[:, np.newaxis] - np.arange(n_residues)[np.newaxis,:])
    nan_mask = np.isnan(coordis[:, 0])

    for r in r_values:
        for k in k_values:
            valid_mask = (D<=r) & (S>k) & ~np.isnan(D)
            unweighted_arity = np.sum(valid_mask, axis=1).astype(float)
            
            weights = np.exp(-((D / r) ** 2))
            weighted_arity = np.nansum(weights * valid_mask, axis=1) 
            
            unweighted_arity[nan_mask] = np.nan
            weighted_arity[nan_mask] = np.nan
            results[(r,k)] = {'weighted':weighted_arity, 'unweighted':unweighted_arity}          
    return results

# execute
targets = pd.read_csv(TARGET_LIST_FILE, sep="\t") 
successfully_processed = 0

for index, row in targets.iterrows():
    uniprot_id = row['UniProt_ID']
    pdb_id = row['PDB_ID'].lower()
    
    af2_filepath = os.path.join(AF2_DIR, f"AF-{uniprot_id}.pdb")
    pdb_filepath = os.path.join(PDB_DIR, f"pdb{pdb_id}.ent")
    if not os.path.exists(pdb_filepath):
        pdb_filepath = os.path.join(PDB_DIR, f"{pdb_id}.pdb")
        
    if not (os.path.exists(af2_filepath) and os.path.exists(pdb_filepath)):
        continue
        
    print(f"Processing pair: {uniprot_id} | {pdb_id}", end=" ")
    
    try: 
        # Parses AlphaFold where it always grabs chain A automatically
        af2_chains = parse_all_chains(af2_filepath)
        af2_sequence, af2_coordinates, af2_plddt = list(af2_chains.values())[0]
        
        # Parses all available PDB chains
        pdb_chains = parse_all_chains(pdb_filepath)
        
        if len(af2_sequence) == 0 or len(pdb_chains) == 0:
            raise ValueError("Parsed sequences are empty.")

        # finds the best matching chain
        sync_af2_coords, sync_pdb_coords, sync_af2_plddt, matched_chain = synchronize_best_chain(
            af2_sequence, af2_coordinates, af2_plddt, pdb_chains)
        print(f"--> Mapped to Chain {matched_chain}", end=" ")

        #core directory for structural data for next phases
        core_dir = os.path.join(PROCESSED_BASE_DIR, "core_structural_data")
        os.makedirs(core_dir, exist_ok=True)
        
        # PAE cashing and parsing
        pae_npy_path = os.path.join(core_dir, f"{uniprot_id}_AF2_PAE.npy")
        
        # download if not already processed into a numpy array
        if not os.path.exists(pae_npy_path):
            pae_path = download_pae(uniprot_id)
            if pae_path:
                try:
                    with open(pae_path, 'r') as f:
                        pae_data = json.load(f)
                        pae_matrix = np.array(pae_data[0]["predicted_aligned_error"])
                        np.save(pae_npy_path, pae_matrix)
                    
                    # delete the JSON to save disk space
                    os.remove(pae_path)
                except Exception as e:
                    print(f" | PAE parsing failed: {e}", end="")
            else:
                print(" | Warning: PAE download failed", end="")
        
        # Save 1D pLDDT scores and 3D coordinates Nx3
        np.save(os.path.join(core_dir, f"{uniprot_id}_AF2_pLDDT.npy"), sync_af2_plddt)
        np.save(os.path.join(core_dir, f"{uniprot_id}_AF2_coords.npy"), sync_af2_coords)
        np.save(os.path.join(core_dir, f"{uniprot_id}_PDB_coords.npy"), sync_pdb_coords)
        af2_arities = compute_1d_arities(sync_af2_coords, R_VALUES, K_VALUES)
        pdb_arities = compute_1d_arities(sync_pdb_coords, R_VALUES, K_VALUES)
        
        for (r,k) in af2_arities.keys():
            w_dir = os.path.join(PROCESSED_BASE_DIR, f"r_{r}_k_{k}", "weighted")
            u_dir = os.path.join(PROCESSED_BASE_DIR, f"r_{r}_k_{k}", "unweighted")
            os.makedirs(w_dir, exist_ok=True)
            os.makedirs(u_dir, exist_ok=True)
            
            np.save(os.path.join(w_dir, f"{uniprot_id}_AF2.npy"), af2_arities[(r, k)]['weighted'])
            np.save(os.path.join(u_dir, f"{uniprot_id}_AF2.npy"), af2_arities[(r, k)]['unweighted'])
            np.save(os.path.join(w_dir, f"{uniprot_id}_PDB.npy"), pdb_arities[(r, k)]['weighted'])
            np.save(os.path.join(u_dir, f"{uniprot_id}_PDB.npy"), pdb_arities[(r, k)]['unweighted'])
            
        successfully_processed += 1
        
        # record of which chain was the true match for Phase 3 use
        log_file = os.path.join("Data", "mapping_log.csv")
        file_exists = os.path.exists(log_file)
        
        with open(log_file, "a") as log:
            if not file_exists or os.path.getsize(log_file) == 0:
                log.write("UniProt_ID,PDB_ID,Matched_Chain\n")
            log.write(f"{uniprot_id},{pdb_id},{matched_chain}\n")         
        print()# next line  
    except Exception as e:
        print(f"\nFailed: {e}")
print(f"\nPhase 1 complete for {successfully_processed} pairs")