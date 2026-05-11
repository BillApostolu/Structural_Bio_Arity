import os
import requests
import time

# --- CONFIGURATION ---
RAW_DATA_DIR = os.path.join("Data", "raw")
AF2_DIR = os.path.join(RAW_DATA_DIR, "AlphaFold2")
PDB_DIR = os.path.join(RAW_DATA_DIR, "PDB")
TARGET_LIST_FILE = os.path.join("Data", "target_proteins.txt")

# How many successful PAIRS (AF2 + PDB) do you want in total?
TARGET_PAIRS_NEEDED = 10000 

# Create directories if they don't exist
os.makedirs(AF2_DIR, exist_ok=True)
os.makedirs(PDB_DIR, exist_ok=True)

# Ensure the target list file exists (we will append to it)
if not os.path.exists(TARGET_LIST_FILE):
    with open(TARGET_LIST_FILE, "w") as f:
        f.write("UniProt_ID\tPDB_ID\n") # Header

# --- HELPER FUNCTIONS ---

def download_alphafold(uniprot_id):
    """
    Queries the AF2 API. Returns True if downloaded successfully, False otherwise.
    """
    filepath = os.path.join(AF2_DIR, f"AF-{uniprot_id}.pdb")
    if os.path.exists(filepath):
        return True  # Already downloaded
        
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    api_response = requests.get(api_url)
    
    if api_response.status_code == 200:
        data = api_response.json()
        if data and len(data) > 0:
            pdb_url = data[0].get("pdbUrl")
            if pdb_url:
                pdb_response = requests.get(pdb_url)
                if pdb_response.status_code == 200:
                    with open(filepath, "w") as f:
                        f.write(pdb_response.text)
                    return True
    return False

def download_pdb(pdb_id):
    """
    Downloads experimental PDB. Returns True if successful, False otherwise.
    """
    filepath = os.path.join(PDB_DIR, f"{pdb_id.lower()}.pdb")
    if os.path.exists(filepath):
        return True  # Already downloaded
        
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    response = requests.get(url)
    
    if response.status_code == 200:
        with open(filepath, "w") as f:
            f.write(response.text)
        return True
    return False

# --- MAIN ACQUISITION LOOP ---

def build_dataset():
    print(f"Initializing pipeline: Targeting {TARGET_PAIRS_NEEDED} human proteins with BOTH AF2 and PDB structures.")
    
    # Trackers
    successful_pairs = 0
    proteins_checked = 0
    
    # Check how many we already have so we can resume a killed script
    with open(TARGET_LIST_FILE, "r") as f:
        existing_lines = f.readlines()[1:] # Skip header
        successful_pairs = len(existing_lines)
        
    if successful_pairs >= TARGET_PAIRS_NEEDED:
        print(f"Already have {successful_pairs} targets. Pipeline complete.")
        return

    # UniProt Initial URL and Params
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": "(taxonomy_id:9606) AND (reviewed:true) AND (database:pdb)",
        "format": "tsv",
        "fields": "accession,xref_pdb",
        "size": 500 # Max allowed per chunk by UniProt
    }

    # Open file in append mode
    with open(TARGET_LIST_FILE, "a") as target_file:
        
        while url and successful_pairs < TARGET_PAIRS_NEEDED:
            print(f"\n--- Fetching a new chunk of candidates from UniProt ---")
            
            # If it's the first request, use params. If it's a paginated URL, params are built-in.
            if params:
                response = requests.get(url, params=params)
                params = None # Clear params for subsequent paginated requests
            else:
                response = requests.get(url)
                
            if response.status_code != 200:
                print(f"Failed to fetch data from UniProt. Error: {response.status_code}")
                break

            # Safely parse the TSV lines to avoid crashing on page turns
            raw_text = response.text.strip()
            if not raw_text:
                print("Received empty chunk. Stopping.")
                break
                
            raw_lines = raw_text.split('\n')
            
            # If the chunk only contains a header and no data, stop
            if len(raw_lines) <= 1:
                print("No more protein data found. Stopping.")
                break
                
            lines = raw_lines[1:] # Safely skip the header row
            
            for line in lines:
                if not line.strip(): # Skip any weird empty lines
                    continue 
                    
                if successful_pairs >= TARGET_PAIRS_NEEDED:
                    break # Stop exactly at the target amount
                    
                parts = line.split('\t')
                if len(parts) >= 2:
                    uniprot_id = parts[0]
                    # Get the first PDB ID in the list
                    pdb_id = parts[1].split(';')[0].strip() 
                    
                    if not pdb_id:
                        continue
                        
                    proteins_checked += 1
                    print(f"[{proteins_checked}] Checking {uniprot_id} (PDB: {pdb_id})... ", end="", flush=True)

                    # 1. Try downloading AlphaFold first (it's the most common bottleneck)
                    has_af2 = download_alphafold(uniprot_id)
                    time.sleep(0.1) # Polite delay
                    
                    if has_af2:
                        # 2. If AF2 exists, try PDB
                        has_pdb = download_pdb(pdb_id)
                        time.sleep(0.1) # Polite delay
                        
                        if has_pdb:
                            # 3. BOTH exist! Save to dataset.
                            target_file.write(f"{uniprot_id}\t{pdb_id}\n")
                            target_file.flush() # Force write to disk immediately
                            successful_pairs += 1
                            print(f"SUCCESS! (Total Pair Count: {successful_pairs})")
                        else:
                            print(f"Failed: Missing Experimental PDB")
                    else:
                        print(f"Failed: Missing AlphaFold Model")

            # Handle UniProt Pagination for the next chunk of 500 (UNTOUCHED)
            url = None # Reset URL
            if "Link" in response.headers:
                links = response.headers["Link"].split(',')
                for link in links:
                    if 'rel="next"' in link:
                        # Extract the URL from inside the < > brackets
                        url = link[link.find("<")+1:link.find(">")]
                        break

    print(f"\nPhase 1 Data Acquisition Complete! Successfully mapped {successful_pairs} paired targets.")

if __name__ == "__main__":
    build_dataset()