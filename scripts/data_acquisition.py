import os
import requests
import time
RAW_DATA_DIR = os.path.join("Data", "raw")
AF2_DIR = os.path.join(RAW_DATA_DIR, "AlphaFold2")
PDB_DIR = os.path.join(RAW_DATA_DIR, "PDB")
TARGET_LIST_FILE = os.path.join("Data", "target_proteins.txt")
TARGET_PAIRS_NEEDED = 20000
REQUEST_TIMEOUT = 30 # seconds per HTTP request, so a stalled server cannot hang the script
os.makedirs(AF2_DIR, exist_ok=True)
os.makedirs(PDB_DIR, exist_ok=True)
if not os.path.exists(TARGET_LIST_FILE):
    with open(TARGET_LIST_FILE, "w") as f:
        f.write("UniProt_ID\tPDB_ID\n")

#Queries the AF2 API
def download_alphafold(uniprot_id):
    filepath = os.path.join(AF2_DIR, f"AF-{uniprot_id}.pdb")
    if os.path.exists(filepath):
        return True

    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        api_response = requests.get(api_url, timeout=REQUEST_TIMEOUT)

        if api_response.status_code == 200:
            data = api_response.json()
            if data and len(data) > 0:
                pdb_url = data[0].get("pdbUrl")
                if pdb_url:
                    pdb_response = requests.get(pdb_url, timeout=REQUEST_TIMEOUT)
                    if pdb_response.status_code == 200:
                        with open(filepath, "w") as f:
                            f.write(pdb_response.text)
                        return True
    except (requests.RequestException, ValueError): # network error, timeout or invalid JSON
        pass
    return False

#queries PDB API
def download_pdb(pdb_id):
    filepath = os.path.join(PDB_DIR, f"{pdb_id.lower()}.pdb")
    if os.path.exists(filepath):
        return True  # Already downloaded

    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    try:
        response = requests.get(url, timeout=REQUEST_TIMEOUT)

        if response.status_code == 200:
            with open(filepath, "w") as f:
                f.write(response.text)
            return True
    except requests.RequestException:
        pass
    return False

#main loop
def build_dataset():
    print(f"Initializing pipeline: Targeting {TARGET_PAIRS_NEEDED} human proteins with BOTH AF2 and PDB structures.")

    # trackers
    proteins_checked = 0

    # load the pairs we already have so a killed script can resume without writing duplicates
    done_ids = set()
    with open(TARGET_LIST_FILE, "r") as f:
        for existing_line in f.readlines()[1:]: # skip header
            existing_parts = existing_line.strip().split('\t')
            if existing_parts and existing_parts[0]:
                done_ids.add(existing_parts[0])
    successful_pairs = len(done_ids)

    if successful_pairs >= TARGET_PAIRS_NEEDED:
        print(f"Already have {successful_pairs} targets. Pipeline complete.")
        return
    if successful_pairs > 0:
        print(f"Resuming: {successful_pairs} pairs already in {TARGET_LIST_FILE} will be skipped.")

    # uniprot url
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": "(taxonomy_id:9606) AND (reviewed:true) AND (database:pdb)",
        "format": "tsv",
        "fields": "accession,xref_pdb",
        "size": 500 } # Max allowed per chunk by UniProt

    # open file
    with open(TARGET_LIST_FILE, "a") as target_file:
        while url and successful_pairs < TARGET_PAIRS_NEEDED:
            print(f"\nFetching a new chunk of candidates from UniProt")
            try:
                if params:
                    response = requests.get(url, params=params, timeout=REQUEST_TIMEOUT)
                    params = None
                else:
                    response = requests.get(url, timeout=REQUEST_TIMEOUT)
            except requests.RequestException as e:
                print(f"Failed to reach UniProt: {e}")
                break
            if response.status_code != 200:
                print(f"Failed to fetch data from UniProt. Error: {response.status_code}")
                break
            raw_text = response.text.strip()
            if not raw_text:
                print("Received empty chunk. Stopping.")
                break
            raw_lines = raw_text.split('\n')
            if len(raw_lines) <= 1:
                print("No more protein data found. Stopping.")
                break
            lines = raw_lines[1:]
            for line in lines:
                if not line.strip():
                    continue
                if successful_pairs >= TARGET_PAIRS_NEEDED:
                    break
                parts = line.split('\t')
                if len(parts) >= 2:
                    uniprot_id = parts[0]
                    if uniprot_id in done_ids:
                        continue # already paired in a previous run
                    pdb_id = parts[1].split(';')[0].strip()
                    if not pdb_id:
                        continue
                    proteins_checked += 1
                    print(f"[{proteins_checked}] Checking {uniprot_id} (PDB: {pdb_id})... ", end="", flush=True)
                    has_af2 = download_alphafold(uniprot_id)
                    time.sleep(0.1)

                    if has_af2:
                        has_pdb = download_pdb(pdb_id)
                        time.sleep(0.1)
                        if has_pdb:
                            target_file.write(f"{uniprot_id}\t{pdb_id}\n")
                            target_file.flush()
                            done_ids.add(uniprot_id)
                            successful_pairs += 1
                            print(f"SUCCESS! (Total Pair Count: {successful_pairs})")
                        else:
                            print(f"Failed: Missing Experimental PDB (Removing orphaned AF2 file)")
                            # ROLLBACK: Delete the AlphaFold file so we don't have unequal counts
                            af2_filepath = os.path.join(AF2_DIR, f"AF-{uniprot_id}.pdb")
                            if os.path.exists(af2_filepath):
                                os.remove(af2_filepath)
                    else:
                        print(f"Failed: Missing AlphaFold Model")

            #'requests' automatically parses standard Link headers safely, ignoring inner commas
            url = response.links.get("next", {}).get("url")

    print(f"\nsuccessfully mapped {successful_pairs} paired targets")
if __name__ == "__main__":
    build_dataset()