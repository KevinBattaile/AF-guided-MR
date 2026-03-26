import os
import json
import urllib.request
import time

def get_uniprot_from_rcsb(pdb_id):
    """Hits the RCSB REST API to get the UniProt ID for Entity 1."""
    # The API endpoint for polymer entity 1 of a given PDB
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id.upper()}/1"
    
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode('utf-8'))
            
            # Navigate the JSON response to grab the first UniProt ID
            uniprot_ids = data.get('rcsb_polymer_entity_container_identifiers', {}).get('uniprot_ids', [])
            if uniprot_ids:
                return uniprot_ids[0]
    except Exception as e:
        print(f"  [!] Failed to fetch for {pdb_id}: {e}")
    return None

def main():
    base_dir = "."
    print("Scanning for PDB directories...")
    
    # Loop through all items in the current directory
    for item in os.listdir(base_dir):
        item_path = os.path.join(base_dir, item)
        
        # If it's a directory and exactly 4 characters long, assume it's a PDB folder
        if os.path.isdir(item_path) and len(item) == 4:
            print(f"Processing {item}...")
            uniprot_id = get_uniprot_from_rcsb(item)
            
            if uniprot_id:
                print(f"  -> Found UniProt: {uniprot_id}")
                # Save it to a text file inside that specific PDB folder
                out_file = os.path.join(item_path, f"{item}_uniprot.txt")
                with open(out_file, "w") as f:
                    f.write(uniprot_id)
            else:
                print(f"  -> No UniProt ID found for {item}.")
            
            # Be polite to the RCSB servers
            time.sleep(0.2)

if __name__ == "__main__":
    main()
