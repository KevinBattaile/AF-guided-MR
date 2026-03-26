import os
import sys
import subprocess
import urllib.request
import urllib.error

def download_test_data(pdb_list_file):
    """Reads a list of PDB IDs, fetches X-ray data via gemmi, and downloads the FASTA."""
    
    # Read and clean up the list of PDB IDs
    if not os.path.isfile(pdb_list_file):
        print(f"Error: Could not find the file '{pdb_list_file}'")
        sys.exit(1)
        
    with open(pdb_list_file, 'r') as f:
        # Strip whitespace, ignore empty lines, and force lowercase for standard naming
        pdb_ids = [line.strip().lower() for line in f if line.strip()]

    print(f"Found {len(pdb_ids)} PDB IDs to process.\n" + "-"*40)

    for pdb_id in pdb_ids:
        print(f"Processing target: {pdb_id.upper()}")
        
        # We will use the PDB ID as the directory name
        out_dir = pdb_id
        os.makedirs(out_dir, exist_ok=True)

        # 1. Fetch the X-ray data using your custom gemmi wrapper
        try:
            print(f"  -> Running fetch_PDB_gemmi...")
            subprocess.run(
                ["fetch_PDB_gemmi", "-d", out_dir, pdb_id], 
                check=True,
                stdout=subprocess.DEVNULL, # Hides the spammy output
                stderr=subprocess.PIPE,
                text=True
            )
        except subprocess.CalledProcessError as e:
            print(f"  [!] Error fetching MTZ data for {pdb_id}: {e.stderr.strip()}")
            continue

        # 2. Fetch the FASTA sequence from the RCSB API
        fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id.upper()}"
        fasta_path = os.path.join(out_dir, f"{pdb_id}.fasta")
        
        try:
            print(f"  -> Downloading FASTA sequence...")
            urllib.request.urlretrieve(fasta_url, fasta_path)
            print(f"  [+] Success! Files saved to ./{out_dir}/")
        except urllib.error.HTTPError as e:
            print(f"  [!] Failed to download FASTA for {pdb_id} (HTTP {e.code})")
        except Exception as e:
            print(f"  [!] Unexpected error downloading FASTA for {pdb_id}: {e}")
            
        print("-" * 40)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python fetch_test_set.py <pdb_list_file.txt>")
        sys.exit(1)
        
    download_test_data(sys.argv[1])
