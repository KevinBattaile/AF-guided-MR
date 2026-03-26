import argparse
import sys
import os
import logging

from af_guided_mr.pipeline import run_pipeline

def validate_args(args):
    """Ensure all inputs are valid before starting the heavy pipeline."""
    # 1. Check if files actually exist
    if not os.path.isfile(args.fasta_path):
        logging.error(f"FATAL: FASTA file not found at {args.fasta_path}")
        sys.exit(1)
        
    if not os.path.isfile(args.mtz_path):
        logging.error(f"FATAL: MTZ file not found at {args.mtz_path}")
        sys.exit(1)
        
    # 2. Check if nproc is a valid number
    if args.nproc is not None and args.nproc < 1:
        logging.error(f"FATAL: Number of processors (--nproc) must be at least 1. Got {args.nproc}")
        sys.exit(1)
        
    # 3. Prevent directory traversal bugs if UniProt IDs are provided
    if args.uniprot_ids and ("/" in args.uniprot_ids or "\\" in args.uniprot_ids):
        logging.error(f"FATAL: uniprot_ids cannot contain slashes. Got {args.uniprot_ids}")
        sys.exit(1)

def parse_args():
    parser = argparse.ArgumentParser(description="Molecular replacement for multi-chain protein complexes using protein sequences and X-ray diffraction data.")
    parser.add_argument("--fasta_path", type=str, required=True, help="Path to the input FASTA file containing protein sequences for each chain.")
    parser.add_argument("--mtz_path", type=str, required=True, help="Path to the input MTZ file containing reduced X-ray diffraction data.")
    parser.add_argument("--uniprot_ids", type=str, required=False, help="Comma-separated UniProt IDs corresponding to the order of sequences in the CSV file.")
    parser.add_argument("--copy_numbers", type=str, required=False, help="Collon-separated copy numbers corresponding to the order of sequences in the CSV file.")
    parser.add_argument("--solvent_content", type=float, required=False, help="Solvent content of the crystal (default: 0.5).")
    parser.add_argument("--output_path", type=str, required=False, help="Base directory for all output files. If not specified, uses the current directory.")
    parser.add_argument("--nproc", type=int, required=False, default=4, help="Number of processors to use (default: 4).")
    parser.add_argument("--no_timeout", action="store_true", help="Disable the 30 minutes timeout for the phaser run.")
    parser.add_argument("--force_af_cluster", action="store_true", help="Force running AF_cluster.py even if the output directory already exists, or for those shorter than 50 residue sequences.")
    parser.add_argument("--skip_af_cluster", action="store_true", help="Skip running AF_cluster.py, while using the models from pre-existing AF_cluster runs.")
    parser.add_argument('--skip_autobuild', action='store_true', help='Skip the autobuild process.')
    parser.add_argument('--no_waters', action='store_true', default=False, help='Explicitly disable water placement in the autobuild run.')
    parser.add_argument("--reference_model", type=str, required=False, help="Developing use: Path to the reference model for optional map-reference correlation.")
    parser.add_argument("--reference_map", type=str, required=False, help="Developing use: Path to the reference map in mtz format for optional map-reference correlation.")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Run the safety checks before doing anything else
    validate_args(args)
    
    # If it passes, run the pipeline
    run_pipeline(args)

if __name__ == "__main__":
    main()
