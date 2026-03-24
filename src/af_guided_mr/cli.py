import argparse
import sys
import os

from af_guided_mr.pipeline import run_pipeline

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
    run_pipeline(args)

if __name__ == "__main__":
    main()
