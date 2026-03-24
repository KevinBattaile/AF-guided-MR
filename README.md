# AF-guided-MR (AlphaFold-Guided Molecular Replacement)

A robust, fully automated pipeline that bridges the gap between AI structure prediction and X-ray crystallography. AF-guided-MR seamlessly orchestrates ColabFold model generation, Phaser molecular replacement, and Phenix AutoBuild/Refinement to solve crystal structures from sequence and reflection data.

## Features & Recent Upgrades
* **Native FASTA Support:** Simply provide a standard `.fasta` file. The pipeline automatically calculates sequence lengths and Matthews coefficients without requiring cumbersome CSV conversions.
* **Smart MTZ Column Sniffing:** Powered by `gemmi`, the pipeline preemptively scans your MTZ file and automatically selects the highest-priority data arrays (e.g., autoPROC `F_osf`, `SA_F`, or standard `IMEAN`). No more "Multiple equally suitable arrays" crashes from Phenix!
* **Automated Fallbacks:** Includes built-in safeguards for headless server environments (headless matplotlib) and automatic fallbacks to CPU-relax if GPU memory limits are hit during ColabFold prediction.
* **Customizable AutoBuild:** Take control of the refinement process with toggles like `--no_waters` to prevent premature solvent placement.

## Prerequisites & Installation

Ensure you have your Conda/Mamba environment set up with the required dependencies, and that your Phenix environment is sourced.

Dependencies include:
* `colabfold_batch`
* `phenix` (Tested on 2.0+ branch)
* `gemmi`
* `biopython`

## Usage

Run the pipeline using the main entry point `run_afmr.py`. 

### Basic Command:

    python run_afmr.py \
      --fasta_path /path/to/your/sequence.fasta \
      --mtz_path /path/to/your/data.mtz \
      --uniprot_id TARGET_ID \
      --nproc 8 

### Advanced Command (Skipping Waters):

    python run_afmr.py \
      --fasta_path /path/to/your/sequence.fasta \
      --mtz_path /path/to/your/data.mtz \
      --uniprot_id TARGET_ID \
      --nproc 8 \
      --no_waters

## Command Line Arguments

| Argument | Description |
| :--- | :--- |
| `--fasta_path` | **(Required)** Path to the input `.fasta` file containing your protein sequence(s). |
| `--mtz_path` | **(Required)** Path to the unphased reflection data `.mtz` file. |
| `--uniprot_id` | **(Required)** Target identifier used for naming output files and tracking. |
| `--nproc` | Number of processors to use for parallelized steps (default: 4). |
| `--no_waters` | **(Optional)** Flag to explicitly disable water placement during the Phenix AutoBuild step. |

## Pipeline Workflow

1. **Structure Prediction:** Uses `colabfold_batch` to generate high-accuracy multimer models from your FASTA sequence.
2. **Domain Parsing:** Evaluates the predicted models and prepares them for MR.
3. **Molecular Replacement:** Calculates the Matthews coefficient and runs `phaser` to place the models into the asymmetric unit.
4. **AutoBuild & Refinement:** Orchestrates `phenix.autobuild` and `phenix.refine` to rebuild missing regions, apply density modifications, and calculate final R-work/R-free and map-model correlation statistics.

## Coming Soon
* **Interactive WebUI:** A Gradio-based graphical interface for users who prefer a point-and-click experience over the command line.
