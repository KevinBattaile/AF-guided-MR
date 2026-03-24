# AF-guided MR

A Python-based tool that automates molecular replacement using protein sequences and x-ray diffraction data, designed especially to handle difficult cases.

This tool uses protein sequences and reduced x-ray diffraction data to automate the process of molecular replacement. It leverages the power of ColabFold for initial structure prediction and Phaser for MR based on various predefined modes. For high-resolution cases better than 3.5 angstroms, it uses AutoBuild to enhance and build the model. For low-resolution cases worse than 3.5 angstroms, it uses phenix.refine to run default refinement cycles for a brief and quick assessment of the molecular replacement correctness. It is specifically designed to handle difficult cases where the predicted structure varies significantly from the final solution.

## Getting Started

### Prerequisites

You will need the following installed on a Linux machine with an Nvidia GPU (CUDA 12) to run the software:

* **Conda** (or **Mamba**, which is highly recommended for faster dependency resolution)
* **PHENIX** (Python-based Hierarchical ENvironment for Integrated Xtallography) ([https://www.phenix-online.org/](https://www.phenix-online.org/))
* **Nvidia GPU Drivers** with CUDA 12 support.

### Installation

The installation process has been streamlined to use a single Conda environment file that resolves all complex dependencies, including JAX, ColabFold, and the necessary C++ extensions.

#### 1. Install PHENIX
Follow the [official PHENIX installation guide](https://www.phenix-online.org/download/) for detailed instructions. The installation process is straightforward and should be completed within a few minutes.

#### 2. Clone the Repository
Download the AF-guided-MR code to your local machine:
```bash
git clone [https://github.com/KevinBattaile/AF-guided-MR](https://github.com/KevinBattaile/AF-guided-MR)
cd AF-guided-MR
```

#### 3. Build the Conda Environment
Use `mamba` (or `conda`) to provision the complete environment. This will automatically install Python 3.10, the CUDA toolkit, ColabFold, and all necessary scientific libraries.
```bash
mamba env create -f environment.yml
```

#### 4. Activate the Environment
```bash
conda activate automatemr
```

*Note on Execution:* The software runs directly from the cloned repository using a wrapper script (`run_afmr.py`). There is no need to run `pip install` on the codebase itself.

## Usage

Ensure your `automatemr` conda environment is active. You can execute the pipeline by pointing your Python interpreter to the `run_afmr.py` script located at the root of the repository.

### Basic Command:
```bash
python run_afmr.py --csv_path your_protein_sequence.csv --mtz_path your_xray_data.mtz --uniprot_id P12345,P23456 --copy_numbers 2:2 --nproc 8
```

### Convenience Setup (Optional):
If you prefer not to type `python /path/to/run_afmr.py` every time, you can make the script executable and add an alias to your `~/.bashrc` or `~/.bash_aliases`:
```bash
chmod +x run_afmr.py
alias afmr='/absolute/path/to/AF-guided-MR/run_afmr.py'
```
You can then run the software from any directory using:
```bash
afmr --csv_path your_protein_sequence.csv --mtz_path your_xray_data.mtz
```

### Inputs
* `--csv_path`: The protein sequence in a CSV format following the ColabFold instructions (`id,sequence`).
* `--mtz_path`: The path to your reduced x-ray diffraction data.
* `--uniprot_id`: (Optional, highly recommended) The UniProt ID(s) for the protein components.
* `--copy_numbers`: (Optional, highly recommended) The copy numbers for the components you wish to search.

Run `python run_afmr.py -h` for a full list of available options and flags.

### How It Works Under the Hood
The script will invoke `AF_cluster.py` to cluster the MSA from ColabFold and sort the resulting structures according to their RMSD to the top-ranked ColabFold model. This ensures that when Phaser runs on these models, it always tests with the order of high RMSD first. 

*Credit: The AF_cluster script draws upon the fabulous work done by the authors of the original AF_Cluster project ([https://github.com/HWaymentSteele/AF_Cluster](https://github.com/HWaymentSteele/AF_Cluster)). A modified version is included here to focus on DBSCAN clustering of the MSA from ColabFold and subsequent model sorting.*

## webUI
A web-based user interface is under development to provide a more user-friendly experience. It will be available soon.

## Note & Benchmarks
This tool has been benchmarked on 372 PDB entries that are deemed hard problems for MR and has shown a 93% success rate at identifying the right solution (R-factors for AutoBuild/refine in a reasonable range), or a 96% success rate at finding a significant solution that requires further human evaluation. It is not a replacement for human evaluation; however, it is a proper tool to greatly speed up the process of MR.

The full tested entry list is available at `resources/tested_cases.txt`. Full test case results are available upon request due to size constraints.

If you experience difficulty with the tool or want to try your case without the full installation process, feel free to contact: `af.guided.mr at gmail.com`.

## License
MIT License.
