# AF-guided-MR (AlphaFold-Guided Molecular Replacement)

A robust, fully automated pipeline that bridges the gap between AI structure prediction and X-ray crystallography. AF-guided-MR seamlessly orchestrates ColabFold model generation, Phaser molecular replacement, and Phenix AutoBuild/Refinement to solve crystal structures from sequence and reflection data.

## Features & Recent Upgrades
* **Native FASTA Support:** Simply provide a standard `.fasta` file. The pipeline automatically calculates sequence lengths and Matthews coefficients without requiring cumbersome CSV conversions.
* **Smart MTZ Column Sniffing:** Powered by `gemmi`, the pipeline preemptively scans your MTZ file and automatically selects data arrays. 
* **Automated Fallbacks:** Includes built-in safeguards for headless server environments (headless matplotlib) and automatic fallbacks to CPU-relax if GPU memory limits are hit during ColabFold prediction.
* **Customizable AutoBuild:** Take control of the refinement process with toggles like `--no_waters` to prevent premature solvent placement.

## Prerequisites & Installation

Ensure you have your Conda/Mamba environment set up with the required dependencies, and that your Phenix environment is sourced.

Dependencies include:
* `colabfold_batch`
* `phenix` (Tested on 2.0+ branch)
* `gemmi`
* `biopython`

### 1. Clone the Repository

First, clone the code to your local machine or workstation:
`git clone https://github.com/KevinBattaile/AF-guided-MR.git`

`cd AF-guided-MR`

### 2. Set Up the Conda Environment

Run:
`conda env create -f environment.yml`

`conda activate automatemr`

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
| `--uniprot_id` | Target identifier used for naming output files and tracking. |
| `--nproc` | Number of processors to use for parallelized steps (default: 4). |
| `--no_waters` | **(Optional)** Flag to explicitly disable water placement during the Phenix AutoBuild step. |

## Pipeline Workflow

1. **Structure Prediction:** Uses `colabfold_batch` to generate high-accuracy multimer models from your FASTA sequence.
2. **Domain Parsing:** Evaluates the predicted models and prepares them for MR.
3. **Molecular Replacement:** Calculates the Matthews coefficient and runs `phaser` to place the models into the asymmetric unit.
4. **AutoBuild & Refinement:** Orchestrates `phenix.autobuild` and `phenix.refine` to rebuild missing regions, apply density modifications, and calculate final R-work/R-free and map-model correlation statistics.

## Graphical User Interface (GUI)

AF-Guided MR includes a web-based GUI powered by Gradio. This interface allows you to run the pipeline, adjust advanced settings, and monitor execution logs in real-time directly from your web browser.

### Launching the GUI

1. Activate your Conda environment:
   ```bash
   conda activate afmr
   ```
2. Start the web server from the root of the repository:
   ```bash
   python gui.py
   ```
3. Open your web browser and navigate to the local URL provided in the terminal (typically `http://127.0.0.1:7860`). If running on a remote cluster or server, navigate to `http://<server_ip_address>:7860`.

### GUI Features
* **Drag-and-Drop Uploads:** Easily supply your target `.fasta` and diffraction `.mtz` files.
* **Live Logging:** Monitor standard output, warnings, and underlying subprocesses (like Phenix and ColabFold) in real-time.
* **Process Management:** Use the **Abort Job** button to safely kill underlying threads and subprocesses if you need to cancel a run early.
* **Advanced Settings:** Toggle specific pipeline behaviors (like skipping AutoBuild or forcing AlphaFold clustering) via the dropdown menu.

### Running the GUI as a Background Service

If you are running AF-Guided MR on a remote server or computing cluster, you will likely want to keep the GUI running continuously even after you close your SSH terminal. You can do this using either `nohup` (for a quick, temporary background process) or `systemd` (for a robust, permanent production service).

#### Option A: The Quick Method (`nohup`)
This method runs the server in the background and ignores terminal disconnects. All terminal output is saved to a log file.

1. Activate your Conda environment:
   ```bash
   conda activate afmr
   ```
2. Start the application with `nohup`:
   ```bash
   nohup python gui.py > gradio_gui.log 2>&1 &
   ```
3. To stop the GUI later, find the process ID (PID) and kill it:
   ```bash
   ps -ef | grep app.py
   kill <PID>
   ```

#### Option B: The Production Method (`systemd`)
This method creates a permanent background service that automatically starts when the server boots and restarts itself if it crashes. *(Note: Requires sudo/root privileges).*

1. Create a new service file:
   ```bash
   sudo nano /etc/systemd/system/afmr-gui.service
   ```
2. Paste the following configuration. **Make sure to update** the `User`, `WorkingDirectory`, and the `ExecStart` path to point to the Python executable inside your specific Conda environment:
   ```ini
   [Unit]
   Description=AF-Guided MR Gradio WebUI
   After=network.target

   [Service]
   User=your_username
   WorkingDirectory=/location/of/install
   ExecStart=/path/to/your/conda/envs/afmr/bin/python gui.py
   Restart=always
   RestartSec=3

   [Install]
   WantedBy=multi-user.target
   ```
3. Reload the systemd daemon so it recognizes your new service:
   ```bash
   sudo systemctl daemon-reload
   ```
4. Enable the service (so it starts on boot) and start it immediately:
   ```bash
   sudo systemctl enable afmr-gui
   sudo systemctl start afmr-gui
   ```
5. You can check the status and live logs of the service at any time using:
   ```bash
   sudo systemctl status afmr-gui
   journalctl -u afmr-gui.service -f
   ```
   
## Reference

This is a fork of software described below and 'should' maintain all of the described functionality in the publication and original repository.

Original repository: https://github.com/ww2283/AF-guided-MR/tree/main

    Wang W, Gong Z, Hendrickson WA. AlphaFold-guided molecular replacement for solving challenging crystal structures. 
    Acta Crystallogr D Struct Biol. 2025 Jan 1;81(Pt 1):4-21. doi: 10.1107/S2059798324011999. Epub 2025 Jan 1. 
    PMID: 39711199; PMCID: PMC11740581.
