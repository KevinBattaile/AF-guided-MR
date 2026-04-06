# utilities.py
import logging
import os
import csv
import shutil
import subprocess
import time
import numpy as np
import glob
import nvidia_smi
import psutil
import gemmi

class CustomFormatter(logging.Formatter):
    def format(self, record):
        if record.levelname == "SUCCESS":
            record.levelname = "\033[1;32mSUCCESS\033[0m"
        elif record.levelname == "FAIL":
            record.levelname = "\033[1;31mFAIL\033[0m"
        elif record.levelname == "ERROR":
            record.levelname = "\033[1;31mERROR\033[0m"
        elif record.levelname == "WARNING":
            record.levelname = "\033[1;33mWARNING\033[0m"
        return super().format(record)

def setup_custom_logger(dir):
    logging.addLevelName(25, "SUCCESS")
    logging.addLevelName(45, "FAIL")

    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    formatter = CustomFormatter(log_format, datefmt='%Y-%m-%d %H:%M:%S')
    
    # Force the directory to exist (does nothing if it already exists)
    os.makedirs(dir, exist_ok=True)

    file_handler = logging.FileHandler(os.path.join(dir, 'automated_structure_solvation.log'))
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    logging.basicConfig(
        level=logging.INFO,
        handlers=[file_handler, stream_handler]
    )

    def success(msg, *args, **kwargs):
        logging.log(25, msg, *args, **kwargs)

    def fail(msg, *args, **kwargs):
        logging.log(45, msg, *args, **kwargs)

    logging.success = success
    logging.fail = fail

def remove_readonly(func, path, excinfo):
    import stat
    os.chmod(path, stat.S_IWRITE)
    func(path)
    
def create_clean_log_copy():
    """
    remove the custom formatting from the automated_structure_solvation.log file and create a cleaned copy
    """
    def remove_ansi_codes(text):
        ansi_escape = re.compile(r'\x1B[@-_][0-?]*[ -/]*[@-~]')
        return ansi_escape.sub('', text)
    
    original_log = 'automated_structure_solvation.log'
    cleaned_log = 'cleaned_automated_structure_solvation.log'

    with open(original_log, 'r') as infile, open(cleaned_log, 'w') as outfile:
        for line in infile:
            clean_line = remove_ansi_codes(line)
            outfile.write(clean_line)

def save_csv_report(output_file, num_sequences, sequence_length, run_time, resolution, tfz_score, successful_phaser_run, successful_phaser_output_dir, reference_model_map_cc, phaser_model_map_cc, r_work, r_free, ):
    with open(output_file, mode='w', newline='') as csvfile:
        fieldnames = ['num_sequences', 'sequence_length', 'run_time', 'resolution', 'tfz_score', 'successful_phaser_run', 'successful_phaser_output_dir', 'reference_model_map_cc', 'phaser_model_map_cc', 'r_work', 'r_free']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({
            'num_sequences': num_sequences,
            'sequence_length': sequence_length,
            'run_time': f'{run_time:.2f}',
            'resolution': f'{resolution:.2f}',
            'tfz_score': tfz_score,
            'successful_phaser_run': successful_phaser_run,
            'successful_phaser_output_dir': successful_phaser_output_dir,
            'reference_model_map_cc': reference_model_map_cc,
            'phaser_model_map_cc': phaser_model_map_cc,
            'r_work': r_work,
            'r_free': r_free
        })


def get_available_cores():
    """This function returns a recommended number of CPU cores to use for a process, e.g. Phaser."""
    total_cores = os.cpu_count()
    
    while True:
        # Get overall CPU usage percentage 
        cpu_percent = psutil.cpu_percent(interval=1)
        
        if cpu_percent > 80:
            logging.info(f"CPU usage at {cpu_percent}%. Waiting 60 seconds before checking again...")
            time.sleep(60)
            continue
            
        used_cores = int(total_cores * (cpu_percent / 100))
        available_cores = total_cores - used_cores
        # Ensure at least 4 cores are used
        recommended_cores = max(4, int(0.25 * available_cores))
        return recommended_cores

def get_cpu_usage(pid):
    cmd = f"ps -p {pid} -o %cpu"
    output = subprocess.check_output(cmd, shell=True, text=True)
    cpu_usage = float(output.splitlines()[1].strip())
    return cpu_usage

def get_available_gpu():
    nvidia_smi.nvmlInit()
    device_count = nvidia_smi.nvmlDeviceGetCount()
    logging.info(f"Total GPUs: {device_count}")

    for device_id in range(device_count):
        handle = nvidia_smi.nvmlDeviceGetHandleByIndex(device_id)
        info = nvidia_smi.nvmlDeviceGetMemoryInfo(handle)
        used_memory_fraction = info.used / info.total
        logging.info(f"GPU {device_id}: {info.used / (1024**2)} MiB used ({used_memory_fraction * 100:.2f}% used)")

        if info.used <= 1000 * 1024 * 1024:
            logging.info(f"GPU {device_id} is available with {info.used / (1024**2)} MiB in use.")
            return device_id

    logging.info("No available GPUs found")
    return None

def wait_for_available_gpu():
    while True:
        logging.info("Checking for available GPUs...")
        device_id = get_available_gpu()
        if device_id is not None:
            logging.info(f"Using GPU {device_id}")
            return device_id
        logging.info("Waiting for 1 minute before checking again...")
        time.sleep(60)

def create_structure_directory(structure_name):
    os.makedirs(structure_name, exist_ok=True)
    return os.path.abspath(structure_name)

def calculate_mean_plddt(pdb_file):
    vals = []

    try:
        with open(pdb_file, "r") as f:
            for line in f:
                if line.startswith("ATOM") and line[13:16].strip() == "CA":  # Check if line is for a C-alpha atom
                    try:
                        b_factor = float(line[60:66].strip())  # Extract B-factor (pLDDT score)
                        vals.append(b_factor)
                    except ValueError:
                        # Handle the case where conversion to float fails
                        continue

        if not vals:
            raise ValueError("No pLDDT scores found in the file.")

        mean_pLDDT = np.mean(vals)
        return mean_pLDDT
    except IOError:
        # Handle file reading errors
        print(f"Error: Unable to read the file {pdb_file}")
        return None
    except Exception as e:
        # Handle other unforeseen errors
        print(f"Error: {str(e)}")
        return None
    

import re



def extract_rfactors(sub_folder_path):

    """
    This extract_rfactors is used to extract R-factors from refinement logs and autobuild logs.
    The input is a folder containing subfolders with refinement logs and/or autobuild logs.
    *For autobuild folders, this sub_folder_path is ./autobuild/AutoBuild_run_?_, 
    *and for refinement folders, this sub_folder_path is ./refine/refine_???.
    """
    def find_refinement_log_file(autobuild_log_path):
        try:
            with open(autobuild_log_path, 'r') as file:
                lines = file.read()
                match = re.search(r"log_refine: (.+\.log_refine)", lines)
                if match:
                    return match.group(1)
                else:
                    logging.debug(f"Could not find 'log_refine' path in {autobuild_log_path}. Falling back to main log.")
                    return None
        except Exception as e:
            logging.error(f"Error finding refinement log file in {autobuild_log_path}: {e}")
            return None

    def extract_rfactors_from_refinement_log(log_file_path):
        try:
            with open(log_file_path, 'r') as file:
                lines = file.readlines()
                last_line = lines[-1]
                r_work, r_free = re.findall(r"R\(work\) = ([\d.]+), R\(free\) = ([\d.]+)", last_line)[0]
                return float(r_work), float(r_free), True
        except Exception as e:
            print(f"Error extracting R-factors from {log_file_path}: {e}")
            return 0.00, 0.00, False

    def extract_rfactors_from_autobuild_log(autobuild_log_path):
        try:
            with open(autobuild_log_path, 'r') as file:
                lines = file.readlines()
                for line in reversed(lines):
                    if "New values of R/Rfree:" in line:
                        r_work, r_free = re.findall(r"R/Rfree:\s*([\d.]+)/\s*([\d.]+)", line)[0]
                        return float(r_work), float(r_free), True
        except Exception as e:
            print(f"Error extracting R-factors from {autobuild_log_path}: {e}")
        return 0.00, 0.00, False

    def extract_rfactors_from_refine_log(refine_folder_path): # folder from phenix.refine
        try:
            refine_log_files = glob.glob(os.path.join(refine_folder_path, '*_refine_001.log'))
            if refine_log_files:
                with open(refine_log_files[0], 'r') as file:
                    lines = file.readlines()
                    for line in reversed(lines):
                        if "Final R-work =" in line:
                            r_work, r_free = re.findall(r"Final R-work = ([\d.]+), R-free = ([\d.]+)", line)[0]
                            return float(r_work), float(r_free), True
        except Exception as e:
            print(f"Error extracting R-factors from refine log: {e}")
        return 0.00, 0.00, False

    r_work, r_free, extracted = 0.00, 0.00, False  # Default values
    """
    end of function definitions; start of main logic for extract_rfactors
    """
    # Check for autobuild log
    autobuild_log_path = os.path.join(sub_folder_path, 'AutoBuild_run_1_1.log')
    if os.path.exists(autobuild_log_path):
        log_refine_path = find_refinement_log_file(autobuild_log_path)
        if log_refine_path:
            r_work, r_free, extracted = extract_rfactors_from_refinement_log(log_refine_path)
        if not extracted:
            r_work, r_free, extracted = extract_rfactors_from_autobuild_log(autobuild_log_path)
    
    # Check for refinement log
    if not extracted:
        refine_folder_path = os.path.join(sub_folder_path)
        r_work, r_free, extracted = extract_rfactors_from_refine_log(refine_folder_path)
    
    return r_work, r_free    


def get_autobuild_results_paths(autobuild_working_path):
    overall_best_pdb = os.path.join(autobuild_working_path, "overall_best.pdb")
    overall_best_refine_map_coeffs = os.path.join(autobuild_working_path, "overall_best_refine_map_coeffs.mtz")
    if not os.path.exists(overall_best_pdb):
        overall_best_pdb = None
    if not os.path.exists(overall_best_refine_map_coeffs):
        overall_best_refine_map_coeffs = None
    return overall_best_pdb, overall_best_refine_map_coeffs

def get_refined_pdb_and_map(refinement_folder_path):
    pdb_file = None
    map_file = None

    # Search for refined pdb file
    pdb_files = glob.glob(os.path.join(refinement_folder_path, '*_refine_001.pdb'))
    if pdb_files:
        pdb_file = pdb_files[0]

    # Search for refined map file
    map_files = glob.glob(os.path.join(refinement_folder_path, '*_refine_001.mtz'))
    if map_files:
        map_file = map_files[0]

    return pdb_file, map_file

def calculate_map_model_correlation(pdb_file, data_file, map_file, solvent_content, output_dir, reference_pdb=None, reference_map=None):
    """
    Calculate the map-model correlation using Phenix, including NCS finding, NCS averaging,
    density modification, and final correlation calculation.

    Parameters:
    pdb_file (str): Path to the PDB file, phaser pdb or autobuild/refinement pdb.
    map_file (str): Path to the map file (MTZ format).
    solvent_content (float): Solvent content percentage (as a decimal, e.g., 0.5 for 50%).
    output_dir (str): Directory where temporary files and results will be stored.

    Returns:
    float: The calculated map-model correlation value.
    """
    correlation_value = None

    try:
        # Ensure the output directory exists
        os.makedirs(output_dir, exist_ok=True)
        # Delete the contents of the output directory
        for file_name in os.listdir(output_dir):
            file_path = os.path.join(output_dir, file_name)
            if os.path.isfile(file_path):
                os.remove(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)

        # Run phenix.get_cc_mtz_pdb to calculate correlation
        try:
            # Determine which PDB to use
            target_pdb = reference_pdb if reference_pdb is not None else pdb_file
            
            # Form the modern Phenix 2.0 command
            get_cc_command = f"phenix.map_correlations {map_file} {target_pdb}"
            
            # Run the command and capture the output directly from the terminal
            result = subprocess.run(get_cc_command, shell=True, check=True, capture_output=True, text=True)
            
            # Scan the terminal output for the CC_mask value
            for line in result.stdout.split('\n'):
                if "CC_mask  :" in line:
                    correlation_value = float(line.split(':')[-1].strip())
                    break
                    
        except subprocess.CalledProcessError as e:
            logging.error(f"Phenix correlation calculation failed: {e.stderr}")
        except Exception as e:
            logging.error(f"An error occurred parsing map-model correlation: {e}")

    # Only attempt the MTZ-to-MTZ correlation if a reference map was actually provided
        if reference_map is not None and (reference_pdb is None or not os.path.exists(cc_log_path)):
            """this is for cc calculation using reference map, in case the reference pdb is not available,
            or for unknown reason the cc.log file is not generated"""
            get_cc_command = f"phenix.get_cc_mtz_mtz mtz_1={map_file} mtz_2={reference_map} output_dir={output_dir}"
            subprocess.run(get_cc_command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            cc_log_path = os.path.join(output_dir, "offset.log")
            if os.path.exists(cc_log_path):
                with open(cc_log_path, "r") as cc_log_file:
                    for line in cc_log_file:
                        if line.startswith("Final CC of maps:"):
                            correlation_value = float(line.split()[-1])
                            break

    except Exception as e:
        # Log the exception if needed
        logging.error(f"An error occurred during map-model/map correlation calculation: {e}")

    return correlation_value

def get_mtz_labels(mtz_path):
    """
    Scans MTZ file to preemptively select the best amplitude/intensity and Free R flag columns.
    Excludes SA_flag to prevent ambiguous array issues in Phenix.
    """
    mtz = gemmi.read_mtz_file(mtz_path)
    columns = [col.label for col in mtz.columns]

    # Priorities for Data labels
    priorities = [
        ("IMEAN", "SIGIMEAN"),
        ("I", "SIGI"),
        ("F", "SIGF"),
        ("FP", "SIGFP"),
    ]

    selected_data_labels = None
    for prio in priorities:
        if prio[0] in columns and prio[1] in columns:
                selected_data_labels = f"{prio[0]},{prio[1]}"
                break

    # Priorities for Free R labels
    free_r_priorities = ["FreeR_flag", "R-free-flags", "FREE"]
    selected_free_r_label = None
    for prio in free_r_priorities:
        if prio in columns and prio != "SA_flag":
            selected_free_r_label = prio
            break

    if selected_data_labels is None:
        logging.warning("Could not find suitable data labels in MTZ file. Refinement may fail.")
    if selected_free_r_label is None:
        logging.warning("Could not find suitable Free R label in MTZ file. Refinement may fail.")

    return selected_data_labels, selected_free_r_label
