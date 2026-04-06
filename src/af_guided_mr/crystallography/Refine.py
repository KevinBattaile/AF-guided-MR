import os
import glob
import subprocess
import re
import shutil
import logging
import concurrent.futures
from queue import Queue, Empty
from dataclasses import dataclass
from typing import Dict, Optional

# Import the lightweight helper functions that remain in utilities
from af_guided_mr.utils.utilities import extract_rfactors
from af_guided_mr.data_management.DataManager import DataManager

@dataclass
class RefinementResult:
    cluster_number: int
    r_work: float 
    r_free: float
    refinement_folder: str
    partial_pdb_path: str
    phaser_output_dir: str
    mode: str
    tfz_score: float
    is_complete: bool = False
    process: Optional[subprocess.Popen] = None 

    @property
    def phaser_output_map(self) -> str:
        """Get the path to the phaser output map file"""
        return glob.glob(os.path.join(self.refinement_folder, '*.mtz'))[0]

def rfactors_from_phenix_refine(pdb_path, data_path, refine_output_root, nproc):
    """Runs phenix.refine and extracts R-factors, handling common errors automatically."""
    # List all existing folders in the refine_output_root directory
    existing_folders = [f for f in os.listdir(refine_output_root) if os.path.isdir(os.path.join(refine_output_root, f))]

    # Filter out folders that match the pattern refine_???
    refine_folders = [f for f in existing_folders if re.match(r'refine_\d{3}', f)]

    # Extract the numeric part of the folder names and find the maximum number
    max_num = 0
    for folder in refine_folders:
        num = int(folder.split('_')[1])
        if num > max_num:
            max_num = num

    # Increment the maximum number by 1 to get the next folder number
    next_num = max_num + 1

    # Format the new folder number to be three digits
    new_folder_name = f"refine_{next_num:03d}"

    # Create the new folder
    refinement_folder = os.path.join(refine_output_root, new_folder_name)
    os.makedirs(refinement_folder, exist_ok=True)

    # check for the existence of the refinement_data.mtz file in the refine_output_root
    if os.path.exists(os.path.join(refine_output_root, "refinement_data.mtz")):
        data_path = os.path.join(refine_output_root, "refinement_data.mtz")

    # Preemptively determine data labels
    selected_data_labels, selected_free_r_label = DataManager.get_mtz_labels(data_path)

    # Initialize the phenix_refine_cmd with base parameters
    phenix_refine_cmd = [
        "phenix.refine",
        pdb_path,
        data_path,
        "strategy=rigid_body+individual_sites+individual_adp",
        "main.number_of_macro_cycles=8",
        f"main.nproc={nproc}",
        "tncs_correction=True",
        "ncs_search.enabled=True",
        "pdb_interpretation.allow_polymer_cross_special_position=True",
        "pdb_interpretation.clash_guard.nonbonded_distance_threshold=None",
        "output.write_eff_file=False",
        "output.write_def_file=False",
        "output.write_geo_file=False",
    ]
    
    if selected_data_labels:
        phenix_refine_cmd.append(f"miller_array.labels.name={selected_data_labels}")
        
    # Inject Free R explicitly to bypass Phenix ambiguity (e.g., Staraniso SA_flags)
    if selected_free_r_label:
        phenix_refine_cmd.append(f"miller_array.labels.name={selected_free_r_label}")

    def run_phenix_refine(cmd):
        formatted_cmd = " ".join(cmd)
        logging.info(f"Running Phenix refine with the following command into {refinement_folder}: {formatted_cmd}")
        # Verify that none of the elements in the command are None
        if any(elem is None for elem in cmd):
            raise ValueError("One or more elements in phenix_refine_cmd are None.")
        process = subprocess.Popen(cmd, cwd=refinement_folder, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate()
        return process.returncode, stdout, stderr, process

    fixed_errors = set()
    max_iterations = 5
    last_process = None

    for iteration in range(max_iterations):
        returncode, stdout, stderr, process = run_phenix_refine(phenix_refine_cmd)
        last_process = process
        if returncode == 0:
            # Success
            break

        re_run = False

        # Check for errors and adjust the command accordingly
        combined_output = stdout + stderr

        if ("Atoms at special positions are within rigid groups" in combined_output) and ('atoms_special_positions' not in fixed_errors):
            logging.warning("Atoms at special positions are within rigid groups. Changing strategy to individual_sites+individual_adp and adding --overwrite.")
            # Modify the strategy argument
            for i, arg in enumerate(phenix_refine_cmd):
                if arg.startswith("strategy="):
                    phenix_refine_cmd[i] = "strategy=individual_sites+individual_adp"
                    break
            else:
                # Strategy argument not found, add it
                phenix_refine_cmd.append("strategy=individual_sites+individual_adp")
            if "--overwrite" not in phenix_refine_cmd:
                phenix_refine_cmd.append("--overwrite")
            fixed_errors.add('atoms_special_positions')
            re_run = True

        if (("R-free flags not compatible" in combined_output) or ("missing flag" in combined_output)) and ('r_free_flags' not in fixed_errors):
            logging.warning("R-free flags not compatible or missing flag. Generating new flags with fraction=0.05 and max_free=500.")
            phenix_refine_cmd.extend([
                "xray_data.r_free_flags.generate=True",
                "xray_data.r_free_flags.fraction=0.05",
                "xray_data.r_free_flags.max_free=500",
            ])
            if "--overwrite" not in phenix_refine_cmd:
                phenix_refine_cmd.append("--overwrite")
            fixed_errors.add('r_free_flags')
            re_run = True

        if re_run:
            continue  # Retry with the adjusted command

        # If no known errors can be fixed, raise an error
        logging.error(f"Phenix refine failed with return code {returncode}.")
        logging.error(stdout)
        logging.error(stderr)
        raise RuntimeError("Phenix refine failed.")

    else:
        # Exceeded maximum iterations
        logging.error("Exceeded maximum number of iterations.")
        raise RuntimeError("Phenix refine failed after maximum attempts.")

    logging.info(f"Phenix refine finished for {os.path.join(os.path.basename(refinement_folder), os.path.basename(pdb_path))}.")
    logging.info(f"Refinement output directory: {os.path.basename(refinement_folder)}")

    r_work, r_free = extract_rfactors(refinement_folder)
    logging.info(f"R_work: {r_work}, R_free: {r_free} for {os.path.basename(pdb_path)}")

    mtz_files = glob.glob(os.path.join(refinement_folder, "*_data.mtz"))
    if mtz_files and os.path.exists(mtz_files[0]):
        shutil.move(mtz_files[0], os.path.join(refine_output_root, "refinement_data.mtz"))

    return r_work, r_free, refinement_folder, last_process


class AsyncRefinementManager:
    def __init__(self, r_free_threshold: float):
        self.refinement_futures: Dict[int, concurrent.futures.Future] = {}
        self.refinement_results: Dict[int, RefinementResult] = {}
        self.executor = concurrent.futures.ThreadPoolExecutor(max_workers=4)
        self.r_free_threshold = r_free_threshold
        self.success_queue = Queue()

    def start_refinement(self, cluster_number: int, partial_pdb_path: str, mtz_path: str, 
                        refine_output_root: str, nproc: int, phaser_output_dir: str,
                        mode: str = "AF_cluster_mode", tfz_score: float = 0.0):
        future = self.executor.submit(
            self._run_refinement,
            cluster_number,
            partial_pdb_path,
            mtz_path,
            refine_output_root,
            nproc,
            phaser_output_dir,
            mode,
            tfz_score
        )
        self.refinement_futures[cluster_number] = future

    def _run_refinement(self, cluster_number: int, partial_pdb_path: str, mtz_path: str, 
                        refine_output_root: str, nproc: int, phaser_output_dir: str,
                        mode: str = "AF_cluster_mode", tfz_score: float = 0.0) -> RefinementResult:
        try:
            # Run refinement locally now, instead of using utilities
            r_work, r_free, refinement_folder, process = rfactors_from_phenix_refine(
                partial_pdb_path, mtz_path, refine_output_root, nproc=nproc
            )

            # Create refinement result if valid R-factors exist
            if r_work > 0 and r_free > 0:
                result = RefinementResult(
                    cluster_number=cluster_number,
                    r_work=r_work,
                    r_free=r_free,
                    refinement_folder=refinement_folder,
                    partial_pdb_path=partial_pdb_path,
                    phaser_output_dir=phaser_output_dir,
                    mode=mode,
                    tfz_score=tfz_score,
                    is_complete=True,
                    process=process
                )

                # Always store valid results
                self.refinement_results[cluster_number] = result

                # Only put in success queue if meets threshold
                if r_free < self.r_free_threshold:
                    logging.info(f"Refinement for cluster {cluster_number} successful with R-free: {r_free:.4f}")
                    self.success_queue.put(result)
                else:
                    logging.info(f"Refinement for cluster {cluster_number} completed with R-free: {r_free:.4f} (above threshold)")

                return result

            return None

        except Exception as e:
            logging.error(f"Refinement for cluster {cluster_number} failed: {e}")
            raise

    def check_completed_refinements(self) -> Optional[RefinementResult]:
        try:
            return self.success_queue.get_nowait()
        except Empty:
            pass

        completed = []
        for cluster_number, future in self.refinement_futures.items():
            if future.done():
                completed.append(cluster_number)
                try:
                    result = future.result()
                    self.refinement_results[cluster_number] = result
                except Exception as e:
                    logging.error(f"Refinement for cluster {cluster_number} failed: {e}")

        for cluster_number in completed:
            del self.refinement_futures[cluster_number]

        return None

    def cleanup(self):
        self.executor.shutdown(wait=False)

    def terminate_all_refinements(self):
        """Terminate all running refinements"""
        for cluster_number, result in self.refinement_results.items():
            if result and result.process and result.process.poll() is None:
                result.process.terminate()
                logging.info(f"Terminated refinement process for cluster {cluster_number}")
        
        for cluster_number, future in self.refinement_futures.items():
            if not future.done():
                future.cancel()
                logging.info(f"Cancelled refinement future for cluster {cluster_number}")
        
        self.refinement_futures.clear()
        self.refinement_results.clear()
        
        while True:
            try:
                self.success_queue.get_nowait()
            except Empty:
                break
