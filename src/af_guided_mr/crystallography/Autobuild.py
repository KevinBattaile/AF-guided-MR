import os
import subprocess
import threading
import logging
import time
import shutil
from af_guided_mr.utils import utilities

class AutobuildManager:
    def __init__(self, job_monitor):
        self.job_monitor = job_monitor

    def run_autobuild(self, data_path, autobuild_input_model, successful_phaser_map, 
                      sequences, solvent_content, nproc, no_waters, skip_autobuild, 
                      output_root, best_result=None):
        
        master_fasta_filename = os.path.join(output_root, "master.fasta")
        with open(master_fasta_filename, "w") as f:
            for sequence_id, sequence in sequences:
                f.write(f">{sequence_id}\n{sequence}\n")

        autobuild_folder = os.path.join(output_root, "autobuild")
        os.makedirs(autobuild_folder, exist_ok=True)
        
        # Save original directory
        original_dir = os.getcwd()
        os.chdir(autobuild_folder)

        phenix_autobuild_cmd = [
            "phenix.autobuild", 
            f"data={data_path}", 
            autobuild_input_model, 
            master_fasta_filename, 
            "rebuild_in_place=False", 
            "include_input_model=True", 
            "n_cycle_rebuild_max=10", 
            f"crystal_info.solvent_fraction={solvent_content}", 
            "use_hl_if_present=False",
            f"nproc={nproc}",
            "thoroughness.ncycle_refine=5",
        ]

        if no_waters:
            phenix_autobuild_cmd.append("refinement.place_waters=False")

        if successful_phaser_map and os.path.exists(successful_phaser_map):
            phenix_autobuild_cmd.append(f"input_map_file={successful_phaser_map}")

        formatted_cmd = ' '.join(phenix_autobuild_cmd)
        logging.info("Formatted phenix_autobuild_cmd for shell execution:")
        logging.info(formatted_cmd)

        if skip_autobuild:
            with open(os.path.join(autobuild_folder, "AUTOBUILD_COMMAND.txt"), "w") as cmd_file:
                cmd_file.write(formatted_cmd)
            logging.info("Skipping autobuild process as --skip_autobuild is specified.")
            
            r_work, r_free, r_factor_folder = None, None, None
            if best_result:
                r_work = best_result.r_work
                r_free = best_result.r_free
                r_factor_folder = best_result.refinement_folder
                
            cc_input_pdb = autobuild_input_model 
            cc_input_map_coeffs = successful_phaser_map
            os.chdir(original_dir)
            return cc_input_pdb, cc_input_map_coeffs, r_work, r_free, r_factor_folder

        if any(elem is None for elem in phenix_autobuild_cmd):
            os.chdir(original_dir)
            raise ValueError("One or more elements in phenix_autobuild_cmd are None.")

        phenix_autobuild_process = subprocess.Popen(phenix_autobuild_cmd)
        main_autobuild_pid = phenix_autobuild_process.pid
        
        # Start tracking this autobuild run
        self.job_monitor.start_autobuild_tracking(main_autobuild_pid, autobuild_folder)
        autobuild_log_path = os.path.join(autobuild_folder, "AutoBuild_run_1_/AutoBuild_run_1_1.log")
        
        # Update subjob PIDs periodically
        def update_tracking():
            while phenix_autobuild_process.poll() is None:
                self.job_monitor.update_subjob_pids()
                time.sleep(60)
                
        tracking_thread = threading.Thread(target=update_tracking, daemon=True)
        tracking_thread.start()

        # Start the monitoring threads
        while not os.path.exists(autobuild_log_path):
            time.sleep(10) 
            
        monitor_autobuild_hanging_thread = threading.Thread(
            target=self.job_monitor.monitor_and_resolve_hangs, 
            args=(autobuild_log_path, phenix_autobuild_process),
            daemon=True
        )
        monitor_autobuild_hanging_thread.start()
        logging.info(f"Monitoring autobuild hanging thread started for {autobuild_log_path}.")

        monitor_autobuild_memory_leaking_thread = threading.Thread(
            target=self.job_monitor.monitor_and_resolve_memory_leaks, 
            args=(autobuild_log_path, phenix_autobuild_process),
            daemon=True
        )
        monitor_autobuild_memory_leaking_thread.start()
        logging.info(f"Monitoring autobuild memory leaking thread started for {autobuild_log_path}.")

        phenix_autobuild_process.wait()
        autobuild_temp_dir = os.path.join(autobuild_folder, "AutoBuild_run_1_/TEMP0")
        time.sleep(60) # wait for the settlement of the TEMP0 folder

        # Terminate only tracked processes for this run
        self.job_monitor.terminate_tracked_processes()

        if os.path.exists(autobuild_temp_dir):
            shutil.rmtree(autobuild_temp_dir, onerror=utilities.remove_readonly)
            
        os.chdir(original_dir)
        logging.info("Autobuild process finished.")
        
        autobuild_working_path = os.path.join(autobuild_folder, "AutoBuild_run_1_")
        cc_input_pdb, cc_input_map_coeffs = utilities.get_autobuild_results_paths(autobuild_working_path)
        logging.info(f"Overall best pdb [Autobuild]: {cc_input_pdb}")
        logging.info(f"Overall best refine map coeffs [Autobuild]: {cc_input_map_coeffs}")
        
        r_factor_folder = autobuild_working_path
        
        return cc_input_pdb, cc_input_map_coeffs, None, None, r_factor_folder
