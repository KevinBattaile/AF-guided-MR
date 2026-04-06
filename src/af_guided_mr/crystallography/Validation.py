import os
import shutil
import subprocess
import logging

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
