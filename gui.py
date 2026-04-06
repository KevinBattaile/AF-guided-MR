import gradio as gr
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "src")))
import traceback
import time
import logging
import multiprocessing
from af_guided_mr.pipeline import run_pipeline

# Global variable to track the running process so we can cleanly kill it
current_process = None

def pipeline_worker(args_dict, queue):
    """Runs the pipeline in an isolated process to capture ALL terminal output."""
    
    # 1. Rebuild the arguments object for the pipeline
    class PipelineArgs:
        pass
    args = PipelineArgs()
    for k, v in args_dict.items():
        setattr(args, k, v)
        
    # 2. Intercept raw print() and subprocess outputs (stdout/stderr)
    class QueueStream:
        def __init__(self, q, prefix=""):
            self.q = q
            self.prefix = prefix
        def write(self, msg):
            if msg.strip():
                self.q.put(f"{self.prefix}{msg.strip()}\n")
        def flush(self):
            pass
            
    sys.stdout = QueueStream(queue)
    sys.stderr = QueueStream(queue, prefix="ERROR: ")
    
    # 3. Intercept standard Python logging
    class QueueHandler(logging.Handler):
        def emit(self, record):
            queue.put(self.format(record) + '\n')
            
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.handlers = [] # Clear old handlers to prevent duplicate terminal prints
    qh = QueueHandler()
    qh.setFormatter(logging.Formatter('%(asctime)s - %(message)s', datefmt='%H:%M:%S'))
    logger.addHandler(qh)

    # 4. Run the pipeline safely
    try:
        run_pipeline(args)
        queue.put("\n✅ PIPELINE COMPLETED SUCCESSFULLY!")
    except Exception as e:
        queue.put(f"\n❌ PIPELINE FAILED:\n{traceback.format_exc()}")


def execute_pipeline(fasta_file, mtz_file, output_dir, uniprot_ids, copy_numbers, solvent_content, nproc, 
                     no_timeout, force_af, skip_af, skip_autobuild, no_waters):
    global current_process
    
    if not fasta_file or not mtz_file:
        yield "Error: You must upload both a FASTA file and an MTZ file."
        return

    # Prepare arguments as a simple dictionary to pass across processes safely
    args_dict = {
        "fasta_path": fasta_file.name,
        "mtz_path": mtz_file.name,
        "output_path": os.path.abspath(output_dir.strip()) if output_dir else os.path.abspath("results_gui"),
        "uniprot_ids": uniprot_ids.strip() if uniprot_ids else None,
        "copy_numbers": copy_numbers.strip() if copy_numbers else None,
        "solvent_content": float(solvent_content) if solvent_content else None,
        "nproc": int(nproc),
        "no_timeout": no_timeout,
        "force_af_cluster": force_af,
        "skip_af_cluster": skip_af,
        "skip_autobuild": skip_autobuild,
        "no_waters": no_waters,
        "reference_model": None,
        "reference_map": None
    }

    os.makedirs(args_dict["output_path"], exist_ok=True)

    # Create a multiprocessing Queue and start the worker process
    log_queue = multiprocessing.Queue()
    current_process = multiprocessing.Process(target=pipeline_worker, args=(args_dict, log_queue))
    current_process.start()

    accumulated_logs = "🚀 Starting pipeline...\n\n"
    yield accumulated_logs

    # Stream UI updates every 0.5 seconds
    while current_process.is_alive() or not log_queue.empty():
        while not log_queue.empty():
            msg = log_queue.get()
            accumulated_logs += msg
        yield accumulated_logs
        time.sleep(0.5)

    # Safety check in case the process was killed externally
    if current_process.exitcode != 0 and current_process.exitcode != -15 and "✅ PIPELINE" not in accumulated_logs and "❌ PIPELINE FAILED" not in accumulated_logs:
        accumulated_logs += f"\n\n⚠️ Process exited unexpectedly with code {current_process.exitcode}"
        yield accumulated_logs


def abort_pipeline(current_logs):
    """Kills the background process and updates the log console."""
    global current_process
    if current_process and current_process.is_alive():
        current_process.terminate() # Force kill the underlying Python process
        current_process.join()
        return current_logs + "\n\n🛑 PIPELINE ABORTED BY USER."
    return current_logs + "\n\n⚠️ No pipeline is currently running."


# --- Build the Gradio Web Interface ---
readable_theme = gr.themes.Soft(
    font=[gr.themes.GoogleFont("Inter"), "ui-sans-serif", "system-ui", "sans-serif"],
    font_mono=[gr.themes.GoogleFont("JetBrains Mono"), "ui-monospace", "Consolas", "monospace"],
    text_size=gr.themes.sizes.text_lg
)

with gr.Blocks(theme=readable_theme, title="AF-Guided MR Pipeline") as ui:
    gr.Markdown("# AlphaFold-Guided Molecular Replacement")
    gr.Markdown("Upload your target sequences and diffraction data to begin the automated structure solution pipeline.")
    
    with gr.Row():
        with gr.Column():
            gr.Markdown("### 1. Input Files & Destinations")
            fasta_input = gr.File(label="Target Sequences (.fasta)", file_types=[".fasta", ".fa"])
            mtz_input = gr.File(label="Diffraction Data (.mtz)", file_types=[".mtz"])
            output_dir_input = gr.Textbox(
                label="Output Directory Path", 
                value=os.path.abspath("results_gui"), 
                info="Specify where the pipeline should save the final model and map."
            )
            
            gr.Markdown("### 2. Protein Information")
            uniprot_input = gr.Textbox(label="UniProt IDs", placeholder="e.g., P00519, P12345 (comma-separated)")
            copy_input = gr.Textbox(label="Copy Numbers", placeholder="e.g., 2:1 (colon-separated)")
            
        with gr.Column():
            gr.Markdown("### 3. Execution Settings")
            nproc_input = gr.Slider(minimum=1, maximum=64, value=4, step=1, label="Number of Processors (nproc)")
            solvent_input = gr.Number(label="Solvent Content (Optional)", info="Leave blank to auto-calculate")
            
            with gr.Accordion("Advanced Settings", open=False):
                no_timeout_toggle = gr.Checkbox(label="Disable Phaser Timeout (no_timeout)", value=False)
                force_af_toggle = gr.Checkbox(label="Force AF Cluster (force_af_cluster)", value=False)
                skip_af_toggle = gr.Checkbox(label="Skip AF Cluster (skip_af_cluster)", value=False)
                skip_autobuild_toggle = gr.Checkbox(label="Skip AutoBuild (skip_autobuild)", value=False)
                no_waters_toggle = gr.Checkbox(label="No Waters in AutoBuild (no_waters)", value=False)

    with gr.Row():
        run_btn = gr.Button("🚀 Run Pipeline", variant="primary")
        abort_btn = gr.Button("🛑 Abort Job", variant="stop")
    
    gr.Markdown("### Status / Output")
    output_console = gr.Textbox(label="Execution Logs", lines=20, max_lines=40, interactive=False)

    # Wire up the Run Button (save the event so we can cancel it!)
    run_event = run_btn.click(
        fn=execute_pipeline,
        inputs=[fasta_input, mtz_input, output_dir_input, uniprot_input, copy_input, solvent_input, nproc_input, 
                no_timeout_toggle, force_af_toggle, skip_af_toggle, skip_autobuild_toggle, no_waters_toggle],
        outputs=output_console
    )

    # Wire up the Abort Button
    abort_btn.click(
        fn=abort_pipeline,
        inputs=[output_console],
        outputs=[output_console],
        cancels=[run_event] # This stops the UI from waiting for the generator!
    )

if __name__ == "__main__":
    ui.launch(server_name="0.0.0.0", server_port=7860, share=False)
