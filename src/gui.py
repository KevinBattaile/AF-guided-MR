import gradio as gr
import os
import argparse
import traceback
from af_guided_mr.pipeline import run_pipeline

def execute_pipeline(fasta_file, mtz_file, output_dir, uniprot_ids, copy_numbers, solvent_content, nproc, 
                     no_timeout, force_af, skip_af, skip_autobuild, no_waters):
    
    # 1. Validate mandatory files
    if not fasta_file or not mtz_file:
        return "Error: You must upload both a FASTA file and an MTZ file."

    # 2. Construct the mock 'args' object that pipeline.py expects
    class PipelineArgs:
        pass
    
    args = PipelineArgs()
    args.fasta_path = fasta_file.name
    args.mtz_path = mtz_file.name
    
    # Set the output path based on user input, default to 'results_gui' if left blank
    args.output_path = os.path.abspath(output_dir.strip()) if output_dir else os.path.abspath("results_gui")
    
    args.uniprot_ids = uniprot_ids.strip() if uniprot_ids else None
    args.copy_numbers = copy_numbers.strip() if copy_numbers else None
    args.solvent_content = float(solvent_content) if solvent_content else None
    args.nproc = int(nproc)
    args.no_timeout = no_timeout
    args.force_af_cluster = force_af
    args.skip_af_cluster = skip_af
    args.skip_autobuild = skip_autobuild
    args.no_waters = no_waters
    args.reference_model = None
    args.reference_map = None

    # Ensure the output directory exists
    os.makedirs(args.output_path, exist_ok=True)

    # 3. Run the pipeline and catch any errors
    try:
        run_pipeline(args)
        return f"Pipeline completed successfully! \nResults saved to: {args.output_path}"
    except Exception as e:
        error_trace = traceback.format_exc()
        return f"Pipeline failed with error:\n{str(e)}\n\nTraceback:\n{error_trace}"

# --- Build the Gradio Web Interface ---
with gr.Blocks(theme=gr.themes.Soft(), title="AF-Guided MR Pipeline") as ui:
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

    run_btn = gr.Button("🚀 Run Pipeline", variant="primary")
    
    gr.Markdown("### Status / Output")
    output_console = gr.Textbox(label="Execution Logs", lines=6, interactive=False)

    # Wire the button to the function
    run_btn.click(
        fn=execute_pipeline,
        inputs=[fasta_input, mtz_input, output_dir_input, uniprot_input, copy_input, solvent_input, nproc_input, 
                no_timeout_toggle, force_af_toggle, skip_af_toggle, skip_autobuild_toggle, no_waters_toggle],
        outputs=output_console
    )

if __name__ == "__main__":
    ui.launch(server_name="0.0.0.0", server_port=7860, share=False)
