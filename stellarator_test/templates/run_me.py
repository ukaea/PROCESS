from process.main import SingleRun, VaryRun
from process.io import plot_proc

from pdf2image import convert_from_path
from pathlib import Path
import argparse
import subprocess
import os, sys

if __name__ == "__main__":

    script_dir = os.path.dirname(os.path.realpath(__file__))

    parser = argparse.ArgumentParser(
        prog='run_PROCESS',
        description="Run PROCESS with IN.DAT file present in the same directory",
    )
    parser.add_argument("-n", "--input_name", "-v")
    args = parser.parse_args()

    if args.input_name is not None:
        prefix = args.input_name
    else:
        prefix = "squid"
        # prefix = "transition"
        # prefix = "updated"
        # prefix = "rebuild"
        # prefix = "helias5"
        

    # Run process on an input file

    single_run = SingleRun(script_dir+'/'+prefix+".IN.DAT")
    single_run.run()

    # vary_run = VaryRun(script_dir+'/'+prefix+".IN.DAT")
    # vary_run.run()

    # Generate pdf with results
    # postprocess(single_run)

def postprocess(single_run):
    # Postprocess the results
    #print(single_run.mfile_path)

    # plot_proc uses command line arguments of the current process. Running plot proc in its own process isolates it from the command line arguments
    subprocess.run(["python", plot_proc.__file__, "-f", str(single_run.mfile_path)])

    # Create a summary PDF
    # Convert PDF to PNG in order to display in notebook
    summary_pdf = str(single_run.mfile_path) + "SUMMARY.pdf"
    print(summary_pdf)
    pages = convert_from_path(summary_pdf)
    for page_no, page_image in enumerate(pages):
        png_path = script_dir / f"plot_proc_{page_no + 1}.png"
        page_image.save(png_path, "PNG")
