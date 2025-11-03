'''
Start PROCESS run - select input files by prefix'''

from process.main import SingleRun

import subprocess
from pdf2image import convert_from_path
from process.io import plot_proc
import os

script_dir = os.path.dirname(os.path.realpath(__file__))

prefix = "/squid"



def postprocess(single_run):
    # Postprocess the results
    #print(single_run.mfile_path)

    """ plot_proc uses command line arguments of the current process. 
    Jupyter adds command line arguments under the hood causing plot_proc to fail. 
    Running plot proc in its own process isolates it from the jupyter command line arguments """
    subprocess.run(["python", plot_proc.__file__, "-f", str(single_run.mfile_path)])

    # Create a summary PDF
    # Convert PDF to PNG in order to display in notebook
    summary_pdf = str(single_run.mfile_path) + "SUMMARY.pdf"
    print(summary_pdf)
    pages = convert_from_path(summary_pdf)
    for page_no, page_image in enumerate(pages):
        png_path = script_dir / f"plot_proc_{page_no + 1}.png"
        page_image.save(png_path, "PNG")


if __name__ == "__main__":
    # Run process on an input file
    single_run = SingleRun(script_dir+prefix+".IN.DAT")
    single_run.run()

    # vary_run = VaryRun(script_dir+prefix+".IN.DAT")
    # vary_run.run()

    # Generate pdf with results
    # postprocess(single_run)