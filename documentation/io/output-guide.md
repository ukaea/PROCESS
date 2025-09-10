
# Output Files

The default names of the output files are given below.  The filenames can be modified as 
described in [File Naming Convention](input-guide.md#file-naming-convention)

The main output from the code is sent to a text file `OUT.DAT` in the working directory. 
It is essential to check that the code reports that it has found a *feasible solution*.  

A second text file, `MFILE.DAT`, is also produced in the working directory. This file contains 
most of the same data as `OUT.DAT` but in a different format and has been designed to be 
"machine-readable" by some of the [utility programs](utilities.md) for 
post-processing and graphical output.

An additional diagnostic output text file `VFILE.DAT` is generated if the `verbose` flag is set 
to 1 in the input file.   

(If there is already a file with any of the above filenames it will be overwritten.)

Utilities generate additional output files.
