"""Convert an existing IN.DAT file to YAML

This script uses the in_dat module, which parses an existing IN.DAT file and
creates an input file data object, containing an input data dict. This dict is
useful within in_dat for manipulating the data and writing new IN.DAT files.

create_formatted_input_data_dict() then processes and restructures that dict to
create formatted_input_data_dict, which contains only the information extracted
from IN.DAT. It is formatted in a similar manner to the original input file, but
with some section heirarchy applied. This dict is then dumped to file as YAML
(input.yml).

The motivation for this is that this input YAML file is easy to load as a dict,
which can readily be used in Python input file checks. Perhaps this more
structured format will begin to be used over the IN.DAT format, which is more
complex to parse and check.
"""

import yaml
from process_io_lib import in_dat

# Input file to read (.DAT) and output to write (.yml)
# If filename is provided, look in current working directory
# Absolute and relative file paths also work
INPUT_FILE_PATH = "IN.DAT"
OUTPUT_FILE_PATH = "input.yml"

# Create InDat object by reading in .DAT input file
input_data = in_dat.StructuredInputData(filename=INPUT_FILE_PATH)

# Dump YAML into output file (which is a YAML input file)
output_file = open(OUTPUT_FILE_PATH, "w")

try:
    yaml.dump(input_data, output_file)
    print("Successfully converted input file to YAML, saved as " f"{OUTPUT_FILE_PATH}")
except Exception:
    print("Failed to write YAML input file")

output_file.close()
