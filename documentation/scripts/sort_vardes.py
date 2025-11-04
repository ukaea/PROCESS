"""

Script to tidy up vardes.md for the GitLab Page

J. Morris
10/08/19
UKAEA

"""

import sys

with open(sys.argv[1]) as vardes:
    lines = vardes.readlines()

new_lines = []

for counter, line in enumerate(lines):
    if "PROCESS Variable Descriptor File : dated" in line:
        date = line.split("dated")[-1].replace(" ", "")
        new_date = date[:4] + "." + date[4:6] + "." + date[6:]
        new_line = f"# PROCESS Variable Descriptions {new_date}\n"
    elif counter == 3:
        new_line = "## Introduction"
    elif "###" in line:
        new_line = line.split("]")[0].replace("[", "").replace(r"\_", " ")
    else:
        new_line = line

    new_lines.append(new_line)

with open(sys.argv[1], "w") as new_vardes:
    new_vardes.writelines(new_lines)
    new_vardes.close()
