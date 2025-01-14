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
counter = 0

for line in lines:
    if "PROCESS Variable Descriptor File : dated" in line:
        date = line.split("dated")[-1].replace(" ", "")
        new_date = date[:4] + "." + date[4:6] + "." + date[6:]
        new_line = "# PROCESS Variable Descriptions {0}\n".format(new_date)
    elif counter == 3:
        new_line = "## Introduction"
    elif "###" in line:
        new_line = line.split("]")[0].replace("[", "").replace(r"\_", " ")
    else:
        new_line = line

    new_lines.append(new_line)
    counter += 1

new_vardes = open(sys.argv[1], "w")
for item in new_lines:
    new_vardes.write(item)
new_vardes.close()
