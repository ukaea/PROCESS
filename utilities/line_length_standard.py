"""

    Script to check line length of repository files

    Files include:
        - *.f90
        - *.py
        - *.tex

    James Morris
    UKAEA
    26.09.19

"""

import os

# Get repository root directory
import pathlib
import time

timeout = time.time() + 10  # 10 seconds
found_root = False
back = ""
while not found_root:
    if time.time() > timeout:
        print(
            "Can't find repository root. Make sure utility is being run "
            "inside a PROCESS repository clone"
        )
        break
    else:
        my_file = pathlib.Path(back + ".gitignore")
        if my_file.is_file():
            found_root = True
            if back == "":
                REPO_ROOT = "."
            else:
                REPO_ROOT = back
        back += "../"

# Set line limit
LINE_LIMIT = 100

# Paths and files to check
CHECK_LIST = {
    "fortran": {"path": REPO_ROOT + "/source/fortran/", "extension": ".f90"},
    "python": {"path": REPO_ROOT + "/utilities/", "extension": ".py"},
    "latex": {"path": REPO_ROOT + "/documentation/", "extension": ".tex"},
}

if __name__ == "__main__":

    # Intro message

    intro = """

    Checking line length standard (<100) lines.

    For the following locations

        - /source/fortran/*.f90
        - /source/utilities/*.py
        - /source/documentation/*.tex

    J. Morris
    UKAEA
    26.07.19

    """

    print(intro)

    # Set status
    STATUS = True
    LINE_COUNT = 0
    FILE_COUNT = 0

    # Check each entry in CHECK_LIST
    for key in CHECK_LIST.keys():
        print("")
        print("# Checking {0} files".format(key))
        print("")
        check_path = CHECK_LIST[key]["path"]
        for filename in os.listdir(check_path):
            if CHECK_LIST[key]["extension"] in filename:
                FORT_STATUS = True
                print("")
                print(filename)
                file_lines = open(
                    check_path + filename, "r", encoding="UTF-8"
                ).readlines()
                counter = 0
                for line in file_lines:
                    counter += 1
                    if len(line) > LINE_LIMIT:
                        print(
                            "|-- {0:30} :: line :: {1:<5} :: length = {2:<10}".format(
                                filename, counter, len(line)
                            )
                        )
                        STATUS = False
                        LINE_COUNT += 1
                        if FORT_STATUS:
                            FILE_COUNT += 1
                        FORT_STATUS = False

    # if STATUS:
    #     sys.exit(0)
    # else:
    #     sys.exit("\nERROR: Line length exceeded {0} times in {1} files".
    #              format(LINE_COUNT, FILE_COUNT))
