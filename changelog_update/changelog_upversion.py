import argparse
import sys
import yaml
import jinja2
from pathlib import Path
from datetime import date


CURRENTDIR = Path(__file__).resolve().parent
# List of all .yaml files in the current directory.

if __name__ == "__main__":

    yaml_files = CURRENTDIR.glob("*.yaml")

    parser = argparse.ArgumentParser(
        description="Create the Changelog markdown for a new Process release"
    )
    parser.add_argument("version", type=str, help="New version of process")
    args = parser.parse_args()

    # Creating the version argument to be used in the command line so the user can define which version of Process is being implemented

    changelog_dict = {
        "Added": [],
        "Changed": [],
        "Deprecated": [],
        "Fixed": [],
        "Removed": [],
    }

    # Defining the changloge dictionary that is the same structure as the yaml files that users have generated using the
    # create_changelog.py script. Each header is able to have mutliple changes.

    for file in yaml_files:
        with open(file, "r") as stream:
            try:
                yaml_file = yaml.safe_load(stream)

                for header in changelog_dict:
                    changelog_dict[header] += yaml_file[header]

            except yaml.YAMLError:
                print(
                    f"Failed to open and include content of {file} - please check headers match changelog_dict and that the structure of your yaml file matches that in `README.md` (i.e. space between hyphen and entry)."
                )
                sys.exit()

    # Looping through and reading YAML files in CURRENTDIR
    # Here, headers are the different values under which an entry to the changelog can be made i.e 'Fixed', 'Removed' etc. The script will exit
    # if the structure of one of the .yaml files is incorrect.

    for header in changelog_dict:
        changelog_dict[header] = [
            entry for entry in changelog_dict[header] if entry is not None
        ]

    changelog_dict = {
        header: entries for header, entries in changelog_dict.items() if entries
    }

    # Only keeping the headers in the changelog dictionary if they have a changelog entry associated with them.
    # So if the 'Added' header has no entries between the .yaml files, it is not included in the changelog dict for use in the changelog.

    environment = jinja2.Environment(
        loader=jinja2.FileSystemLoader(str(Path(__file__).resolve().parent))
    )

    template = environment.get_template("changelog_template.jinja2")

    context = {
        "date_of_upversion": date.today(),
        "process_version": args.version,
        "changelog_dict": changelog_dict,
    }

    # The changelog_template.jinja2 file holds the structure for producing the markdown file that is the output of this script.
    # It creates a .md file with the headers and values from users' yaml scripts and is then available for copy and pasting into the
    # changelog upon process upversioning. It takes todays date so you know when the script/changelog was updated and asks the user
    # for the version of process that is being utilised in the upgrade via the command line ArgParse method.

    with open(CURRENTDIR / "changelog_update.md", "w") as f:
        f.write(template.render(**context))

    # Using the above logic, a file named 'changelog_update.md' is output in the CURRENTDIR- a markdown file.

    files_to_remove = CURRENTDIR.glob("*.yaml")
    for f in files_to_remove:
        f.unlink(missing_ok=True)

    # Checks the files in the CURRENTDIR after creating the changelog_update.md file for this upversion. It then finds which
    # are .yaml files and removes them so that once they are uploaded, they are not stored and get uploaded twice in the next merge.
    # This is safe as if the .yaml files did need to be accessed they would be in a previous commit.
