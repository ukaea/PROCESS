import sys
import jinja2
from pathlib import Path
from subprocess import run

CURRENTDIR = Path(__file__).resolve().parent
# Defining the currentdir such that the script is always executed within the correct 'ChangeLog_Update' directory and
# that the YAML files are output here.


if __name__ == "__main__":

    branch_name = run(
        "git rev-parse --abbrev-ref HEAD",
        shell=True,
        capture_output=True,
        cwd=CURRENTDIR,
        check=True,
    )

    git_branch_name = branch_name.stdout.decode()
    git_branch_no = branch_name.stdout.decode().split("-")[0]

    if git_branch_no.isdigit() is True:
        git_branch_no == git_branch_no
    else:
        git_branch_no == git_branch_name

    # Finding the name of the current git branch which shulde have the follwowing structure
    # 1234-issue-name-from-merge-request
    # Then the name is split by using the '-' as the partition.
    # The number is wanted to use as the name for the .yaml file so this is taken as the 'git_branch_no'

    generated_yaml_file = Path(CURRENTDIR / f"{git_branch_no}.yaml")

    # This is just defining the yaml file that will be created by the running of this script.
    # Showing that the generated_yaml_file will always be in the currentdir and named as git_branch_no.yaml e.g 1234.yaml

    if generated_yaml_file.is_file():
        print(
            "This YAML file already exists - please locate it in the current directory and edit there."
        )
        sys.exit()
    else:
        print(
            f"Your file will be created in the changeLog_update directory and called: {git_branch_no}.yaml"
        )

    # The if statement checks if this script has already been run in the current branch and therefore if the yaml file you are about
    # to create already exists. If this is the case, then the yaml file should be located in currentdir and edited with any changes
    # there. This prevents overwriting of any changes that have already been documented in the yaml file.
    # If the yaml file hasnt been created, the console will print the name your yaml file will have so you are able to locate it in
    # the current directory.

    environment = jinja2.Environment(
        loader=jinja2.FileSystemLoader(str(Path(__file__).resolve().parent))
    )

    yaml_template = environment.get_template("yaml_template.jinja2")

    # Using jinja2 to template the standard .yaml file that is to be generated (empty) for user to fill out with their changes.

    with open(CURRENTDIR / f"{git_branch_no}.yaml", "w") as f:
        f.write(yaml_template.render())

    # Using the jinja2 template to write out a .yaml file in the current directory with the name being the users
    # git branch number.
