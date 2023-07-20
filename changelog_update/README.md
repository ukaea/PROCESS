# Instructions for the RSE updating the CHANGELOG on Process Upversioning

### Step 1
*This process should only be followed upon PROCESS upgrade i.e. from `2.4.0` to `2.5.0`.*

The Python script uses `ArgParse` with the argument `version` to specify the version being unversioned to. This requires an extra argument to be included in the command when running the script. It takes a string as the `version` argument. For example, if you were up versioning to Process 2.5.1, you would use the following to run the Python script:

From the terminal, run:
```bash
python changelog_update/changelog_upversion.py 2.5.1
```
This script will take the .yaml files in the `changelog_update` directory and using the jinja2 template, will combine the changes since the last Process upgrade by header as defined in the changelog dict (the headers seen above).

A markdown file- `changelog_update.md`, will be created in the same directory and all the YAML files will be removed so they cannot be uploaded more than once.

### Step 2
Navigate to the `changelog_update.md` file which should include the new version number (as input by you as an argumument in the command line), the date of the running of the script (should be today), and then the content of the changes since the last upgrade.

Now, the content of this file should be copied and pasted into the CHANGELOG.md file that already exists in the process root directory. 