# Pre-commit

!!! Info "TL;DR: the quality job fails"
    If the quality job fails, run pre-commit using `pre-commit run --all-files` and fix any failures (some issues may be fixed automatically). Add and commit the changes until pre-commit no longer fails. If this doesn't solve
    the problem then see below for more details.


Pre-commit is PROCESS' way of ensuring all code pushed to our repository is of a certain quality 
and style. One common style ensures PROCESS can be easily navigated and understood by all users. It 
also reduces the "diffs" when reviewing code changes as often small style changes (such as the 
addition of a new-line) will litter the "diffs" and make reviewing harder.

[pre-commit website](https://pre-commit.com/#top_level-files){ .md-button }

Pre-commit works on staged files (ie those that have been `git add`ed) after you issue the 
`git commit` command but before the commit is actually made. If any pre-commit job **failed** then 
the commit will not be made. On a failure, one of two things can happen:

1. Pre-commit plugins will rectify the mistakes made. This will happen with code formatters 
   (whose job it is to edit your files to the correct style). The files the plugins have changed 
    can then be `git add`ed again and the `git commit` command re-issued.
2. A pre-commit plugin will identify the mistake but will NOT fix it. This could happen with 
   `ruff` which is a linter and warns of code mistakes but does not correct them. You will need 
   to fix these mistakes manually. Now, the changed files can be `git add`ed and the `git commit` command re-issued.

!!! Info "VSCode GUI users"
    On the VSCode GUI, you will see that files pre-commit modifies will change from green as 
    there are changes that need to made. Follow your usual workflow to re-add these files and re-do the commit as usual.

!!! Tip
    I would advise you become familiar with the *What does pre-commit check for?* section of this 
    document. This will allow you to understand whether a mistake has been automatically fixed or not.


## Installation

Check if pre-commit is installed in your environment by running:

```bash
pre-commit -h
```

If pre-commit is not installed, then it can be installed by running:

```bash
pip install -r requirements_dev.txt
```

in your environment, which will install all development dependencies including pre-commit. 
Alternatively, pre-commit can be installed on its own by running 

```bash
pip install pre-commit
```

You should then run
```bash
pre-commit install
```
which will install pre-commit as a **Git pre-commit hook**. When pre-commit is a hook, it will run on all of the 
files you have changed before allowing a commit -- if any problems are identified, they will be 
fixed and you will need to re-add the files that pre-commit has changed.

!!! example "Adding two files"
    Consider that two files are being `git add`ed.
    One of the files, `foo.py` has stylistic changes which **ruff** objects to.

    ```
    > git add foo.py bar.py
    > git commit -m "Adding foo and bar"
        Trim Trailing Whitespace.................................................Passed
        Check for merge conflicts................................................Passed
        Debug Statements (Python)................................................Passed
        ruff.....................................................................Failed
            - hook id: ruff-format
            - exit code: 1
            - files were modified by this hook

            Fixing foo.py
        Format YAML files....................................(no files to check)Skipped

    > git add foo.py # since ruff has modified foo.py
    > git commit -m "Adding foo and bar"
    ```

## Pre-commit and the `quality` CI stage
The Process continuous integration system (used to check Process builds and passes tests) also has 
a `quality` stage. This is where several jobs will run to ensure the quality of your code. If all 
your commits pass through pre-commit, then these jobs should not fail as your code will be of a high quality.

## Using ruff with VSCode
Although not required, the `ruff` VSCode extension will ensure that all the Python files you save 
will be ruff-compliant and, as such, won't need to modified by pre-commit.

Open or create the file `.vscode/settings.json` and add/modify the following settings:
```json
{
    "[python]": {
        "editor.defaultFormatter": "charliermarsh.ruff",
        "editor.formatOnSave": true
    }
}
```

## What does pre-commit check for?

### General style on all files

Pre-commit performs a few checks on each and every file you add, regardless of type.

* `end-of-file-fixer` will check that each file ends with exactly one new line. **This plugin will automatically fix any mistakes it finds**.
* `trailing-whitespace` will check that there is no excess whitespace (spaces or tabs) at the 
  end of any line. **This plugin will automatically fix any mistakes it finds**.
* `check-merge-conflict` checks that all merge conflict's have been resolved. **This plugin will NOT automatically fix any mistakes it finds**.

### Python checks

* [`ruff`](https://github.com/astral-sh/ruff) formats and lints code. It will identify formatting and stylistic errors in Python code. **This plugin will automatically fix any mistakes it finds**.

### Other checks

* [`yamlfmt`](https://github.com/jumanjihouse/pre-commit-hook-yamlfmt) formats YAML code (similar 
  to what `ruff` does for Python code). **This plugin will automatically fix any mistakes it finds**.
