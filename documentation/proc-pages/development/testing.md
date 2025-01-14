# Testing

## Understanding testing

The purpose of testing is to check that the existing code behaves as it is expected to, and that any 
code changes don't produce unexpected results, such as breaking pre-existing functionality or 
creating an unwanted side-effect. It allows changes to be made with confidence and it increases 
confidence in the results produced by the code. 

Tests run part of the code with a given input and assert that the observed output is equal to an 
expected output. If the assertion is true the test will pass, if it is false, it will fail. A 
collection of tests in known as a test suite, and they can be classified into three different testing scopes:

### Unit tests

Unit tests test the smallest possible piece of code, like a function or method. It might call a 
function with some arguments and assert that it returns what is expected.

### Integration tests

Integration tests test larger pieces of code that make up parts of the program, such as classes, an 
input file reader or a plotting utility. It might create an instance of a class with an input 
filename, then call a method to create a plot file, then check that the plot was created without throwing any exceptions.

### Regression tests

Regression tests detect changes in the entire program's behaviour by checking that for a given 
input, it produces exactly the same output as before those changes. It detects changes in the 
program's output. Therefore if your code changes the program's output, it will fail the regression 
test. In this case that output difference will need to be reviewed and if accepted, the expected (or 
reference) will be updated.

Regression tests compare the output of PROCESS locally when running an input file to a reference output stored in a 
[repository](https://github.com/timothy-nunn/process-tracking-data). The test suite will download the reference output for the commit where the your current branch branched off of main. This means **each branch is accountable for only its changes since it branched off of main**. 

## pytest

Process uses the `pytest` testing framework in its test suite. `pytest` tests are modular, quick to 
write with little code and produce helpful information when they fail. It is used widely in the Python world.

### Running pytest

`pytest` can be run locally by running `pytest` in the project root directory. This will run all 
tests. `pytest` can also be configured to run in the sidebar of VS Code. 

Individual test collections can be run by specifying the test directory to run, e.g.

```BASH
pytest tests/unit 
```
will only run unit tests.

Furthermore, the `-k` can be used to match tests within a test collection, e.g.

```BASH
pytest tests/regression/ -k "large_tokamak"
```

The Continuous Integration (CI) system also runs the `pytest` test suite in the `testing` stage of 
the pipeline. Unit, integration and regression tests are run as separate jobs to make it easier to see where failures lie.

### Custom pytest options

`pytest` can also be run with various options custom to PROCESS.

```bash
pytest --solver=<"legacy-vmcon"|"new-vmcon">
```

Allows different solvers to be used.

```bash
pytest --opt-params-only
```

Only compares final optimisation parameters (solution vector) and the objective function value 
(figure of merit) in regression tests. `pytest -k baseline_jan --opt-params-only` will do so only 
for the "baseline_jan_2017" test. This is particularly useful when comparing solutions obtained 
from different solvers, for example.

### Test failures

Unit and integration tests must pass. Any that fail require source or test modification. Regression 
tests may fail, however; their purpose is to make you aware of significant changes to results as a 
result of your source changes. A 5% tolerance is applied by default to the regression tests: if any 
values differ by >5% from the reference, the regression test will fail for that scenario. This new 
value may be the desired result of the changes, however. Optionally, a 0% tolerance regression 
test can be run using `pytest --reg-tolerance=0`.

By default, pytest will only display warnings and above that are logged in PROCESS (not info or debug logs). However, lower level logging may be appropriate when tests are failing. To capture lower level logs, run pytest with either of the following options `--log-level=INFO` or `--log-level=DEBUG`.

It is incumbent on the author to check the test results created by their code changes, and modify 
source or tests if required. Are the regression changes expected and acceptable?

For a guide on contributing code to PROCESS, see `CONTRIBUTING.md`.

### Speeding up tests

Running the entire test suite can be time consuming, as by default it runs on a single core. 
`pytest-xdist` allows `pytest` tests to be distributed across multiple cores to speed up testing.

`pytest-xdist` should be installed already (included in `requirements.txt`), but if not it can be installed manually with:

```bash
pip install pytest-xdist
```

To run tests on as many processes as your computer has CPU cores, use:

```bash
pytest -n auto
```

This can result in considerable speedups when running tests. Normal `pytest` commands can also be combined, e.g.

```bash
pytest -n auto -k regression
```

runs just the regression tests on all available cores.

## Test coverage

Test coverage (in Python only) is provided in a badge on the repository homepage. A report can 
also be generated locally. A development (editable) pip installation (which is run by default by the `cmake` build script) ensures that `pytest` and `pytest-cov` will use the same installed location of Process. Then:

```bash
pytest --cov=process tests/unit/
```

runs the unit tests and produces a coverage report for them.

## pytest failures on older OS's

As discussed in the Installation guide, PROCESS is dependant on a number of dynamically linked 
libraries. The versions of these libraries are different on different versions of OS's. This 
introduces floating-point differences in the code which can propogate and show tests failing by 
~0.70%. The cause of such issues has been isolated and will be highlighted by a warning message 
when running pytest:

```
You are running the PROCESS test suite on an outdated system.
This can cause floating point rounding errors in regression tests.

Please see documentation for information on running PROCESS (and tests)
using a Docker/Singularity container.
```

It is suggested that PROCESS is run, built, and tested via a container when not using Ubuntu 20.

## Reasoning behind the CONTRIBUTING.md method

The method in the CONTRIBUTING.md is standard apart from how regression tests are handled; this is 
explained below.

When changing or reviewing code, it is important to understand the effect those changes will have 
on the results of various reference input files, called regression scenarios. A regression test 
compares the current observed results with previous ones called references; in this case the references 
are the results when the branch was first created. Comparing the observed and reference results is 
important when assessing the magnitude of intended changes to the output in certain scenarios, as 
well as detecting large unintended changes or failures to solve. The impact of the code changes on 
the regression results needs to be visible to the author and reviewer of the merge request.

Typically the solver will arrive at a very slightly different solution when any change is made, 
and so there are typically large numbers of very small changes with a few more significant changes 
amongst them. It is therefore useful for the reviewer to be able to filter out those more significant 
changes, say >5%, in order to understand the more significant effects of the code changes. This is 
why a 5% tolerance regression test job is used in the CI system; it will fail and report if any 
values in any scenario differ from the reference values by >5%.

### Running tests locally

When running `pytest` locally, by default the 5% tolerance regression tests are run. These will 
fail if the code changes on that branch cause a regression scenario result value to change by >5%. 
If this is intended, this is fine; a regression test failure is not wrong, it informs you that 
something has changed, in this case by >5%.

### Running the CI on a branch

When those local changes are committed and pushed, the CI system for the branch runs. This runs 5% 
and 0% tolerance regression jobs, which are allowed to fail. This shows the author and reviewer 
what the changes to the regression results are as a result of the code changes on that branch.

## Drawbacks to this approach

!!! note
  

- In time, it may be better to use a data repository for the regression references.
- It's possible that two regression-acceptable branches can be merged to make a regression-unacceptable 
  develop branch, but as regression tests aren't run on develop, only overwritten, this wouldn't produce 
  a failure. Perhaps checking for regression failures to solve or significant changes in certain key 
  variables on develop would help catch these cases. Monitoring of the 
  [tracker](http://process.gitpages.ccfe.ac.uk/process/tracking.html) will help detect these in the 
  meantime. Otherwise it is felt that this vulnerability could only be addressed by being 
  unnecessarily cautious when merging, and is outweighed by the ease of using a regular git flow as 
  outlined in the `CONTRIBUTING.md`.
