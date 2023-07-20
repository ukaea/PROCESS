# Coverage

Code coverage is tracked using `gcov` and the data files are generated when PROCESS is run.
The coverage report is generated using `lcov` and is stored as webpages in the `html` folder inside
the `lcov_results` folder, to access it simply open `index.html`. The coverage data inside the `.gcda` files
is accumulative, that is to say that if PROCESS is run multiple times the `.gcda` files will have
coverage data for all of the PROCESS runs combined. There is no way to separate these, so if coverage data for
individual runs needs to be recovered then the best method is to clear existing `.gcda` files, run PROCESS,
generate the coverage report and copy the contents of the `lcov_results` folder out of the `process` directory
and store it separately.

<h2>Generate coverage report</h2>

To generate a coverage report, use the command:

```bash
cmake --build build --target coverage
```

This must be done after running PROCESS at least once.

<h2>Clear existing coverage data</h2>

To clear existing coverage data (`.gcda` files) before a new run of PROCESS, use:

```bash
cmake --build build --target coverage_cleanup
```

Note that all existing coverage data is cleared on every re-build of PROCESS.
