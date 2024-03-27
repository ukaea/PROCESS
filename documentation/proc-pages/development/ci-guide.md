
# GitHub Actions

Our GitHub actions Continuous Integration (CI) pipeline serves to ensure each branch and pull request conforms to our testing and style requirements. Due to the requirement of all stages on the built PROCESS artefacts and Docker image, all of our CI can currently be found in one workflow file: `process.yml`. A summary of each job within this workflow can be found below:

| Name | Functionality |
| ---- | ------------- |
| docker | Checks if the `process-ci` Docker container is up-to-date and builds it if not. Only runs on the **main** branch. |
| make-py38 | Builds and archives the PROCESS build artefacts |
| unit-py38 | Installs PROCESS and runs the unit tests. The job will fail if any of the unit tests fail. |
| integration-py38 | Installs PROCESS and runs the integration tests. The job will fail if any of the integration tests fail. |
| regression-py38 | Installs PROCESS and runs the regression tests with a 0% and 5% tolerance, respectively. The job will fail if any of the regression tests fail. |
| large-tokamak-py38 | Installs PROCESS and runs the `large-tokamak` input file, archiving the output MFILE. Only runs on the **main** branch. |
| flake8 | Runs the flake8 Python linter and fails if any lint errors occur. |
| black | Runs the black Python formatter and fails if any formatting issues are detected. |
| tracking | Collects MFILEs for input files of interest and creates a dashboard of changes in key values over time (one datapoint for each commit on main). |
| docs | Builds and deploys the documentation onto GitHub pages. |