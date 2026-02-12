
# GitHub Actions

Our GitHub actions Continuous Integration (CI) pipeline serves to ensure each branch and pull request conforms to our testing and style requirements. Due to the requirement of all stages on the built PROCESS artefacts and Docker image, all of our CI can currently be found in one workflow file: `process.yml`. A summary of each job within this workflow can be found below:

| Name | Functionality |
| ---- | ------------- |
| unit-test | Installs PROCESS and runs the unit tests. The job will fail if any of the unit tests fail. |
| integration-test | Installs PROCESS and runs the integration tests. The job will fail if any of the integration tests fail. |
| regression-test | Installs PROCESS and runs the regression tests with a 0% and 5% tolerance, respectively. The job will fail if any of the regression tests fail. |
| run-tracking-inputs | Installs PROCESS and runs the regression test input files, archiving the output MFILEs. Only runs on the **main** branch. |
| tracking | Collects MFILEs for input files of interest and creates a dashboard of changes in key values over time (one datapoint for each commit on main). Only runs on the **main** branch. |
| pre-commit-quality-check | ensures the pushed code meets our standards as defined in `.pre-commit-config.yaml`. |
| docs | Builds and deploys the documentation onto GitHub pages. |