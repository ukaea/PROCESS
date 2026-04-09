
# GitHub Actions

Our GitHub actions Continuous Integration (CI) pipeline serves to ensure each branch and pull request conforms to our testing and style requirements. All of our CI can currently be found in one workflow file: `process.yml`. A summary of each job within this workflow can be found below:

| Name | Functionality |
| ---- | ------------- |
| unit-test | Installs PROCESS and runs the unit tests. The job will fail if any of the unit tests fail. |
| integration-test | Installs PROCESS and runs the integration tests. The job will fail if any of the integration tests fail. |
| regression-test | Installs PROCESS and runs the regression tests with a 0.2% and 5% tolerance, respectively. The job will fail if any of the regression tests fail. The job uses tracked MFILEs for regression test comparisons. |
| run-tracking-inputs | Installs PROCESS and runs the regression test input files, archiving the output MFILEs. Only runs on the **main** branch. |
| tracking | Collects MFILEs for input files of interest and creates a dashboard of changes in key values over time (one datapoint for each commit on main). Only runs on the **main** branch. |
| pre-commit-quality-check | ensures the pushed code meets our standards as defined in `.pre-commit-config.yaml`. |
| docs | Builds and deploys the documentation onto GitHub pages. |

!!! Info "Regression Job"
    When creating a PR, the tracked MFILEs from the point at which you branched off main will be used during the regression test job. This job may fail on your PR. If so, these changes to the output will need to be reviewed. If the changes are accepted (i.e. they come from a necessary change to a model), upon merge to main the 
    `run-tracking-inputs` job will run the regression input files to generate the associated MFILEs to be tracked by the 
    `tracking job`, and these will now reflect the changes your PR has introduced. Therefore when the regression job runs upon merge of your code into main, it will use these updated tracked MFILEs, and the regression tests will pass.