# Contributing to `PROCESS`
There are many valuable contributions that can be made to PROCESS:
* Reporting bugs.
* Requesting/ recommending features.
* Implementing features or bugs raised as issues.
* Updating and improving documentation.

When contributing to this repository, please first discuss the change you wish to make via [issue](https://github.com/ukaea/PROCESS/issues), [discussion](https://github.com/ukaea/PROCESS/discussions), [email](https://github.com/ukaea/PROCESS#contacts), or any other method with the owners of this repository before making a change.

Please remember that all contributions and communication regarding PROCESS are subject to our [Code of Conduct](https://github.com/ukaea/PROCESS/blob/main/CODE_OF_CONDUCT.md).

## Creating an issue
Issues can be used to report bugs or request features and improvements. We ask you help us manage our issues by:
* Verifying your issue is not a duplicate of another issue; in this case, we welcome your contribution as a reply to the issue.
* Ensure you completely describe the bug/feature/problem and complete the given templates with appropriate detail.

## Submitting a pull request
Please discuss any feature ideas you have with the developers before submitting them, as you may not be aware of parallel development work taking place, or implementation decisions / subtleties which are relevant. The ideal workflow for submitting a pull request is as follows:

* Discuss the feature with a core PROCESS developer (as mentioned above).
* Submit an issue (if one does not exist for this feature/ bug) that documents the proposed change.
* Fork our repository.
* Create a branch off `main` with an appropriate name (e.g. `feature-abc`).
* Make the relevant changes for the repository (ensuring the changes do not creep away from the scope of the issue). Also ensure that the implementation follows the PROCESS [style guide](https://ukaea.github.io/PROCESS/development/standards/) for names of variables and functions etc.
* Discuss any problems or development choices in the issue and keep everyone updated on the progress of your development.
* If the changes are notable and it would benefit other users to be aware, [create a changelog entry](https://ukaea.github.io/PROCESS/development/versioning/).
* Finally, submit a pull request onto the `main` branch:
    * Link the relevant issue to the pull request.
    * Assign the pull request to a maintainer of the code that will have the correct expertise to review the change.

When making code contributions, we strongly recommend using pre-commit to verify your changes conform to our style requirements, otherwise the pull request will fail the 'quality' section of our GitHub actions. We document how this project uses pre-commit [here](https://ukaea.github.io/PROCESS/development/pre-commit/).

Please remember that all contributions are made under the [MIT license](https://github.com/ukaea/PROCESS/blob/main/LICENSE.txt).

### Testing
PROCESS has unit, integration, and regression tests. Any new functionality must be appropriately tested. Sometimes, changes may require other tests to be changed. These changes should be justified in the pull request description. Tests can be run locally by following [our testing documentation](https://ukaea.github.io/PROCESS/development/testing/). All pull requests will also be run against our GitHub actions which will run all of the tests and report back to the reviewer any failures. **The unit and integration tests must pass on the CI for the changes to be accepted**.

Regression tests, due to the nature of PROCESS, can change as model changes affect the optima which PROCESS converges to. A reviewer will review these changes to ensure they are minor and justified. We recommend justifying how a regression test is changing in the pull request discussion, a reviewer will likely request this anyway. For convenience, the CI system runs a 0% tolerance job that will highlight all differences between the current version of PROCESS on the `main` branch and your modified version of PROCESS; the 5% job excludes all differences under 5% differences between the two versions.

