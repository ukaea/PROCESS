# Contributing to `PROCESS`
There are many valuable contributions that can be made to PROCESS:
* Reporting bugs.
* Requesting/ recommending features.
* Implementing features or bugs raised as issues.
* Updating and improving documentation.

When contributing to this repository, please first discuss the change you wish to make via [issue](https://github.com/ukaea/PROCESS/issues), [discussion](https://github.com/ukaea/PROCESS/discussions), [contact form](https://www.ukaea.org/service/process/), or any other method with the owners of this repository before making a change.

Please remember that all contributions and communication regarding PROCESS are subject to our [Code of Conduct](https://github.com/ukaea/PROCESS/blob/main/CODE_OF_CONDUCT.md).

## Creating an issue
Issues can be used to report bugs or request features and improvements. We ask you help us manage our issues by:
* Verifying your issue is not a duplicate of another issue; in this case, we welcome your contribution as a reply to the issue.
* Ensure you completely describe the bug/feature/problem and complete the given templates with appropriate detail.

## Submitting a pull request
Please discuss any feature ideas you have with the developers before submitting them, as you may not be aware of parallel development work taking place, or implementation decisions / subtleties which are relevant. The ideal workflow for submitting a pull request is as follows:

1. Discuss the feature with a core PROCESS developer (as mentioned above). This will normally involve creating an issue with the bug/feature--please indicate your willingness to contribute at the end of the issue (e.g. "I propose we do x to fix this bug, I am happy to create a PR"). If an issue already exists, please comment indicating your proposed changes and a willingness to contribute. A PROCESS maintainer will assign the issue when they are happy for your to start development and create a PR. 
2. Fork our repository and create a branch off of `main` where you will make your changes. Members of the PROCESS development team do not need to complete this step, but should ensure they create an appropriately named branch.
3. Make the relevant changes for the repository (ensuring the changes do not creep away from the scope of the issue). Please ensure that the implementation follows the PROCESS [style guide](https://ukaea.github.io/PROCESS/development/standards/) for names of variables and functions etc. Ensure that as part of your pull request you include tests for your feature/change and any updates/additions to the documentation that would be appropriate.
4. Discuss any problems or development choices in the issue and keep everyone updated on the progress of your development.
5. Finally, submit a pull request onto the `main` branch:
    * Please indicate which issue(s) this PR closes by writing `Closes #<issue number>` at the top of the PR.
    * Provide an overview of the changes you have made, including a list of new/renamed files, functions, and variables. 
    * Please highlight and justify any changes you have made to the test suite as these will need to be scrutinised. 
    * You may also find it appropriate to discuss any development/design choices you have made, including things you have chosen NOT to do. 

When making code contributions, we strongly recommend using pre-commit to verify your changes conform to our style requirements, otherwise the pull request will fail the 'quality' section of our GitHub actions. We document how this project uses pre-commit [here](https://ukaea.github.io/PROCESS/development/pre-commit/).

Please remember that all contributions are made under the [MIT license](https://github.com/ukaea/PROCESS/blob/main/LICENSE.txt).

### Reviewing a pull request
There are three parties involved in a pull request:
**The reviewers**
The reviewers are responsible for:
* Ensuring that the pull request implements the feature or fixes the bug from the issue(s) it is closing.
* Verifying that the contributions are semantically correct and of an appropriate quality.
* Ensuring that a pull request is complete (e.g. not missing tests or documentation updates).

This will often involve making comments on specific areas of the code or asking general questions about the contributions. Reviewers will often review the pull request several times as the author makes changes based upon their comments. If a reviewer "requests changes" they must approve the pull request before it can be merged. 

A pull request must be reviewed by at least one of the codeowners (@ukaea/process-maintainers) before it can be merged. Codeowners are the maintainers of (a section of) the PROCESS repository. Their review will be requested automatically on each pull request. Other reviewers can still be added by maintainers/developers if they believe the pull request will benefit from that specific person's review.

**The assignee**
This is the 'lead reviewer' of the pull request and will be one of the codeowners. Codeowners will assign themselves/eachother to a pull request based upon availability and expertise. The assignee is ultimately responsible for ensuring a pull request is ready to be merged:
* The changes are appropriately tested and documented.
* The changes pass all existing and new tests (see below for more detail about tests that can 'fail').
* The pull request fully implements the feature/bugfix from the issue.
* The code is of sufficient quality.
* The changes will not conflict with recent/upcoming merges onto main.

The assignee will be the person to merge the pull request onto main. They may also need to create new issues (e.g. for future work that arises from the pull request).

**The author**
The author is the user that creates the pull request after being assigned to an issue and making changes on a branch. They are responsible for ensuring they make correct and high-quality contributions to PROCESS. This will engaging with the reviewers and assignee in the comments of the pull request and ammending/correcting the pull request according to comments from the reviewers. 

### Testing
PROCESS has unit, integration, and regression tests. Any new functionality must be appropriately tested. Sometimes, changes may require other tests to be changed. These changes should be justified in the pull request description. Tests can be run locally by following [our testing documentation](https://ukaea.github.io/PROCESS/development/testing/). All pull requests will also be run against our GitHub actions which will run all of the tests and report back to the reviewer any failures. **The unit and integration tests must pass on the CI for the changes to be accepted**.

Regression tests, due to the nature of PROCESS, can change as model changes affect the optima which PROCESS converges to. A reviewer will review these changes to ensure they are minor and justified. We recommend justifying how a regression test is changing in the pull request discussion, a reviewer will likely request this anyway. For convenience, the CI system runs a 0.2% tolerance job that will highlight all differences between the current version of PROCESS on the `main` branch and your modified version of PROCESS; the 5% job excludes all differences under 5% differences between the two versions.

