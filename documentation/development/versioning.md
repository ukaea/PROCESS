# Versioning

## Semantic Versioning
Process attempts to use [semantic versioning](https://semver.org/), which takes the form `MAJOR.MINOR.PATCH`. Increment the:

1. MAJOR version when you make incompatible API changes,
2. MINOR version when you add functionality in a backwards compatible manner, and
3. PATCH version when you make backwards compatible bug fixes.

## Upversioning
Process uses a `develop` and `master` release model. All development work is merged into `develop`, whereas users of the code use `master`. When a new version is required to be released to users, the version is incremented on `develop`, `develop` is merged into `master` and the merge commit tagged with the version number on `master`. Users will always use the latest tagged version on `master`.

To upversion Process:

1. Create a branch from `main` and merge request for the upversion
2. Increment the version in `setup.py` according to the semver rules
3. Similarly increment the version in `main_module.f90`:`inform()`. This Fortran version step will be removed in time as the Python conversion progresses
4. Push and merge into `main`
5. Create a release through GitHub and create a tag according to semantic versioning

### Extracting merge commits from git log
A useful command for getting merge commits from git log for writing the changelog is:
```bash
git log v2.1..HEAD --merges --first-parent develop --pretty=format:%Cblue%B%n >> mergeCommits.log
```
`v2.1..HEAD` is the log range: from the last tag (`v2.1`) to `HEAD`

`--merges` filters for merge commits

`--first-parent develop` filters for the branch you're merging into (`develop`): i.e. only include merges into `develop`

`%B` is subject and body together: full merge commit content

`>> mergeCommits.log` output to file (optional)

## Tags
### Tagging a commit
```
git tag -a vX.Y.Z -m "Version X.Y.Z"
```
`git tag -a vX.Y.Z` tags the current commit with an annotated (`-a`) tag of "vX.Y.Z". A message is specified with `-m`. `git show vX.Y.Z` will show the commit that has just been tagged.

The commit is only tagged locally, and `git push` won't transfer the tag to the remote.
```
git push origin vX.Y.Z
```
is required to transfer the tag to the remote.

### Working with tags
| Command               | Description                    |
| --------------------- | ------------------------------ |
| `git describe`        | show the current tag           |
| `git tag`             | list all tags                  |
| `git tag -l "1.0.*"`  | list tags contained in `1.0.z` |
| `git checkout vX.Y.Z` | checkout a specific tag        |

Between user tags `git` will create tags in the following format:
```
1.0.12-11-g3f1b433
```

- `1.0.12` is the last manually entered tag by the user
- `11` is the number of commits since that tag
- `g3f1b433` is a unique identifier for this specific commit

This allows the user to checkout a specific commit between tagged versions. PROCESS now outputs this information into the `OUT.DAT` and `MFILE.DAT` and is 
updated upon compilation. This way each output file is trackable to a specific commit.