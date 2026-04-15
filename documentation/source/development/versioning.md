# Versioning

## Semantic Versioning
Process attempts to use [semantic versioning](https://semver.org/), which takes the form `MAJOR.MINOR.PATCH`. Increment the:

1. MAJOR version when you make incompatible API changes,
2. MINOR version when you add functionality in a backwards compatible manner, and
3. PATCH version when you make backwards compatible bug fixes.

## Tags
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

!!! Warning "PROCESS Version 3"

PROCESS is still undergoing a significant restructure and, as such, PROCESS version 3 is unstable and does not guarantee backward compatibility. PROCESS version 4 will be the first major version to enforce backward-compatible API changes and will be released following a refactor of the data structure.