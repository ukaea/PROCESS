# git Usage
## Commit Logs
To see the commit messages you can use the `git log` command. There are various options
described below.

| Command | Description | Example |
| -------- | -------- | --------- |
| `-(n)`   | show the last `n` commits  | `git log -5` |
| `--since or --after` |  limits the logs to be from date given | `git log --since "21-01-15"` |
| | can use `number.scale` where scale=year, month, week, day and minute | `git log --since 2.weeks` |
| `--until or --before` |  limits the logs to be up to date given | `git log --until "22-01-15"` |
| `--author` | only shows commits from given author | `git log --author "morrisj"` |
| `--grep` | only show commits with a commit message containing the string given | `git log --grep "magnet"`  |
| `--stat` | if you want to see some abbreviated stats for each commit | `git log --stat` |
| `--oneline` | Outputs commit number, date and message to a single line | `git log --oneline` |
| `--graph` | display commits in a ASCI graph/timeline | `git log --graph` |
| `-S` | only show commits adding or removing code matching the string | `git log -S "find_me"` |

!!! Info "Log --> File"
    to output the log to a file add `>> file_name.log` to the end of the command

## Changing the code

- Create an issue on the PROCESS GitHub repo
- Create a branch associated with that issue on GitHub 
- Clone the repository `git clone https://github.com/ukaea/PROCESS`
- Swap to your branch `git checkout [new_branch_name]`

## Committing changes

Make your changes to the code and at suitable stages commit locally:

  - `git add file_changed_1 file_changed_2`
  - `git commit -m "COMMIT MESSAGE"`

The commit message should be informative and give useful information for future development. Such as:

```
Made changes to the TF coil magnet model. 
  - Updated the allowable stress in the coils to be 600MPa. 
  - Remove side-wall case. Ran test suite and everything OK.

```  
not

```
Update to magnet model
```

When you wish to push your branch back to the repository enter `git push`

## Merging
### Develop into your branch

1. Make sure you have committed all of your changes to your local branch.
2. Update your local repo with `git pull`
3. Checkout the development branch `git checkout develop`
4. Check remote repo again `git pull`
5. Checkout your new branch `git checkout my_branch_name`
6. Merge develop into your branch `git merge develop`
7. If there are conflicts check the files listed for the following:
  ```
  This line was edited in develop branch
  ```
8. Resolve any conflicts then `git add file_1 file_2` where file_1 and file_2 are files that had conflicts.
9. Commit the changes `git commit`
10. Push the branch back to the remote repo `git push`

### Your branch into develop

1. Check your repo is up to date `git pull`
2. `git checkout my_branch_name`
3. `git checkout develop`
4. `git merge my_branch_name`
5. Resolve conflicts in similar manner to section above
6. `git push`