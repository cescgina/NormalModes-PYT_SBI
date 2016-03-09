## NormalModes-PYT_SBI

.. contents::

### AUTHORS

    David Mas Ponte(@davidmasp) & Joan Francesc Gilabert Navarro(@cescgina)

### GIT USAGE

1) You can propose changes (add it to the Index) using

...`git add <filename>`
...`git add *`

2) This is the first step in the basic git workflow. To actually commit these changes use
  `git commit -m "Commit message"`
Now the file is committed to the HEAD, but not in your remote repository yet.

3) Your changes are now in the HEAD of your local working copy. To send those changes to your remote repository, execute
  `git push origin master`

...B) Change master to whatever branch you want to push your changes to.
create a new branch named "feature_x" and switch to it using

...`git checkout -b feature_x`

...switch back to master
...`git checkout master`
...and delete the branch again
...`git branch -d feature_x`
...a branch is not available to others unless you push the branch to your remote repository
...`git push origin <branch>`

...C) update & merge
...to update your local repository to the newest commit, execute
...`git pull`
...in your working directory to fetch and merge remote changes.
...to merge another branch into your active branch (e.g. master), use
...`git merge <branch>`
...in both cases git tries to auto-merge changes. Unfortunately, this is not always possible and results in conflicts. You are responsible to merge those conflicts manually by editing the files shown by git. After changing, you need to mark them as merged with
...`git add <filename>`
...before merging changes, you can also preview them by using
...`git diff <source_branch> <target_branch>`

### TO-DO

* Write header and trailer to superimposed PDB file (it seems like it is not
  easy, maybe we should not write the new coordinates)
* Add possibility of doing the essential dynamics from MD trajectories
