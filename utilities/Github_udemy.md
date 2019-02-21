# INTRODUCTION AND CORE CONCEPTS

Git is decentralized and works mainly offline, making it extremely fast.

The **repository** is the collection of files managed by Git.

Git works by saving the current state of all its files it manages into snapshopts called **commits**. A commit can contain one or mant more file changes (Git version files, not folders).

As you make changes, commits are saved onto the timeline, known as **branch**.

**GitHub** is a repository hosting service provided bt GitHub Inc.

# Git CONFIGURATION

Once installed, setup username and email (using same ones as in Github):

```{bash}
git config --global user.name "idelvalle"
git config --global user.email "ignacio.delvalle.torres@gmail.com"
git config --list # To show global config
```

We can also add VS code as default editor for git:

```{bash}
git config --global core.editor "code --wait"
git config --global -e # To test it
```


Is also recommended to install a graphical tool  for visual comparison and merging (e.g P4Merge: Visual Merge Tool)

```{bash}
git config --global diff.tool p4merge
git config --global difftool.prompt false
git config --global merge.tool p4merge
git config --global difftool.prompt false
```

# THE BASICS

## Initialization

To initiate a repository we use

```{bash}
git init . # The . specifies current folder
```


## Git States

**Locally** There are:

* The working Directory

* The staging area - Used to prepare for the next commit

* Git repository or commit history - `.git folder` Contains all committed or saved changes to Git repository

There is also a **Remote** state


## First Commit

The `git status` command tells you the situation in the repository

We add files using the `add` command, and from there to commit history by using `commit`

```
git add README.md # Moves README.md from wd into staging area
git commit -m "First file in demo repo" # To add comment -m"
```

## Repository

The `.git` folder is a special directory that Git manages internally

We can add comments to `commit` by typing `git commit`

## Commit details with Log and Show

`git log` Shows our listing of commits

`git show` shows the last commit and a **diff** containing all the changes. Press `q` to exit.

* To tell what files Git is tracking we use 

`git ls-files`

## Express Commit

We use the `git commit` command passing the `-a` parameter

* The `-a` parameter tells git tells Git to first add modified files to the Git staging area and then directly proceed to commiting.

`git commit -am` # The `m` allows us to add comment


## Backing Out Changes

If we add changes to staging area and we want to revert them we can use the `git reset` command:

`git reset HEAD <fileName>`

In this case, the file has been **unstaged**

If we can to revert all changes entirely (even deleting original file changes) we use `git checkout`

`git checkout -- <fileName> 


## History and Making New Commands with `alias`

A Git alias is used to shorten a command down

* `git log --oneline` will provide a simplified commit entry

* `git log --graph` will provide an asterisk based graph denoting our branching directory

* `git log --decorate` will provide us which commits are part of which branches

* `git log --all` will provide the history for all of the branches available in the repository

We can also create an **alias**, a shorter version of an existing command:

`git config --global alias.hist "log --oneline --graph --decorate --all"`

We can confirm by listing all git config entries:

`git config --global --list`

And test it by:

`git hist`

## Rename and Delete Files

Once a file has been commited,  can be renamed by using `git mv example.txt newname.txt`

**This will change ALSO the original filename**

This change is in the staging are and has to be commited

To delete a file:

`git rm <filename>` This removes it from our wd, but we have to commit it to remove it from the repository:

`git commit -m "removing <filename>"`

## Managing Files Outside Git

If we modify or delete files outside Git (e.g - delete or rename a file in our wd) we can use 

`git add -u` To stage deletions

`git add -A` To stage all types of possible modifications in the wd and accordingly update them in the git index

These changes should be **commited**


## Excluding Unwanted Files

Git provides a facility for excluding files and folders we don't want in our repository: a`.gitignore **file**`

The syntax for this file is only one pattern or expression per line

```
code .gitignore
*.log #Excludes every file ending in .log
git status
git add .gitignore
git commit -m "Adding ignore file"
```


# BEYOND THE BASICS

## Comparing Differences

We compare differences with the `diff` command

If we make a modification to a file and is in the staging area:

`git diff`

`git difftool` # Using p4merge

OR:

```
git hist
git diff c0c8183 HEAD # Difference between last commit on current branch and that commit
git diftool c0c8183 HEAD
```

## Branching and Merge Types

A branch is a timeline of commits; branches are the names or labels we give timelines in Git.

* **Fast-Forward Merges**: Simplest cases, when no additional work has been detected on the parental branch (master)

* **Automatic Merges**: Git detects non-conflicting changes in the parent branch. A new merge commit is created to show the merging of the two branches.

* **Manual Merge**: Happens when Git is unable to automatically resolve any conflicts. Git enters a special conflicting merge state, all merge conflicts must be resolved prior to moving forward with a commit.


## Special Markers

Git has special markers or pointers:

* **HEAD**: Points to last commit of current branch. Can also be moved.


## Branching Example

` git branch -a` shows us the branches. We can create a new branch to the master branch by : 

`git checkout -b <branchName>`

By doing so, we create a new branch called "branchname". If there are pending modifications, it carries those modifications forward into the new branch.

Once we commit the modifications on the new branch, we can observe differences between two branches by:

`git diff <master_branch> <branchName>`

To integrate changes on master branch we have to switch first to main branch (master):

`git checkout master` We switch branches by using `checkout` command

We then merge changes:

`git merge <branchName>`

Once we have updated the master branch we can delete the new branch by:

`git branch -d <branchName>`

## Conflict Resolution

We can, for example, create a new branch, update a file and modify the same file in the master branch in the same position without commiting it.

Wen we merge branches, we obtainan error. In this case, the auto-merging is not able to resolve the conflict. We have entered the **merging state**. If we `cat` the file implicated in the conflict, we can see the problem. To solve it we can open our mergetool:

`git mergetool`

Then, p4merge opens with all the possible versions of the file. We can select the version we want (select shape in bottom panel) and save it.

To complete the merge we need to commit what we have saved:

`git commit -m "Resolving Conflict"`

An `.orig` file appear, which we can add to `.gitignore` 

The original `.orig` file can then be deleted from our wd

## Marking Special Events with Tagging

Tags are **labels** you can put at any arbitrary commit point. There are two types:

* **Lightweight tags**: you can just give a name. There's no associated information with it. You create them with a name by:

    `git tag mytag`
    

And to see the list of tags: 

    `git tag --list` # To see the list

We can delete a tag by:

`git tag -d mytag`
    
* **Annotated tags**: they have information associated with it.
  

` git tag -a v1.0 -m "Release 1.0"` We also associate a commit message with the new tag `v1.0`

We can see the information associated with the tag:

`git show v1.0`

## Saving Work in Progress with Stashing

For example, we modify a file but we do not add it or commit it, we can **stash** it:

`git stash`

`git stash list`

If we come back to `git status`, we see we come back on a clean working directory

If we want to come back to our stash:

`git stash pop` This applies our last stash and drops the stash that was applied.

We can then commit that changes.

## Time Travel with Reset and Reflog

To "come back" and modify a previous commit we can use `reset`

`git reset b7b12cb --soft`

There are three ways of resetting:

* **Soft resetting**: Least destructive. All it does is change where HEAD is pointing. Preserves the staging area and working directory. We can back out our changes, make minir modifications to them, and commit where head is currently pointing.

`git reset b7b12cb --soft`

* **Mixed**: Default, unstages the number of changes.

`git reset 2a56a98 --mixed` 

* **Hard**: The most destructive of all reset modes. The working directory is now clean, any pending changes have been wiped out, along with anything that was in the staging area.

` git reset 7b98a87 --hard`

`git log` shows us our commit ids , and `git reflog` shows us all different actions taken while in the repository. We can use it to move between commits

# WELCOME TO GitHub

## Creating Repositories

Is the most popular hosting service for Git repositories.

* Click on green button "create new repository"

* Add repository name (demo)

* Add description

* Press Create Repository

To link (push) an existing repository from the command line, navigate to "demo" Git repository and:

`git remote add origin https://github.com/idelvalle/demo.git` 

`remote` sets up connections; can be listed by using:

`git remote -v`

The `add` command takes two parameters:

* Name of remote reference we want to create ("origin"). ** by convention, the first and primary remote repository is named origin**

* Full URL to remote repository

`git remote -v` shows two URLs: one for fetching and one for pushing.

To synchronize changes between local and remote Git repositories we use `push` command, to push commits made on your local branch to a remote repository:

The `git push` command takes two arguments:

A `remote name`, for example, origin
A `branch name`, for example, master

`git push -u origin master --tags`

When we do pushes in the future we will not need `-u`

`--tags` sends all tags currently in our Git repository to GitHub


# SSH AUTHENTICATION

## Overview

Takes time to set up, but helps to save it. Do it in own computer or in the one you use regularly.

## SSH vs HTTPS

Avoids to introduce username and password each time we want to `push` into repository:

`git push origin master`

## Generating an SSH Key

* From home directory, look for `.ssh` directory

* Create it if we do not have it: `mkdir .ssh` and move in to it `cd .ssh`

* Create key `ssh-keygen -t rsa -C "ignacio.delvalle.torres@gmail.com"`

* Press enter to save in same file and add passphrase.

* List all files with `ls -al`

* The `.pub` file is the public key and the other newly generated file is the private key. Open the `.pub` file with editor and select and copy everything.

## Verify SSH Authentication with GitHub

* Add ssh key from clipboards to GitHub profile in personal settings section. 

* Come back to terminal, inside the `.ssh` folder and type:

`ssh -T git@github.com`

* Introduce yes and confirm


# GitHub REPOSITORY

## Starting Remote

Create repository online:

* Give name

* Initialize with a README

* Add .gitignore: Node

* Add icense: Apache 2.0


## Create a Local Copy with Clone

* Select Clone with SSH and copy to clipboard

* Go to desired folder in terimanl and type:

* `git clone git@github.com:idelvalle/<name of repository.git>`

** By default GitHub will clone the repository into a project flder named under the repository name**, we can also specify it if we want to change it:

* `git clone git@github.com:idelvalle/<name of repository.git> <foldername>`


## Seeding the Repository with Sample Content

* `cp -R ~/Downloads/initializr/* .`

* `git status`

* `git add .`

* `git status`

* `git push origin master` origin: name of remote to use master: branch to push up


## Fetch and Pull

If we make changes on our remote repository (GitHub) we have to update our local repository. If then we update a file in our local repository and try to push it we get an error message.

* **Pull** is two commands in one; it is a **fetch** and also a **merge**. Git will first fetch all the updates from the remote repository, and then it will merge those changes into our current repository. Fetches are safer than pulls, because they are non-destructive.

The order of commands would be:

* `git fetch`

* `git status`

* `git pull` Accept or edit Merge file

* `git push`

## Repository Features and Settings

We can delete the repository under settings in GitHub

The "download zip"" option is useful when we want to pass our code to someone not in GitHub.


## Updating Repository and Remote References

If we rename our repository in GitHub under settings we need to update it. Any clones from this repository need tobe updated so they will continue to work.

When we do `git remote -v` we see the "old name" of our repository.

In order to update the reference to origin:

* Copy url to clipboard

* `git remote set-url origin git@github.com:idelvalle/website.git` Name of reference we need to update: origin

* Confirm update by `git remote -v`

* We can get additional information about remote reference (origin) by `git remote show origin`


## Directly Editing Files on GitHub

We can edit files directly on GitHub (not really recommended): select a file, edit it and commit changes (either directly to master branch or in a new branch)


## Creating New Files on GitHub

We can create new files. 

Once created we will create a new branch for the commit and start a pull request:

* Rename the branch 

* Propose new file

* Create pull request commit message

* Create pull request

* Merge Pull Request

* Confirm Merge

* Delete Branch

## Create a New File on Master Branch

In the same way we edited files we can also create a new file and then commit it.


## Renaming and Deleting Files on GitHub

We can edit and delete files (trash symbol)


## Synchronizing Changes with Local Repository

If me make changes in our remote repository should be synchronized with our local one 

`git fetch`

`git status` Shows us that we are behind origin/master. 

`git pull` To update local branch


## Reviewing Commits

We can click on `commits` and see our commit list ordered by date; we can also copy the full SHA-1 if we want.

We can also leave comments (plus mark) and edit the comment. This allow us make comments on other people's code.

## GitHub Time Travel

Under commit section, we can click on `<>` and browse the repository at that point in the history.

Instead of pointing to a specific branch (master in our case) we are now pointing to `Tree` and then the **commit id**

We can also `browse the files` at that point of history.

## Using Commit IDs with Local Repository

Under the commits tab, we can click on the `clipboard icon` which will copy the full SHA-1 commit id to my clipboard.

Then we head over to out local Git repository and we use `show` command to get additional information about the commit id we just copied: 

` git show 44336e058ad33f2f4bf02a7360f707bdf4fad8e7`

This should show the commit id listed, along with the author, date and commit message, as well as a diff listing for all the changes that are part of the commit.


# GitHub REPOSITORY BRANCHES

## Creating Branches on GitHub

We can click on `switch branches` and if the branch oes not exist, it will be created. After pressing `enter` we will automatically switch into the new branch (example)


## Local Branches

It is usually best practice to create a branch locally and then synchronize it with GitHub:

`git checkout -b remove-ipsum` Create a new Branch locally and switch into it at same time

`git rm lipsum.txt` Remove a file

`git status` We see the change is staged

`git commit -m "Removing ipsum file"` Commit changes

`git push -u origin remove-ipsum` Push changes with new branch. The new branch is being pushed up into GitHub.

If we go to GitHub and we refresh we can see the changes


## Comparing and Pull Requests in GitHub

We can click the green button `Compare & pull request` in order to merge changes and delete the branch afterwards.


## Merging Locally

We can also merge changes locally from Git. In order to integrate our changes into master, we need to change back to the master branch






# TAGS AND RELEASES

Locally we use tags to mark important events.

## Local Tags

`git log --oneline` Give us a summary of recent commits in our Git repository.


`git tag` Without parameters gives us a list of existing tags

`git tag unstable develop`  `develop` is the reference point to associate the tag

We can also create an annotated tag:

`git tag -a v0.1-alpha -m "Release 0.1 (Alpha)" f25924a`

With annotated tags we have the option of using the `git show` command specifying the tag: 


`git show v0.1-alpha`
















​    
​    






















































































































































