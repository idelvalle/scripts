---
title: "Linux COMMAND LINE AND SHELL SCRIPTING BIBLE"
output: html_notebook
---

# STARTING WITH Linux SHELLS

Four main parts make up a Linux system:

1. The Linux kernel

**The core of the Linux system is the kernel**. The kernel controls all the hardware and software on the computer system, allocating hardware when necessary and executing software when required.

Linus is the person responsible for creating the first Linux kernel software when he was a student at the University of Helsinki. He intended it to be a copy of the Unix system, at the time a popular operating system used at many universities.

The kernel is primarily responsible for four main functions:

* System memory management

One of the primary functions of the operating system kernel is memory management. Not only does the kernel manage the physical memory available on the server, but it can also create and manage virtual memory, or memory that does not actually exist.
It does this by using space on the hard disk, called the swap space.

* Software program management

The Linux operating system calls a running program a process.
The kernel creates the first process, called the init process, to start all other processes on the system.

* Hardware management

Any device that the Linux system must communicate with needs driver code inserted inside the kernel code.
The driver code allows the kernel to pass data back and forth to the device, acting as a middle man between applications and the hardware.

* Filesystem management

The Linux kernel can support different types of filesystems to read and write data to and from hard drives.


2. The GNU utilities

The GNU organization (GNU stands for GNU’s Not Unix) developed a complete set of Unix utilities
The GNU/Linux shell is a special interactive utility. It provides a way for users to start programs, manage files on the filesystem, and manage processes running on the Linux system. The core of the shell is the command prompt, which is the interactive part of the shell.
The default shell used in all Linux distributions is the bash shell. The bash shell was developed by the GNU project as a replacement for the standard Unix shell, called the Bourne shell

3. A graphical desktop environment

There are a plethora of graphical desktops you can choose from in Linux (KDE, GNOME, Unity-UBUNTU).

4. Application software

Each of these parts has a specific job in the Linux system. No part is very useful by itself.

A complete Linux system package is called a **distribution**. Most distributions are customized for a specific user group, such as business users, multimedia enthusiasts, software developers, or average home users.


# GETTING TO THE SHELL

Before the days of graphical desktops, the only way to interact with a Unix system was through a text command line interface (CLI) provided by the shell.


Often, you needed only a simple dumb terminal to interact with the Unix system. A dumb terminal was usually nothing more than a monitor and keyboard connected to the Unix system via a communication cable (usually a multi-wire serial cable). This simple combination provided an easy way to enter text data into the Unix system and view text results.

The alternative to using a virtual console terminal is to use a terminal emulation package from within the Linux graphical desktop environment. A terminal emulation package simulates working on a console terminal, but within a desktop graphical window.

Although many graphical terminal emulator packages are available, often are installed in Linux distributions by default: they are GNOME Terminal, Konsole Terminal, and xterm.


# BASIC SHELL COMMANDS

The GNU bash shell is a program that provides interactive access to the Linux system.

The /etc/passwd file contains a list of all the system user accounts, along with some basic configuration information about each user.


## The prompt

After you start a terminal emulation package or log in to a Linux virtual console, you get access to the shell CLI prompt. The prompt is your gateway to the shell. This is the place where you enter shell commands.

The default prompt symbol for the bash shell is the dollar sign (```$``` ). This symbol indicates that the shell is waiting for you to enter text.

## The bash manual

The ```man``` command provides access to the manual pages stored on the Linux system.

When you are finished with the man pages, press the ```q``` key to quit.

You can search the man pages using **keywords**. The syntax is ```man -k keyword```. For example, to find commands
dealing with the terminals, you type ```man -k terminal```.

## Navigating File System

The two special characters used for relative directory references are:

* The single dot (```.```) to represent the current directory

* The double dot (```..```) to represent the parent directory


## Listing Files and Directories


The ```-F``` parameter flags the directories with a forward slash (```/```), to help identify them in the
listing.

```{bash}
ls -F
```

To display hidden files along with normal files and directories, use the ```-a``` parameter.

```{bash}
ls -F -R
```
The ```-R``` parameter is another option the ls command can use. Called the recursive option, it shows files that are contained within subdirectories in the current directory.

Option parameters don’t have to be entered separately as shown in the nearby example: ```ls -F -R```. They can often be combined as follows: ```ls -FR```

For listing additional information, another popular parameter is ```-l``` . The ```-l``` parameter produces a long listing format, providing more information about each file in the directory:

The ```ls``` command also recognizes standard *wildcard* characters and uses them to match patterns within the filter:

The ```-d``` option lists a directory’s information but not its contents

* A question mark ```?``` to represent **one** character

* An asterisk ```*``` to represent **any number** of characters (zero or more)

Using the asterisk and question mark in the filter is called **file globbing**. 

File globbing is the processing of pattern matching using wildcards. The wildcards are officially called metacharacter wildcards.

You can list choices of characters and you can specify a range of characters, such as an alphabetic range ```[a - i]```

Also, you can specify what should not be included in the pattern match by using the exclamation point ```!```


## Creating Files

The ```touch``` command creates the new file you specify and assigns your username as the file owner.


## Copying Files

In its most basic form, the ```cp``` command uses two parameters — the source object and thedestination object: ```cp source destination```.

When both the source and destination parameters are filenames, the ```cp``` command copies the source file to a new destination file.

It is best to add the ```-i``` option to force the shell to ask whether you want to overwrite a file

The ```-R``` parameter is a powerful cp command option. It allows you to recursively copy the contents of an entire directory in one command.


## Linking Files

If you need to maintain two (or more) copies of the same file on the system, instead of having separate physical
copies, you can use one physical copy and multiple virtual copies, called **links**.

A **symbolic link*** is simply a physical file that points to another file somewhere in the virtual directory structure. The two symbolically linked together files do not share the same contents.

To create a symbolic link to a file, the original file must pre-exist. We can then use the ```ln``` command with the ```-s``` option to create the symbolic link

Another way to tell that these linked files are separate physical files is by viewing their **inode** number. The inode number of a file or directory is a unique identification number that the kernel assigns to each object in the filesystem. 
To view a file or directory’s inode number, add the ```-i``` parameter to the ls command:

```ls -i <filename>```

A **hard link** creates a separate virtual file that contains information about the original file and where to locate it. However, they are physically the same file.

No parameter is needed on the ```ln``` command:

```ln <filename> <hl_filename>```

When we check the inode numbers, both files share the same inode number. This is because they are physically the same file.

## Renaming Files

Renaming files is called **moving** files using the ```mv``` command.

## Deleting Files

The command to remove files in the bash shell is ```rm```. A good habit is to always tack on the ```-i``` parameter to the ```rm``` command.


## Creating Directories

Creating a new directory in Linux is easy — just use the ````mkdir``` command

**To create several directories and subdirectories at the same time, you need to add the ```-p``` parameter**:

```mkdir -p New_Dir/Sub_Dir/Under_Dir```

The ```-p``` option on the ```mkdir``` command makes any missing **parent directories** as needed.

## Removing Directories

The basic command for removing a directory is ```rmdir```

By default, the ```rmdir``` command works only for removing **empty** directories

You can also use the ```rm``` command on entire non-empty directories. Using the ```-r``` option allows the command to descend into the directory, remove the files, and then remove the directory itself.

The ```tree``` utility nicely displays directories, subdirectories, and their files


## Viewing the file type

The ```file``` command can peek inside of a file and determine just what kind of file it is:


## Viewing the whole file

Linux has three different commands tosee text files:

**1 ```cat``` command**

The ```-n``` parameter numbers all the lines for you

If you just want to number the lines that have text in them, use the ```-b``` parameter.

For large files, the cat command can be somewhat annoying. The text in the file just quickly scrolls off the display without stopping.


**2 ```more``` command**

The ```more``` command displays a text file, but stops after it displays each page of data.

You can use more to navigate through a text file by pressing the ```spacebar``` or you can go forward line by line using the ```Enter``` key. When you are finished navigating through the file using ```more``` , press the ```q``` key to quit.

**3 ```less``` command**

Is an advanced version of the ```more``` command (less is more)

Supports the same command set as the more command, plus many more options


## Viewing file parts

### ```tail``` command

By default, it shows the last 10 lines in the file.

You can change the number of lines shown using ```tail``` by including the ```-n``` parameter.

```tail -n 2 <filename>```

The ```-f``` parameter is a pretty cool feature of the tail command. It allows you to peek inside a file as the file is being used by other processes.

### ```head``` command

By default, it displays the first 10 lines of text.


## Monitoring Programs

When a program runs on the system, it’s referred to as a **process**.

By default, the `ps` command shows only the processes that belong to the current user and that are running on the current terminal.

```{bash}
ps
```


The basic output shows the process ID (`PID`) of the programs, the terminal (`TTY`) that they are running from, and the `CPU` time the process has used.


The GNU `ps` command that’s used in Linux systems supports three different types of command line parameters

The key to using the ps command is not to memorize all the available parameters — only those you find most useful.

### Unix-style parameters (preceded by a dash)

For example, if you need to see everything running on the system, use the -ef parameter combination

`ps -ef`

* **UID**: The user responsible for launching the process

* **PID**: The process ID of the process

* **PPID**: The PID of the parent process (if a process is started by another process)

* **C**: Processor utilization over the lifetime of the process

* **STIME**: The system time when the process started

* **TTY**: The terminal device from which the process was launched

* **TIME**: The cumulative CPU time required to run the process

* **CMD**: The name of the program that was started

For even more information, you can use the -l parameter, which produces the long format output

* F: System flags assigned to the process by the kernel

* S: The state of the process ( O = running on processor; S = sleeping; R = runnable, waiting to run; Z = zombie, process terminated but parent not available; T = process stopped)

* PRI:The priority of the process (higher numbers mean lower priority)

* NI:The nice value, which is used for determining priorities

* ADDR:The memory address of the process

* SZ :Approximate amount of swap space required if the process was swapped out

* WCHAN :Address of the kernel function where the process is sleeping


### BSD-style parameters, which are not preceded by a dash

Unix and BSD types of parameters have lots of overlap. Most of the information you can get from one you can also get from the other.

```{bash}
ps l
```

while many of the output columns are the same as when we used the Unix-style parameters, some different ones appear as well:

* **VSZ**: The size in kilobytes of the process in memory

* **RSS**: The physical memory that a process has used that isn’t swapped out

* **STAT**: A two-character state code representing the current process state

In the STAT column the first character uses the same values as the Unix-style S output column, showing when a process is sleeping, running, or waiting. The second character further defines the process’s status:

* < : The process is running at high priority.

* N : The process is running at low priority.

* L : The process has pages locked in memory.

* s : The process is a session leader.

* l : The process is multi-threaded.

* + : The process is running in the foreground.

### GNU long parameters, which are preceded by a double dash

The `—forest parameter` displays the hierarchical process information, but using ASCII characters to draw charts:

`ps --forest`


### Real-time process monitoring

The `top` command displays process information similarly to the ps command, but it does it in real-time mode

The first section of the output shows general system information.

The **load average** appears as three numbers: the 1-minute, 5-minute, and 15-minute load averages. The higher the values, the more load the system is experiencing. **If the 15-minute load value is high, your system may be in trouble.**


### Stopping processes

In Linux, processes communicate with each other using **signals**. A process signal is a predefined message that processes recognize and may choose to ignore or act on.


**1. The `kill` command**

The `kill` command allows you to send signals to processes based on their process ID (PID). You can only use the process PID instead of its command name, making the kill command difficult to use sometimes.


**2. The `killall` command**

The `killall` command is a powerful way to stop processes by using their names rather than the PID numbers. The `killall` command allows you to use wildcard characters as well.


## Monitoring Disk Space


### Mounting media

The Linux filesystem combines all media disks into a single virtual directory. Before you can use a new media disk on your system, you must place it in the virtual directory. This task is called **mounting**.

The `mount` command displays a list of media devices currently mounted on the system providing:

* The device filename of the media

* The mount point in the virtual directory where the media is mounted

* The filesystem type

* The access status of the mounted media

The basic command for manually mounting a media device is:

`sudo mount -t type device directory`

To remove a removable media device, you should **unmount** it first with `umount`:

`umount [directory | device]`

The `umount` command gives you the choice of defining the media device by either its device location or its mounted directory name.


### Using the `df` command

The `df` command shows each mounted filesystem that contains data

The `-h` command shows the disk space in human-readable form, usually as an M for megabytes or a G for gigabytes

### Using the `du` command

The `du` command shows the disk usage for a specific directory (by default, the current directory).

The number at the left of each line is the number of disk blocks that each file or directory takes.


## Working with Data Files

### Sorting data

The `sort` command sorts the data lines in a text file using standard sorting rules for the language you specify as the default for the session.

The `-n` parameter tells the sort command to recognize numbers as numbers instead of characters and to sort them based on their numerical values


The `-M` is the month-sort parameter

The `-k` and `-t` parameters are handy when sorting data that uses fields. Use the `-t` parameter to specify the field separator character, and use the `-k` parameter to specify which field to sort on.


### Searching for data

The command line format for the `grep` command is:

`grep [options] pattern [file]`

If you want to reverse the search (output lines that don’t match the pattern), use the `-v` parameter

If you need to find the line numbers where the matching patterns are found, use the `-n` parameter

If you just need to see a count of how many lines contain the matching pattern, use the `-c` parameter

If you need to specify more than one matching pattern, use the `-e` parameter to specify each individual pattern:

`grep -e t -e f file1`


### Compressing for data

The gzip package is a creation of the GNU Project, in their attempt to create a free version of the original Unix compress utility. This package includes these files:

* **gzip** for compressing files

* **gzcat** for displaying the contents of compressed text files

* **gunzip** for uncompressing files


### Archiving data

The most popular archiving tool used in Unix and Linux is the `tar` command.

`tar function [options] object1 object2 ...`


### Using command aliases

An **external command**, sometimes called a filesystem command, is a program that exists outside of the bash shell. They are not built into the shell program. An external command program is typically located in /bin , /usr/bin , /sbin , or /usr/sbin .

The `ps` command is an external command. You can find its filename by using both the `which` and the `type` commands.

**Built-in commands** were compiled into the shell and thus are part of the shell’s toolkit. No external program file exists to run them.You can tell a command is built-in by using the type command:

```{bash}
type cd
```

A useful built-in command is the `history` command. The bash shell keeps track of the commands you have used.

To recall and reuse your last command, type `!!` and press the Enter key.

To see a list of the active aliases, use the alias command with the -p parameter:

`alias -p`


# USING Linux ENVIRONMENTAL VARIABLES

## Exploring Environment Variables

### Looking at global environmental variables

Global environment variables are visible from the shell session and from any spawned child subshells.

To view global environment variables, use the `env` or the `printenv` command:

You	can	also	use	the	 `echo` 	 command	 to	 display	 a	 variable’s	 value.	 When	 referencing	 an environment	variable	in	this	case,	you	must	place	a	dollar	sign	(`$`)	before	the	environment variable	name:

```{bash}
echo $HOME
```


### Looking at local environmental variables

Can	be	seen	only	in	the	local	process in	which	they	are	defined.

The `set` command	 displays	 all	 variables defined	for	a	specific	process,	including	both	local	and	global	environment	variables	and user-defined	variables.

The `set` command	displays	both	global	and	local	environment	variables	and	user-defined variables.	It	also	sorts	the	display	alphabetically.	The `env` and `printenv` are	different from `set` 	in	that	they	do	not	sort	the	variables,	nor	do	they	include	local	environment or	local	user-defined	variables.	

The	`env` command	has	additional	functionality	that `printenv` does	not	have,	making	it	the	slightly	more	powerful	command.


## Setting	User-Defined	Variables

### Setting	local	user-defined	variables

You	 can	 assign	 either	 a numeric	or	a	string	value	to	an	environment	variable	by	assigning	the	variable	to	a	value
using	the	equal	sign.

If	you	need	to	assign	a	string	value	that	contains	spaces,	you	need	to	use	a	single	or	double quotation	mark	to	delineate	the	beginning	and	the	end	of	the	string:

```{bash}
my_variable=“Hello	World”
echo	$my_variable
```


**The	standard	bash	shell	convention	is	for	all	environment	variables	to	use	uppercase letters**.

**It’s	 extremely	 important	 that	 you	 not	 use	 spaces	 between	 the	 variable	 name,	 the	 equal sign,	and	the	value.**

The	local	variable	set	within	the	child	shell	doesn’t	exist	after	a	return	to	the	parent	shell.
You	 can	 change	 this	 behavior	 by	 turning	 your	 local	 user-defined	 variable	 into	 a	 global environment	variable.

### Setting global environment	variables

The	method	used	to	create	a	global	environment	variable	is to	first	create	a	local	variable	and	then	export	it	to	the	global	environment.

This	is	done	by	using	the `export` command	and	the	variable	name	minus	the	dollar	sign:

`my_variable=“I	am	Global	now”`
`export	my_variable`


## Removing environment variables

You	can	do	this	with	the `unset` command.	When referencing	 the	 environment	 variable	 in	 the `unset`  command,	 remember	 not	to	 use	 the dollar	sign:

`echo	$my_variable`
`unset my_variable`

**If	you	are	doing	anything	with	the variable,	use	the	dollar	sign.	If	you	are	doing	anything	to	the	variable,	don’t	use	the
dollar	sign.**

## Setting	the	PATH	Environment	Variable

Reference	the original	PATH value and add any new directories	to	the	string.	

`PATH=$PATH:<path>`

## Locating	System	Environment	Variables

### Understanding	the	login	shell	process

When	you	log	in	to	the	Linux	system,	the	bash	shell	starts	as	a	login	shell.	The	login	shell typically	looks	for	five	different	startup	files	to	process	commands	from:
	
* **/etc/profile**

The	/etc/profile file	is	the	main	default	startup	file	for	the	bash	shell.

The	remaining	startup	files	are	all	used	for	the	same	function	—	to	provide	a	user-specific startup	file	for defining	user-specific	environment	variables.

* $HOME/.bash_profile

* $HOME/.bashrc

* $HOME/.bash_login

* $HOME/.profile


**On	 most	 distributions,	 the	 best	 place	 to	 store	 an	 individual	 user’s	 persistent	 bash	 shellvariables	is	in	the	 $HOME/.bashrc file.**

## Variable Arrays 

An	**array** is	a	variable	that	can	hold	multiple	values.	Values	can	be	referenced	either	individually	or as	a	whole	for	the	entire	array.

To	 set	 multiple	 values	 for	 an	 environment	 variable,	 just	 list	 them	 in	 parentheses,	 with values	separated	by	spaces:

```{bash}
mytest=(one	two	three	four	five)
```

To	 reference	 an	 individual	 array	 element,	 you must	 use	 a	 numerical	 index	 value,	 which	 represents	 its	 place	 in	 the	 array.	 The	 numeric value	is	enclosed	in	square	brackets:

```{bash}
echo	${mytest[2]}
```

**Environment	variable	arrays	start	with	an	index	value	of zero.**

To	 display	 an	 entire	 array	 variable,	 you	 use	 the	 asterisk	 wildcard	 character	 as	 the	 index value

```{bash}
echo	${mytest[*]}
#remove entire array
unset	mytest
```


# UNDERSTANDING Linux FILE PERMISSIONS

## Linux Security

The core of the Linux security system is the **user account**. The UID (user ID) is a numerical value, unique for each user.

### The `/etc/passwd/` file

Matchs the login name to a corresponding UID value.


The **root** user account is the administrator for the Linux system and is always assigned UID 0.

A **system account** is a special account that services running on the system use to gain access to resources on the system. Linux reserves UIDs below 500 for system accounts.


### The `/etc/shadow/` file

Provides more control over how the Linux system manages passwords. Only the root user has access to the /etc/shadow file, making it more secure than the /etc/passwd file.

The `/etc/shadow` file contains one record for each user account on the system.


### Adding a new user

The `useradd` command uses a combination of system default valuesand command line parameters to define a user account.

To see the system default values used on your Linux distribution, enter the `useradd` command with the `-D` parameter:


```{bash}
useradd -D
```


### Removing a user

By default, the `userdel` command removes only the user information from the `/etc/passwd file`. It doesn’t remove any files the account owns on the system.

If you use the `-r` parameter, `userdel` removes the user’s HOME directory, along with the user’s mail directory. However, other files owned by the deleted user account may still be on the system. This can be a problem in some environments.


### Modifying a user

#### `usermod`

Provides options for changing most of the fields in the `/etc/passwd` file.


#### `passwd` and `chpasswd`

If you just use the `passwd` command by itself, it changes your own password. Only the root user can change someone else’s password.

The `chpasswd` command reads a list of login name and password pairs (separated by a colon) from the standard input, automatically encrypts the password, and sets it for the user account.


## Using Linux Groups

Group permissions allow multiple users to share a common set of permissions for an object on the system, such as a file, directory, or device


### The `/etc/group` file

Contains information about each group used on the system.

The `groupadd` command allows you to create new groups on your system.

The `groupmod` command allows you to change the GID (using the `-g` parameter) or the group name (using the `-n` parameter) of an existing group


## Decoding File Permissions

The first character in the field defines the type of the object:

* - for files
* d for directories
* l for links
* c for character devices
* b for block devices
* n for network devices

## Changing Security Settings

### Changing permissions

The format of the `chmod` command is:

`chmod <options> <mode> <file>`


* `u` for the user
* `g` for the group
* `o` for others (everyone else)
* `a` for all of the above


* `x` assigns execute permissions only if the object is a directory or if it already had execute permissions.
* `s` sets the UID or GID on execution.t saves program text.
* `u` sets the permissions to the owner’s permissions.
* `g` sets the permissions to the group’s permissions.
* `o` sets the permissions to the other’s permissions.

The `u-x` entry removes the execute permission that the user already had.


### Changing ownership

The format of the `chown` command is:

`chown <options> <owner[.group]> <file>`


The `chgrp` command provides an easy way to change just the default group for a file or directory.



# MANAGING FILESYSTEMS


# INSTALLING SOFTWARE


## Package Management

Each of the major Linux distributions utilizes some form of a Package Management System (PMS) to control installing software applications and libraries.

Software packages are stored on servers, called **repositories**, and are accessed across the Internet via PMS utilities running on your local Linux system.

The two primary PMS base utilities commonly used in the Linux world are **dpkg** and **rpm**.


## The Debian-based Systems

These other tools are included in this PMS:

* `apt-get`
* `apt-cache`
* `aptitude`

### Managing packages with `aptitude`

If you already know the packages on your system and want to quickly display detailed information about a particular package:

`aptitude show <package_name>`



### Installing software packages with `aptitude`

Use the `aptitude` command with the search option (wildcards are implied):

`aptitude search <package_name>`

Before each package name is either a `p` or `i` . If you see an `i` , the package is currently installed on your system. If you see a `p` or `v` , it is available but not installed.

To install a software package on a system from a repository using aptitude:

`aptitude install <package_name>`


### Updating software with `aptitude`

To safely update all the software packages on a system with any new versions in the repository, use the `safe-upgrade` option:

`aptitude safe-upgrade`

The `safe-upgrade` option upgrades all the installed packages to the most recent version available in the repository, which is safer for system stabilization.


### Uninstalling software with `aptitude`

To remove a software package, but not the data and configuration files, use the `remove` option of aptitude.

To remove a software package and the related data and configuration files, use the `purge` option

`sudo aptitude purge <package_name>`


### The `aptitude` repositories

The repository locations are stored in the file `/etc/apt/sources.list`



## The Red Hat-Based Systems

These are the common tools:

* yum : Used in Red Hat and Fedora
* urpm : Used in Mandriva
* zypper :Used in openSUSE


### Listing installed packages

`yum list installed > installed_software`

To find out detailed information for a particular software package:

`yum list <package_name>`

`yum list installed <package_name>`

### Installing software with yum

`yum install <package_name>`

You can also manually download an rpm installation file and install it using yum:

`yum localinstall <package_name.rpm>`


### Updating software with yum

`yum list updates`

`yum update`


### Uninstalling software with yum

To just remove the software package and keep any configuration and data files, use the following command:

`yum remove <package_name>`

To uninstall the software and all its files, use the `erase` option:

`yum erase <package_name>`


## Installing from Source code

To unpack a software tarball, use the standard `tar` command:

`tar -zxvf <software_tarball>`

You should typically see a README or AAAREADME file. The actual instructions you need to finish the software’s installation are in this file.

`./configure`

The `make` command compiles the source code and then the linker to create the final executable files for the package.

`make`

`make install`


# WORKING WITH EDITORS


## The vim Editor

To start the vim editor, just type the `vim` command and the name of the file you want to edit:

`vim myprog.c`

The vim editor has two modes of operation:

* **Normal mode**: the editor interprets keystrokes as commands

* **Insert mode**: vim inserts every key you type at the current cursor location in the buffer

To enter insert mode, press the `i` key. To get out of insert mode and go back into normal mode, press the `Escape` key on the keyboard.

The **command line** mode provides an interactive command line where you can enter additional commands to control the actions in vim.

To get to command line mode, press the `colon :` key in normal mode:

* `q` to quit if no changes have been made to the buffer data

* `q!` to quit and discard any changes made to the buffer data

* `w <filename>` to save the file under a different filename

* `wq` to save the buffer data to the file and quit

The copy command in vim is `y` (for yank).


To enter a search string, press the forward slash (/) key. Enter the text you want to find, and press the Enter key.

The `substitute` command allows you to quickly replace (substitute) one word for another in the text:

`:s/old/new/`

* `:s/old/new/g` to replace all occurrences of old in a line

* `:n,ms/old/new/g` to replace all occurrences of old between line numbers n and m

* `:%s/old/new/g` to replace all occurrences of old in the entire file

* `:%s/old/new/gc`to replace all occurrences of old in the entire file, but prompt for each occurrence


# BASIC SCRIPT BUILDING


## Using Multiple Commands

If you want to run two commands together, you can enter them on the same prompt line, separated with a semicolon:

```{bash}
date; who
```


## Creating a Script File

When creating a shell script file, you must specify the shell you are using in the first line of the file:

`#!/bin/bash`

The pound sign followed by the exclamation point tells the shell what shell to run the script under.

After indicating the shell, commands are entered onto each line of the file, followed by a carriage return.

The next step is to give the file owner permission to execute the file, using the `chmod` command

`chmod u+x <test>`

`./test`


## Displaying Messages

Many times, however, you will want to add your own text messages to help the script user know what is happening within the script. You can do this with the `echo` command.

The `echo` command uses either double or single quotes to delineate text strings. If you use them within your string, you need to use one type of quote within the text and the other type to delineate the string

```{bash}
#!/bin/bash
# This script displays the date and who’s logged on
echo The time and date are:
date
echo “Let’s see who’s logged into the system:”
who
```


If you want to echo a text string on the same line as a command output you can use the `-n` parameter for the echo statement to do that

```{bash}
#!/bin/bash
# This script displays the date and who’s logged on
echo -n "The time and date are: "
date
echo “Let’s see who’s logged into the system:”
who
```


## Using Variables

**Variables** allow you to temporarily store information within the shell script for use with other commands in the script.

### Environment variables

You can display a complete list of active environment variables available by using the `set` command, and access these values from your shell scripts

You can tap into these environment variables from within your scripts by using the environment variable’s name preceded by a dollar sign.

```{bash}
#!/bin/bash
# display user information from the system.
echo "User info for userid: $USER"
echo UID: $UID
echo HOME: $HOME
```


We were able to place the $USER system variable within the double quotation marks in the first string, and the shell script still
figured out what we meant.

However, to display an actual dollar sign, you must precede it with a backslash character:

```{bash}
echo "The cost of the item is \$15"
```


### User variables

User variables can be any text string of up to 20 letters, digits, or an underscore character and are also **case sensitive**.


**Values are assigned to user variables using an equal sign. No spaces can appear between the variable, the equal sign, and the value**

```{bash}
#!/bin/bash
# testing variables
days=10
guest="Katie"
echo "$guest checked in $days days ago"
```

**When referencing a variable value you use the dollar sign, but when referencing the variable to assign a value to it, you do not use the dollar sign.**


### Command substitution

There are two ways to assign the output of a command to a variable:

* **The backtick character (`)**

* **The $() format**

You must either surround the entire command line command with two backtick characters:

testing=`date`

or use:

`testing=$(date)`

The variable `testing` receives the output from the `date` command, and it is used in the `echo` statement to display it:

```{bash}
#!/bin/bash
testing=$(date)
echo "The date and time are: " $testing
```


## Redirecting Input and Output


### Output redirection

The most basic type of redirection is sending output from a command to a file. The bash shell uses the greater-than symbol (`>`) for this.

You can use the double greater-than symbol (`>>`) to append data



### Input redirection

Input redirection (`<`) takes the content of a file and redirects it to a command.


**Inline input redirection** (`<<`) allows you to specify the data for input redirection on the command line instead of in a file.



## Pipes

The **pipe** is put between the commands to redirect the output from one to the other:

`command1 | command2`


## Performing Math

### The `expr` command

`expr 5 \* 2`

Many of the expr command operators have other meanings in the shell!


### Using brackets

This is a much easier way of performing mathematical equations.

When assigning a mathematical value to a variable, you can enclose the mathematical equation using a dollar sign and square brackets (`$[ operation ]`):

```{bash}
var1=$[1 + 5]
echo $var1
var2=$[$var1 * 2]
echo $var2
```


Using brackets makes shell math much easier than with the expr command. This same technique also works in shell scripts:

```{bash}
#!/bin/bash
var1=100
var2=50
var3=45
var4=$[$var1 * ($var2 - $var3)]
echo The final result is $var4
```


### A floating-point solution

The bash shell mathematical operators support only integer arithmetic. The most popular solution uses the built-in bash calculator, called `bc`.

You can access the bash calculator from the shell prompt using the `bc` command. To exit the bash calculator, you must enter `quit`.

The floating-point arithmetic is controlled by a built-in variable called `scale`. You must set this value to the desired number of decimal places you want in your answers:

```
bc -q
scale=4
3.44/5
```

### Using `bc` in scripts

You can use the **command substitution character** to run a `bc` command and assign the output to a variable. The basic format to use is this:

`variable=$(echo “options; expression” | bc)`

The first portion, `options`, allows you to set variables. If you need to set more than one variable, separate them using the **semicolon**. The `expression` parameter defines the mathematical expression to evaluate using bc .

```{bash}
#!/bin/bash
var1=$(echo "scale=4; 3.44 / 5" | bc)
echo The answer is $var1
```

You can also use variables defined in the shell script:

```{bash}
#!/bin/bash
var1=100
var2=45
var3=$(echo "scale=4; $var1 / $var2" | bc)
echo The answer for this is $var3
```


If you have more than just a couple of calculations, it gets confusing trying to list multiple expressions on the same command line. The best method is to use inline input redirection, which allows you to redirect data directly from the command line.

```
variable=$(bc << EOF
options
statements
expressions
EOF
)
```

Now you can place all the individual bash calculator elements on separate lines in the script file:

```
#!/bin/bash
var1=10.46
var2=43.67
var3=33.2
var4=71
var5=$(bc << EOF
scale = 4
a1 = ( $var1 * $var2)
b1 = ($var3 * $var4)
a1 + b1
EOF
)
echo The final answer for this mess is $var5
```

The `EOF` string indicates the start and end of the data to redirect to the bc command.


## Exiting the Script

Every command that runs in the shell uses an exit status to indicate to the shell that it’s finished processing. The exit status is an integer value between 0 and 255 that’s passed by the command to the shell when the command finishes running. You can capture this value
and use it in your scripts.


### Checking the `exit` status

The `$?` special variable holds the exit status value from the last command that executed. You must view or use the $? variable immediately after the command you want to check.

```{bash}
date
echo $?
```


**By convention, the exit status of a command that successfully completes is zero**. If a command completes with an error, then a positive integer value is placed in the exit status.


**By default, your shell script exits with the exit status of the last command in your script**.


The `exit` command allows you to specify an exit status when your script ends:

```{bash}
#!/bin/bash
# testing the exit status
var1=10
var2=30
var3=$[$var1 + $var2]
echo The answer is $var3
exit 5
```

```
echo $?
5
```

# USING STRUCTURED COMMANDS

## The `if-then` Statement

```
if <command>
then
<commands>
fi
```

In other programming languages, the object after the `if` statement is an equation that is evaluated for a `TRUE` or `FALSE` value. 

**That’s not how the bash shell if statement works.**

The bash shell if statement runs the command defined on the `if` line. If the `exit` status of the command is zero (the command completed successfully), the commands listed under the then section are executed. 

If the `exit` status of the command is anything else, the `then` commands aren’t executed, and the bash shell moves on to the next command in the script. The `fi` statement delineates the `if-then` statement’s end.

```{bash}
#!/bin/bash
# testing the if statement
if pwd
then
echo "It worked"
fi
```


You might see an alternative form of the `if-then` statement used in some scripts:

```
if <command>; then
<commands>
fi
```

You are not limited to just one command in the `then` section. You can list commands just as in the rest of the shell script.

```{bash}
#!/bin/bash
# testing multiple commands in the then section
#
testuser=nacho
#
if grep $testuser /etc/passwd
then
echo "This is my first command"
echo "This is my second command"
echo "I can even put in other commands besides echo:"
ls -a /home/$testuser/.b*
fi
```



## The `if-then-else` Statement

In the `if-then` statement, you have only one option for whether a command is successful. If the command returns a non-zero `exit` status code, the bash shell just moves on to the next command in the script.

The `if-then-else` statement provides another group of commands in the statement:

```
if <command>
then
<commands>
else
<commands>
fi
```
When the command in the `if` statement line returns a non-zero `exit` status code, the bash shell executes the commands in the `else` section.


## Nesting `if`s

Sometimes, you must check for several situations in your script code. For these situations, you can nest the `if-then` statements

Instead of having to write separate `if-then` statements, you can use an alternative version of the `else` section, called `elif`. The `elif` continues an `else` section with another `if-then` statement:

```
if <command1>
then
<commands>
elif <command2>
then
more <commands>
fi
```

```{bash}
#!/bin/bash
# Testing nested ifs - use elif
#
testuser=NoSuchUser
#
if grep $testuser /etc/passwd
then
echo “The user $testuser exists on this system.”
#
elif ls -d /home/$testuser
then
echo “The user $testuser does not exist on this system.”
echo “However, $testuser has a directory.”
#
fi
```


With an `elif` statement, any `else` statements immediately following it are for that `elif` code block. They are not part of a preceding `if-then` statement code block.


## Trying the `test` Command

The `test` command provides a way to test different conditions in an `if-then` statement. 

* If the condition listed in the `test` command evaluates to `TRUE`, the `test` command exits with a zero `exit` status code.

* If the condition is `FALSE`, the `test` command exits with a non-zero `exit` status code, which causes the `if-then` statement to exit.

`test <condition>`

```
if test <condition>
then
<commands>
fi
```

If you leave out the `condition` portion of the `test` command statement, it exits with a non-zero `exit` status code and triggers any `else` block statements:

```{bash}
#!/bin/bash
# Testing the test command
#
if test
then
echo "No expression returns a True"
else
echo "No expression returns a False"
fi
```

The bash shell provides an alternative way of testing a condition without declaring the `test` command in an `if-then` statement:

```
if [condition]
then
  <commands>
fi
```

The square brackets define the test condition. Be careful; you **must** have a space after the first bracket and a space before the last bracket, or you’ll get an error message.

The `test` command and `test` conditions can evaluate three classes of conditions:

* Numeric comparisons

* String comparisons

* File comparisons


### Using numeric comparisons

The numeric test conditions can be used to evaluate both numbers and variables.

```{bash}
#!/bin/bash
# Using numeric test evaluations
#
value1=10
value2=11
#
if [ $value1 -gt 5 ]
then
echo "The test value $value1 is greater than 5"
fi
#
if [ $value1 -eq $value2 ]
then
echo "The values are equal"
else
echo "The values are different"
fi
```


**The bottom line is that you cannot use floating-point values for test conditions.**


### Using string comparisons

```{bash}
#!/bin/bash
# mis-using string comparisons
#
val1=baseball
val2=hockey
#
if [ $val1 \> $val2 ] # we have to escape the > symbol
then
echo "$val1 is greater than $val2"
else
echo "$val1 is less than $val2"
fi
```

The `-n` and `-z` comparisons are handy when trying to evaluate whether a variable contains data:

* `-n` str1 Checks if str1 has a length greater than zero.

* `-z` str1 Checks if str1 has a length of zero.

Empty and uninitialized variables can have catastrophic effects on your shell script tests. If you’re not sure of the contents of a variable, it’s always best to test if the variable contains a value using `-n` or `-z` before using it in a numeric or string comparison.


### Using file comparisons

* The `-d` test checks to see if a specified directory exists on the system.

* The `-e` comparison allows you to check if either a file or directory object exists before you attempt to use it in your script

* To be sure that the object specified is a file and not a directory, you must use the `-f` comparison

* The `-r` comparison test wether you can read a file.

* You should use `-s` comparison to check whether a file is empty, especially if you don't want to remove a non-empty file.

* The `-w` comparison determines whether you have permission to write to a file.

* The `-x` comparison is a handy way to determine whether you have execute permission for a specific file.

* The `-O` comparison allows you to easily test whether you’re the owner of a file.

```{bash}
#!/bin/bash
# check file ownership
#
if [ -O /etc/passwd ]
then
  echo "You are the owner of the /etc/passwd file"
else
  echo "Sorry, you are not the owner of the /etc/passwd file"
fi
```

* The `-G` comparison checks the default group of a file, and it succeeds if it matches the group of the default group for the user.

* The `-nt` comparison determines whether a file is newer than another file.

* The `-ot` comparison determines whether a file is older than another file.


## Considering Compound Testing

The `if-then` statement allows you to use Boolean logic to combine tests. You can use these two Boolean operators:

* [ condition1 ] && [ condition2 ] (AND)

* [ condition1 ] || [ condition2 ] (OR)

```{bash}
#!/bin/bash
# testing compound comparisons
#
if [ -d $HOME ] && [ -w $HOME/testing ]
then
echo "The file exists and you can write to it"
else
echo "I cannot write to the file"
fi
```

Using the `AND` Boolean operator, both of the comparisons must be met.


## Advanced `if-then` Features

### Using double parentheses

The double parentheses command allows you to incorporate advanced mathematical formulas in your comparisons.

```((expression))```

The **expression** term can be any mathematical assignment or comparison expression.

```{bash}
#!/bin/bash
# using double parenthesis
#
val1=10
#
if (( $val1 ** 2 > 90 ))
then
(( val2 = $val1 ** 2 ))
echo "The square of $val1 is $val2"
fi
```


### Using double brackets

The double bracket command provides advanced features for string comparisons.

```[[expression]]```

The double bracketed **expression** uses the standard string comparison used in the test evaluations. However, it provides an additional feature that the test evaluations don’t — **pattern matching**.

```{bash}
#!/bin/bash
# using pattern matching
#
if [[ $USER == n* ]]
then
echo "Hello $USER"
else
echo "Sorry, I do not know you"
fi
```


## The `case` Command

The `case` command checks multiple values of a single variable in a list-oriented format:

```
case variable in
pattern1 | pattern2) commands1;;
pattern3) commands2;;
*) default commands;;
esac
```

The asterisk symbol is the catch-all for values that don’t match any of the listed patterns.


```{bash}
#!/bin/bash# 
#using the case command
#
case $USER in
rich | barbara | nacho)
echo "Welcome, $USER"
echo "Please enjoy your visit";;
testing)
echo "Special testing account";;
jessica)
echo "Do not forget to log off when you’re done";;
*)
echo "Sorry, you are not allowed here";;
esac
```


# MORE STRUCTURED COMMANDS

## The `for` Command

The bash shell provides the `for` command to allow you to create a loop that iterates through a series of values.

```
for var in list
do
commands
done
```
The first iteration uses the first item in the list, the second iteration the second item, and so on until all the items in the list have been used.


### Reading values in a list

The most basic use of the `for` command is to iterate through a list of values defined within the `for` command itself:

```{bash}
#!/bin/bash
# basic for command
for test in Alabama Alaska Arizona Arkansas California Colorado
do
echo The next state is $test
done
```

After the last iteration, the `$test` variable remains valid throughout the remainder of the shell script. It retains the last
iteration value (unless you change its value):

```{bash}
#!/bin/bash
# testing the for variable after the looping
for test in Alabama Alaska Arizona Arkansas California Colorado
do
echo “The next state is $test”
done
echo “The last state we visited was $test”
test=Connecticut
echo “Wait, now we’re visiting $test”
```


### Reading complex values in a list

The `for` loop assumes that each value is separated with a space.

```{bash}
#!/bin/bash
# another example of how not to use the for command
for test in "Nevada" "New Hampshire" "New Mexico" "New York" "North Carolina"
do
echo "Now going to $test"
done
```

When you use double quotation marks around a value, the shell doesn’t include the quotation marks as part of the value

### Reading a list from a variable

You accumulate a list of values stored in a variable and then need to iterate through the list.

```{bash}
#!/bin/bash
# using a variable to hold the list
list="Alabama Alaska Arizona Arkansas Colorado"
list=$list" Connecticut"
for state in $list
do
echo "Have you ever visited $state?"
done
```

The code also uses another assignment statement to add (or concatenate) an item to the existing list contained in the $list variable.

### Reading values from a command

You use **command substitution** to execute any command that produces output and then use the output of the command in the `for` command:

```{bash}
#!/bin/bash
# reading values from a file
file="states"
for state in $(cat $file)
do
echo "Visit beautiful $state"
done
```

This example uses the `cat` command in the command substitution to display the contents of the file `states`.


### Changing the field separator

The **IFS** environment variable defines a list of characters the bash shell uses as field separators. By default, the bash shell considers the following characters as field separators:

* A space

* A tab

* A newline

If the bash shell sees any of these characters in the data, it assumes that you’re starting a new data field in the list.

To solve this problem, you can temporarily change the IFS environment variable values in your shell script to restrict the characters the bash shell recognizes as field separators.

For example, if you want to change the IFS value to recognize only the newline character, you need to do this:

`IFS=$’\n’`

Adding this statement to your script tells the bash shell to ignore spaces and tabs in data values.

A **safe practice** is to save the original IFS value before changing it and then restore it when you’re finished.

```
IFS.OLD=$IFS
IFS=$‘∖n’
<use the new IFS value in code>
IFS=$IFS.OLD
```

This ensures that the IFS value is returned to the default value for future operations within the script.

If you want to specify more than one IFS character, just string them together on the assignment line:

`IFS=$’\n’:;”` # This assignment uses the newline, colon, semicolon, and double quotation mark characters as field separators.

### Reading a directory using wildcards

You can use the for command to automatically iterate through a directory of files.

To do this, you must use a wildcard character in the file or pathname. This forces the shell to use **file globbing**. 

File globbing is the process of producing filenames or pathnames that match a specified wildcard character.

```{bash}
#!/bin/bash
# iterate through all the files in a directory
for file in /home/rich/test/*
do
  if [ -d "$file" ]
  then
    echo "$file is a directory"
  elif [ -f "$file" ]
  then
    echo "$file is a file"
  fi
done
```


## THE `while` COMMAND
































































