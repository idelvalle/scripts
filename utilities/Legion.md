# Login using

ssh sejjide@legion.rc.ucl.ac.uk

#  Work in Scratch directory (200GB quote)

module avail #lists available modules

module load #loads a module

module remove #removes a module

# To copy files into Scratch:

scp <path/to/filename> sejjide@legion.rc.ucl.ac.uk:~/Scratch/<dir>

scp -r  <path/to/dirname> sejjide@login05.external.legion.ucl.ac.uk:~/Scratch/<dir>

# And viceversa using dedicated node (faster)

scp -r sejjide@login05.external.legion.ucl.ac.uk:~/Scratch/<dirname> .

# sftp

You can use sftp to log in to the remote machine, navigate through directories and use put and get to copy files from and to your local machine. 

lcd and lls are local equivalents of cd and ls so you can navigate through your local directories as you go.


sftp <remote_user_id>@<remote_hostname>
cd <remote_path>
get <remote_file>
To copy whole directory into local machine use: 
get -r <dir name> 
lcd <local_path>
put <local_file>

# Managing quota

To check your quota, run the command lquota.

Also useful is du, giving you information about your disk usage. 

For example, du -ch <dir> will give you a summary of the sizes of directory tree and subtrees, in human-readable sizes, with a total at the bottom. 

du -h --max-depth=1 will show you the totals for all top-level directories relative to where you are, plus the grand total. 

To check your quota on Legion, please run quota_check. 

# qsub

The qsub command submits your job to the batch queue.

qsub myscript.sh

# qstat

The qstat command shows the status of your jobs. By default, if you run it with no options, it shows only your jobs (and no-one else’s). This makes it easier to keep track of your jobs. If you want to get more information on a particular job, note its job ID and then use the -f and -j flags to get full output about that job

 qstat -f -j 12345

If you see that your job is in Eqw state then an error occurred before your job could begin. You can see a truncated version of the error in the output of qstat -j - this is often enough to tell what the problem is if it is a file or directory not found.

# qdel

The qdel command lets you delete a job from the queue. You need to provide qdel with a job ID like so:

qdel 12345

If you wish to learn about additional commands, please run the command "man qstat"
