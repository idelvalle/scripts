
#!/bin/bash

#Remember to make script executable- run chmod u+x!

set -e #stops the execution of a script if a command or pipeline has an error
set -u #Treat unset variables as an error when performing parameter expansion.
set -o pipefail #Causes a pipeline to return the exit status of the last command in the pipe that returned a non-zero return value.

######### VCALL FILES #########

PT='/home/nacho/Desktop/Federica_Nonacus/fastq/platypus.sh' 

for f in */; do cd "$f" && bash ${PT} ; cd .. ; done


