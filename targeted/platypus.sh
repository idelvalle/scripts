
#!/bin/bash

#Remember to make script executable- run chmod u+x!

set -e #stops the execution of a script if a command or pipeline has an error
set -u #Treat unset variables as an error when performing parameter expansion.
set -o pipefail #Causes a pipeline to return the exit status of the last command in the pipe that returned a non-zero return value.

###### Platypus Variant Calling ##########

#NONACUS=~/src/NonacusTools_1_Read/consensus.py

#BATCH=~/Desktop/Federica_Nonacus/fastq/batch_list.txt

PLATYPUS=~/src/Platypus_0.8.1/Platypus.py

REF=~/src/NonacusTools_1_Read/genomeRef/Homo_sapiens_assembly38.fasta

REGIONS=~/Desktop/Federica_Nonacus/Design_Files/Buonocore_Design_v1.0_target_merged.bed

# pyhon ${NONACUS} -i ${BATCH} -t32

for i in F*.bam;

do 

$PLATYPUS callVariants --bamFiles=${i} --refFile=${REF} --filterDuplicates=0 --regions=${REGIONS} --output=${i%.*}.vcf 

done

