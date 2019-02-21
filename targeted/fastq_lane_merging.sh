
#!/bin/bash

for i in $(ls *.fastq.gz | rev | cut -c 22- | rev | uniq)
	
	do echo "Merging R1"

cat "$i"_L00*_R1_001.fastq.gz > "$i"_ME_L001_R1_001.fastq.gz

       echo "Merging R2"

cat "$i"_L00*_R2_001.fastq.gz > "$i"_ME_L001_R2_001.fastq.gz

	echo "Merging R3"

cat "$i"_L00*_R3_001.fastq.gz > "$i"_ME_L001_R3_001.fastq.gz

mv *ME*.gz ../../Analysis

done


