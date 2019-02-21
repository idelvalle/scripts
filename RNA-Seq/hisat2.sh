#Rename Files

cd <directory with files>
rename -n 's/_S[0-9]{1,2}//' *.fastq.gz # To view changes
rename -v 's/_S[0-9]{1,2}//' *.fastq.gz # To visualize changes
rename -v 's/_001//' *.fastq.gz

#Run Alignment with HISAT2

for i in $(ls *.fastq.gz | rev | cut -c 13- | rev | uniq)

do 

hisat2 -p 8 --dta -x ~/Desktop/Genomes/GRCh38-HISAT2/grch38_tran/genome_tran \
-1 ${i}_R1.fastq.gz -2 ${i}_R2.fastq.gz --rna-strandness RF -t | samtools sort > ${i}.bam
samtools index ${i}.bam

done