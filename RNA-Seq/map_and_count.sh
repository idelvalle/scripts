
### MAP AND COUNT ###

for i in *.bam
do
    samtools sort ${i} -o ${i}.sorted.bam;
    samtools index ${i}.sorted.bam;
    htseq-count -f bam -s no -a 10 ${i}.sorted.bam ~/Desktop/Genomes/genes.gtf > ${i}_htseq_counts.txt;
done
