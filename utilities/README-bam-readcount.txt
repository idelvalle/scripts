
#Necessary .bai file for each .bam file

#Analysis done with bam-readcount: https://github.com/genome/bam-readcount

 for i in *.bam; do bam-readcount -f '/home/nacho/Desktop/Genomes/hg19_surecall/hg19.fasta' ${i} chr11:2905900-2905906 > ${i}_exon1.txt; done

for i in *.bam; do bam-readcount -f '/home/nacho/Desktop/Genomes/hg19_surecall/hg19.fasta' ${i} chr11:2905342-2905364 > ${i}_exon2.txt; done
