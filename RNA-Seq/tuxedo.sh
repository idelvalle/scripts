#!/bin/bash
set -e
set -x

cd
cp -r ~/Desktop/RNA-seq ~/Desktop/RNA-seq_TEST
cd ~/Desktop/RNA-seq_TEST
head -n 20 data/2cells_1.fastq 

#TOPHAT
tophat --solexa-quals -g 2 --library-type fr-unstranded -j annotation/Danio_rerio.Zv9.66.spliceSites -o tophat/ZV9_XXX genome/ZV9 data/6h_1.fastq data/6h_2.fastq 
samtools sort tophat/ZV9_XXX/accepted_hits.bam tophat/ZV9_XXX/accepted_hits.sorted
samtools index tophat/ZV9_XXX/accepted_hits.sorted.bam 

#CUFFLINKS
cufflinks -o cufflinks/ZV9_XXX_gff -G annotation/Danio_rerio.Zv9.66.gtf -b genome/Danio_rerio.Zv9.66.dna.fa -u --library-type fr-unstranded tophat/ZV9_XXX/accepted_hits.bam 
wc -l cufflinks/ZV9_XXX_gff/transcripts.gtf

#CUFFDIFF
cuffdiff -o cuffdiff/ -L ZV9_2cells,ZV9_XXX -T -b genome/Danio_rerio.Zv9.66.dna.fa -u --library-type fr-unstranded annotation/Danio_rerio.Zv9.66.gtf tophat/ZV9_2cells/accepted_hits.bam tophat/ZV9_XXX/accepted_hits.bam 
head -n 20 cuffdiff/gene_exp.diff 
sort -t$'\t' -g -k 13 cuffdiff/gene_exp.diff > cuffdiff/gene_exp_qval.sorted.diff
head -n 20 cuffdiff/gene_exp_qval.sorted.diff

egrep 'yes$' cuffdiff/gene_exp_qval.sorted.diff | wc -l

#rm -rf tophat/ZV9_XXX/ cufflinks/ZV9_XXX_gff/ cuffdiff/*

echo "Test successful. RUN:"
echo "rm -rf ~/Desktop/RNA-seq_TEST"
