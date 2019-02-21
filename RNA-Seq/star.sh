for i in $(ls *.fastq.gz | rev | cut -c 17- | rev | uniq);

do

STAR --genomeDir ~/Scratch/genome/star/43bp --sjdbGTFfile ~/Scratch/genome/genes.gtf --sjdbOverhang 42 --readFilesIn ${i}_R1_001.fa$

done;
