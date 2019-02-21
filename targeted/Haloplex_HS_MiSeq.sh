
#!/bin/bash

#Remember to make script executable- run chmod u+x!

set -e #stops the execution of a script if a command or pipeline has an error
set -u #Treat unset variables as an error when performing parameter expansion.
set -o pipefail #Causes a pipeline to return the exit status of the last command in the pipe that returned a non-zero return value.


#### FASTQC Analysis ########

mkdir -p QC_fastqc #makes dir only if it does not exists

for file in *R*.gz
do
	fastqc -t 12 ${file} -o QC_fastqc
done


###### Surecall Trimmer ######## 

# Added to .profile: export SURECALLTRIMMER='/home/nacho/src/agent/SurecallTrimmer_v4.0.1.jar'

mkdir -p vcalling

for file in $(ls *.fastq.gz | rev | cut -c 17- | rev | uniq)
do
	R1=${file}_R1*.gz
	R2=${file}_R2*.gz
	
	java -jar $SURECALLTRIMMER -fq1 ${R1} -fq2 ${R2} -halo -hs
	
done


mv *Cut*.gz vcalling
cd vcalling
rename -v 's/\.[[:alnum:]]+\_[[:alnum:]]+\_[[:alnum:]]/_cut/' *.fastq.gz


######## BWA-MEM ALIGNMENT #########

IDX='/home/nacho/Desktop/Genomes/hg19_surecall/hg19.fasta'
#IDX='/home/nacho/Desktop/Genomes/broad_hg19/hg19.genome/ucsc.hg19.fasta'

for file in $(ls *.fastq.gz | rev | cut -c 29- | rev | uniq)
do
	R1=${file}*_R1*.gz
	R2=${file}*_R2*.gz

	bwa mem -t 16 -M $IDX ${R1} ${R2} > ${file}.sam
	
done

rm -r *.fastq.gz
cd ..
mv *I2*.gz ./vcalling
cd ./vcalling

######### LocatIt ################
# Added to .profile: export SURECALLTRIMMER='/home/nacho/src/agent/LocatIt_v4.0.1.jar'

AMPLICON='/home/nacho/Desktop/Haloplex_Federica/bed/DSD_HS_Amplicons.bed'


for file in *.sam

do
	
java -jar $LOCATIT  -U -IS -OB -r -c 2500 -l ${AMPLICON} -o ${file%.*}.lc.bam ${file} ${file%.*}*I2*.fastq.gz

samtools sort ${file%.*}.lc.bam >  ${file%.*}.bam

done

mv *I2*.fastq.gz ..
rm *lc*
rm *.sam

################ VARIANT CALLING WITH GATK4 ###############

#### Convert aligned BAM to uBAM and discard problematic records using RevertSam ####

for i in *.bam
do
gatk RevertSam \
-I ${i} \
-O ${i%.*}.revert.bam \
--SANITIZE true \
--MAX_DISCARD_FRACTION 0.005 \
--ATTRIBUTE_TO_CLEAR NM \
--ATTRIBUTE_TO_CLEAR MD \
--ATTRIBUTE_TO_CLEAR AS \
--ATTRIBUTE_TO_CLEAR XS \
--ATTRIBUTE_TO_CLEAR RG \
--ATTRIBUTE_TO_CLEAR XA \
--ATTRIBUTE_TO_CLEAR SA;
#
gatk AddOrReplaceReadGroups \
-I ${i%.*}.revert.bam \
-O ${i%.*}.added.bam \
-SM ${i%.*} \
-ID ${i%.*} \
-PU ${i%.*} \
-PL ILLUMINA \
-LB DSD_HS \
-CN Barclay_House;

rm *.revert.bam

gatk MarkIlluminaAdapters \
-I ${i%.*}.added.bam \
-O ${i%.*}.markadapters.bam \
-M ${i%.*}.revertsam_metrics.txt;
#
gatk SamToFastq \
-I ${i%.*}.markadapters.bam \
-F ${i%.*}.fq \
--CLIPPING_ATTRIBUTE XT \
--CLIPPING_ACTION 2 \
--INTERLEAVE true \
--NON_PF true;


rm *.markadapters.bam


bwa mem -M -t 16 -p ${IDX} ${i%.*}.fq > ${i%.*}.bwa.sam;
rm *.fq
#
gatk MergeBamAlignment \
-R ${IDX} \
--UNMAPPED_BAM ${i%.*}.added.bam \
--ALIGNED_BAM ${i%.*}.bwa.sam \
-O ${i%.*}.bwa.bam \
--CREATE_INDEX true \
--CLIP_ADAPTERS false \
--MAX_INSERTIONS_OR_DELETIONS -1 \
--PRIMARY_ALIGNMENT_STRATEGY MostDistant \
--ATTRIBUTES_TO_RETAIN XS;

rm *.added.bam
rm *.bwa.sam

gatk BaseRecalibrator \
-R ${IDX} \
-I ${i%.*}.bwa.bam \
--known-sites '/home/nacho/Desktop/Genomes/broad_hg19/dbsnp_138.hg19.vcf' \
--known-sites '/home/nacho/Desktop/Genomes/broad_hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf' \
--known-sites '/home/nacho/Desktop/Genomes/broad_hg19/1000G_phase1.indels.hg19.vcf' \
-O ${i}.table;
#
gatk ApplyBQSR \
-R ${IDX} \
-I ${i%.*}.bwa.bam \
--bqsr-recal-file ${i}.table \
-O ${i}.recal.bam;

rm *.bwa.bam

gatk HaplotypeCaller \
-R ${IDX} \
-I ${i}.recal.bam \
-L '/home/nacho/Desktop/Haloplex_Federica/bed/DSD_HS_Covered.bed' \
--dbsnp '/home/nacho/Desktop/Genomes/broad_hg19/dbsnp_138.hg19.vcf' \
-O ${i}.vcf
done

mkdir -p bam

mv *.recal* bam

rm *.bam
find . \( -name '*bwa.bai' -or -name '*.table'  -or -name '*.properties'  -or -name '*.txt' \) -type f -delete



##### APPLY HARD FILTERS ###################################

##### Extract SNPS  AND INDELS ########

mkdir -p vcf

mv *.vcf ./vcf
mv *.idx ./vcf

cd vcf

for i in *.vcf

do

gatk SelectVariants \
-R ${IDX} \
-V ${i} \
-select-type SNP \
-O ${i%%.*}.snps.vcf;

gatk VariantFiltration \
-R ${IDX} \
-V ${i%%.*}.snps.vcf \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filter-name "my_snp_filter" \
-O ${i%%.*}.filtered.snps.vcf;

gatk SelectVariants \
-R ${IDX} \
-V ${i} \
-select-type INDEL \
-O ${i%%.*}.indels.vcf;

gatk VariantFiltration \
-R ${IDX} \
-V ${i%%.*}.indels.vcf \
--filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
--filter-name "my_indel_filter" \
-O ${i%%.*}.filtered.indels.vcf

done


######### ANNOVAR #############

mkdir -p annovar

for i in *filtered.snps.vcf

do 

table_annovar.pl ${i} '/home/nacho/src/annovar/humandb/'  -buildver hg19 -out annovar/${i%%.*}.snps -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation g,r,f,f,f -nastring . -vcfinput

done

for i in *filtered.indels.vcf

do 

table_annovar.pl ${i} '/home/nacho/src/annovar/humandb/'  -buildver hg19 -out annovar/${i%%.*}.indels -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation g,r,f,f,f -nastring . -vcfinput

done

#################################


#for i in *filtered.snps.vcf

#do

#vt decompose -s ${i} | vt normalize -r ${IDX} - > ${i%.*}.vt.vcf

#done

#for i in *filtered.indels.vcf

#do

#vt decompose -s ${i} | vt normalize -r ${IDX} - > ${i%.*}.vt.vcf

#done

#for i in *.vt.vcf

#do

#java -Xmx16G -jar '/home/nacho/src/snpEff/snpEff.jar' -c '/home/nacho/src/snpEff/snpEff.config' hg19 ${i} -classic -formatEff > ${i%.*}.snpEff.vcf

#done




















