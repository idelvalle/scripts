# GATK4 BEST PRACTICES



## 1. GENERATE UNMAPPED BAM FROM FASTQ

Adapted from https://gatkforums.broadinstitute.org/gatk/discussion/6483

**REQUIRED**: *.fasta, .fai and .dict* files. Obtained from GATK4 bundle and stored in Genomes/broad_hg38

https://software.broadinstitute.org/gatk/download/bundle

If your reads are mapped, or in BCL or FASTQ format, then generate an unmapped BAM according to the following instructions.

* To convert FASTQ or revert aligned BAM files, follow following directions (**Sections A and B**), adapted from  [Tutorial#6484](http://gatkforums.broadinstitute.org/discussion/6484/#latest#top). The resulting uBAM needs to have its adapter sequences marked as outlined in the next step (step 2).
* To convert an Illumina Base Call files (BCL) use [IlluminaBasecallsToSam](http://broadinstitute.github.io/picard/command-line-overview.html#IlluminaBasecallsToSam). The tool marks adapter sequences at the same time. The resulting uBAMXT has adapter sequences marked with the XT tag so you can skip step 2 of this workflow and go directly to step 3.

### (A) Convert FASTQ to uBAM and add read group information using FastqToSam

Picard's [FastqToSam](https://broadinstitute.github.io/picard/command-line-overview.html#FastqToSam) transforms a FASTQ file to an unmapped BAM, requires two read group fields and makes optional specification of other read group fields. 

**1)** Check the header of fastq files with `head read_1.fastq` 

Each entry in an Illumina FASTQ file consists of four lines: (see  [IlluminaFormat](http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm) ):

* Sequence Identifier **@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<sample number>**
* Sequence
* Quality score identifier line (consisting only of one + )
* Quality score

**2)** Add read group field required by GATK (see [ReadGroups](https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups#latest) for how-to and discussion :

`ID`: <flowcell ID>.<sample number>.<lane>

`PU`: Same as `ID`

`SM`: Sample Name

`PL`: ILLUMINA

`LB`: Library Name

We use the following command to add fields required for GATK Best Practices Workflows:

```bash
gatk FastqToSam --help \
gatk FastqToSam \
-F1 Rokia-MS_S10_L001_R1_001.fastq.gz \
-F2 Rokia-MS_S10_L001_R2_001.fastq.gz \
-O Rokia-MS.bam \
-SM Rokia-MS \
-RG 000000000-AL5PT.10.1 \
-PU 000000000-AL5PT.10.1 \
-PL ILLUMINA \
-LB tulay_advanced_0515_v2 \
-CN Barclay_House
```

We can confirm the addition of the fields using:

`samtools view -H Rokia-MS.bam | grep '@RG'`

### **(B) Convert aligned BAM to uBAM and discard problematic records using RevertSam**

First we list all tags within BAM file:

```bash
samtools view Rokia-MS.raw.bam \
| cut -f 12- | tr '\t' '\n' | cut -d ':' -f 1 | awk '{ if(!x[$1]++) { print }}' 
```

We use Picard's [RevertSam](https://broadinstitute.github.io/picard/command-line-overview.html#RevertSam) to remove alignment information and generate an unmapped BAM (uBAM). 

```bash
gatk RevertSam --help \
gatk RevertSam \
-I Rokia-MS.raw.bam \
-O Rokia-MS.raw.revertsam.bam \
--SANITIZE true \
--MAX_DISCARD_FRACTION 0.005 \
--ATTRIBUTE_TO_CLEAR NM \
--ATTRIBUTE_TO_CLEAR MD \
--ATTRIBUTE_TO_CLEAR AS \
--ATTRIBUTE_TO_CLEAR XS \
--ATTRIBUTE_TO_CLEAR RG \
--ATTRIBUTE_TO_CLEAR XA \
--ATTRIBUTE_TO_CLEAR SA \
```

Next we add groups according to the BAM file  [Fix BAM](https://gatkforums.broadinstitute.org/gatk/discussion/2909) using Picard's [AddOrReplaceReadGroups](http://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups) to appropriately label read group (`@RG`) fields, coordinate sort and index a BAM file.

```bash
gatk AddOrReplaceReadGroups --help \
gatk AddOrReplaceReadGroups \
-I Rokia-MS.raw.revertsam.bam \
-O Rokia-MS.raw.revertsam.added.bam \
-SM Rokia-MS \ #sample name
-ID 000000000-AL5PT.10.1 \
-PU 000000000-AL5PT.10.1 \ #platform unit
-PL ILLUMINA \ #platform used
-LB tulay_advanced_0515_v2 \ #library used
-CN Barclay_House #sequencing center
```

We confirm the addition by using:

```bash
samtools view -H Rokia-MS.raw.revertsam.added.bam | grep '@RG'
```



## 2. Mark adapter sequences using MarkIlluminaAdapters

[MarkIlluminaAdapters](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.1.2/picard_illumina_MarkIlluminaAdapters.php) adds the XT tag to a read record to mark the 5' start position of the specified adapter sequence and produces a metrics file. Tools such as `SamToFastq` use the XT tag in various ways to effectively remove adapter sequence contribution to read alignment and alignment scoring metrics. Depending on your library preparation, insert size distribution and read length, expect varying amounts of such marked reads.

The tool reads a SAM or BAM file and rewrites it with new adapter-trimming tags. Clears any existing adapter-trimming tags (XT:i:) in the optional tag region of a SAM file. The SAM/BAM file must be sorted by query name. Outputs a metrics file histogram showing counts of bases_clipped per read.

### A) Use with FASTQ Converted to uBAM

```bash
gatk MarkIlluminaAdapters \
-I Rokia-MS.bam \
-O Rokia-MS.markilluminaadapters.bam \
-M Rokia-MS.markilluminaadapters_metrics.txt
```

- By default, the tool uses Illumina adapter sequences.
- Adjust the default standard Illumina adapter sequences to any adapter sequence using the `FIVE_PRIME_ADAPTER` and `THREE_PRIME_ADAPTER` parameters. To clear and add new adapter sequences first set `ADAPTERS` to 'null' then specify each sequence with the parameter.

We can plot the metrics data (is in [GATKReport file format](https://www.broadinstitute.org/gatk/guide/article?id=1244) using RStudio), where marked bases vary in size up to the full length of reads.

### B) Use with BAM Converted to uBAM

```bash
gatk MarkIlluminaAdapters \
-I Rokia-MS.raw.revertsam.added.bam \
-O Rokia-MS.raw.revertsam.added.markilluminaadapters.bam \
-M Rokia-MS.raw.revertsam_metrics.txt
```



## 3. Align reads with BWA-MEM and merge with uBAM using MergeBamAlignment

This step actually pipes three processes, performed by three different tools. For larger data, however, using [Unix pipelines](https://en.wikipedia.org/wiki/Pipeline_(Unix)) can add up to significant savings in processing time and storage.

The three tools we pipe are SamToFastq, BWA-MEM and MergeBamAlignment. By piping these we bypass storage of larger intermediate FASTQ and SAM files. 

### 3.A Convert BAM to FASTQ and discount adapter sequences using SamToFastq 

Picard's [SamToFastq](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.2.0/picard_sam_SamToFastq.php) takes read identifiers, read sequences, and base quality scores to write a Sanger FASTQ format file. We use additional options to effectively remove previously marked adapter sequences, in this example marked with an XT tag. By specifying `CLIPPING_ATTRIBUTE`=XT and `CLIPPING_ACTION`=2, SamToFastq changes the quality scores of bases marked by XT to two--a rather low score in the Phred scale. This effectively removes the adapter portion of sequences from contributing to downstream read alignment and alignment scoring metrics.

```bash
gatk SamToFastq \
-I Rokia-MS.markilluminaadapters.bam \
-F Rokia-MS.markilluminaadapters.fq \
--CLIPPING_ATTRIBUTE XT \
--CLIPPING_ACTION 2 \
--INTERLEAVE true \
--NON_PF true
```

```bash
gatk SamToFastq \
-I Rokia-MS.raw.revertsam.added.markilluminaadapters.bam \
-F Rokia-MS.raw.revertsam.added.markilluminaadapters.fq \
--CLIPPING_ATTRIBUTE XT \
--CLIPPING_ACTION 2 \
--INTERLEAVE true \
--NON_PF true
```

This step produces a FASTQ file in which all extant meta data, i.e. read group information, alignment information, flags and tags are purged. What remains are the read query names prefaced with the `@`symbol, read sequences and read base quality scores.

- For our paired reads example file we set SamToFastq's `INTERLEAVE` to true. During the conversion to FASTQ format, the query name of the reads in a pair are marked with /1 or /2 and paired reads are retained in the same FASTQ file. [BWA aligner](http://bio-bwa.sourceforge.net/bwa.shtml) accepts interleaved FASTQ files given the `-p` option.
- We change the `NON_PF`, aka `INCLUDE_NON_PF_READS`, option from default to true. SamToFastq will then retain reads marked by what [some consider an archaic 0x200 flag bit](https://github.com/samtools/hts-specs/issues/85) that denotes reads that do not pass quality controls, aka reads failing platform or vendor quality checks. Our tutorial data do not contain such reads and we call out this option for illustration only.
- Other CLIPPING_ACTION options include (1) X to hard-clip, (2) N to change bases to Ns or (3) another number to change the base qualities of those positions to the given value.

### 3B. Align reads and flag secondary hits using BWA-MEM

 Alignment is the most compute intensive and will take the longest time. GATK's variant discovery workflow recommends Burrows-Wheeler Aligner's maximal exact matches (BWA-MEM) algorithm.

BWA-MEM is suitable for aligning high-quality long reads ranging from 70 bp to 1 Mbp against a large reference genome such as the human genome.

 BWA alignment requires an indexed reference genome file. Indexing is specific to algorithms. To index the human genome for BWA, we apply BWA's `index` function on the reference genome file, e.g. `Homo_sapiens_assembly38.fasta`, obtained from the FTP Server GTAK Resource Bundle: [GATK_Bundle](https://software.broadinstitute.org/gatk/download/bundle). This produces five index files with the extensions `amb`, `ann`, `bwt`, `pac` and `sa`.

```
bwa index -a bwtsw ucsc.hg19.fasta.gz
```

To align data against reference genome we apply:

```bash
bwa mem -M -t 16 -p /path/ucsc.hg19.fasta.gz \
Rokia-MS.markilluminaadapters.fq > Rokia-MS_bwa.sam
```

```bash
bwa mem -M -t 16 -p /path/Homo_sapiens_assembly38.fasta \
Rokia-MS.raw.revertsam.added.markilluminaadapters.fq > Rokia-MS.raw_bwa.sam
```

We invoke three options in the command.

- `-M` to flag shorter split hits as secondary. 

  This is optional for Picard compatibility as MarkDuplicates can directly process BWA's alignment, whether or not the alignment marks secondary hits. However, if we want MergeBamAlignment to reassign proper pair alignments, to generate data comparable to that produced by the Broad Genomics Platform, then we must mark secondary alignments.

- `-p` to indicate the given file contains interleaved paired reads.

- `-t` followed by a number for the number of processor threads to use concurrently. Check your server or system's total number of threads with the following command:

  `getconf _NPROCESSORS_ONLN`

### 3C. Restore altered data and apply & adjust meta information using MergeBamAlignment 

Broadly, the tool merges defined information from the unmapped BAM (uBAM, step 1) with that of the aligned BAM (step 3) to conserve read data, e.g. original read information and base quality scores. 

 [MergeBamAkignment](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/picard_sam_MergeBamAlignment.php) also generates additional meta information based on the information generated by the aligner, which may alter aligner-generated designations, e.g. mate information and secondary alignment flags. The tool then makes adjustments so that all meta information is congruent, e.g. read and mate strand information based on proper mate designations. We ascribe the resulting BAM as *clean*.

As the tool name implies, MergeBamAlignment applies read group information from the uBAM and retains the program group information from the aligned BAM.

```bash
gatk MergeBamAlignment \
-R/path/Homo_sapiens_assembly38.fasta \
--UNMAPPED_BAM Rokia-MS.bam \
--ALIGNED_BAM Rokia-MS_bwa.sam \
-O Rokia-MS_mergebamalignment.bam \
--CREATE_INDEX true \
--CLIP_ADAPTERS false \
--MAX_INSERTIONS_OR_DELETIONS -1 \
--PRIMARY_ALIGNMENT_STRATEGY MostDistant \
--ATTRIBUTES_TO_RETAIN XS
```

```bash
gatk MergeBamAlignment \
-R/path/Homo_sapiens_assembly38.fasta \
--UNMAPPED_BAM Rokia-MS.raw.revertsam.added.bam \
--ALIGNED_BAM Rokia-MS.raw_bwa.sam \
-O Rokia-MS.raw_mergebamalignment.bam \
--CREATE_INDEX true \
--CLIP_ADAPTERS false \
--MAX_INSERTIONS_OR_DELETIONS -1 \
--PRIMARY_ALIGNMENT_STRATEGY MostDistant \
--ATTRIBUTES_TO_RETAIN XS
```

This generates a coordinate-sorted and *clean* BAM, and corresponding `.bai` index. These are ready for further analyses starting with MarkDuplicates.

 **Comments on select parameters**

- Setting `PRIMARY_ALIGNMENT_STRATEGY`to MostDistant marks primary alignments based on the alignment *pair* with the largest insert size. This strategy is based on the premise that of chimeric sections of a read aligning to consecutive regions, the alignment giving the largest insert size with the mate gives the most information.
- Setting `MAX_INSERTIONS_OR_DELETIONS` to -1 retains reads irregardless of the number of insertions and deletions. The default is 1.
- `ATTRIBUTES_TO_RETAIN` is specified to carryover the XS tag from the alignment, which reports BWA-MEM's suboptimal alignment scores. 
- Setting `CLIP_ADAPTERS` to false leaves reads unclipped.
- By default the merged file is coordinate sorted. We set `CREATE_INDEX` to true to additionally create the `bai` index.

### 3D. PIPE TO COMBINE EVERY STEP

In the piped command, the commands for the three processes are given together, separated by a vertical bar `|`. We also replace each intermediate output and input file name with a symbolic path to the system's output and input devices, here `/dev/stdout` and `/dev/stdin`, respectively. We need only provide the first input file and name the last output file.

We should [ask UNIX to stop the piped command](https://sipb.mit.edu/doc/safe-shell/) if any step of the pipe should error and also return to us the error messages:

```bash
set -o pipefail
gatk SamToFastq \
-I Rokia-MS.markilluminaadapters.bam \
-F /dev/stdout \
--CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 --INTERLEAVE true --NON_PF true \
--TMP_DIR /tmp/ | \ 
bwa mem -M -t 16 \
-p '/home/nacho/Desktop/Genomes/broad_hg38/Homo_sapiens_assembly38.fasta' /dev/stdin | \  gatk MergeBamAlignment \
--ALIGNED_BAM /dev/stdin \
--UNMAPPED_BAM Rokia-MS.bam \ 
--OUTPUT Rokia-MS_mergebamalignment.bam \
-R '/home/nacho/Desktop/Genomes/broad_hg38/Homo_sapiens_assembly38.fasta' \
--CREATE_INDEX true --ADD_MATE_CIGAR true \
--CLIP_ADAPTERS false --CLIP_OVERLAPPING_READS true \
--INCLUDE_SECONDARY_ALIGNMENTS true --MAX_INSERTIONS_OR_DELETIONS -1 \
--PRIMARY_ALIGNMENT_STRATEGY MostDistant --ATTRIBUTES_TO_RETAIN XS \
--TMP_DIR /tmp/
```



We have produced a ***clean*** BAM that is coordinate-sorted and indexed, in an efficient manner that minimizes processing time and storage needs. The file is ready for further steps, such as marking duplicates.

For **multiplexed samples**, first perform the workflow steps on a file representing one sample and one lane. Then mark duplicates. Later, after some steps in the GATK's variant discovery workflow, and after aggregating files from the same sample from across lanes into a single file, mark duplicates again. These two marking steps ensure you find both optical and PCR duplicates.



## 4. MARK DUPLICATES TO MITIGATE DUPLICATION ARTIFACTS

For Amplicon-based such as Haloplex samples-we are NOT marking duplicates.

Otherwise (e.g SureSelect-capture hybridization) follow  the following [tutorial](https://gatkforums.broadinstitute.org/gatk/discussion/2799#latest).



## 5. BASE RECALIBRATION CORRECTS FOR INSTRUMENT ERRORS



```bash
gatk BaseRecalibrator \
-R '/home/nacho/Desktop/Genomes/broad_hg19/hg19.genome/ucsc.hg19.fasta' \
-I Rokia-MS_mergebamalignment.bam \
--known-sites '/home/nacho/Desktop/Genomes/broad_hg19/dbsnp_138.hg19.vcf' \
--known-sites '/home/nacho/Desktop/Genomes/broad_hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf' \
--known-sites '/home/nacho/Desktop/Genomes/broad_hg19/1000G_phase1.indels.hg19.vcf' \
-O Rokia-MS_mergebamalignment.table
```



## 6. APPLY BASE QUALITY SCORE RECALIBRATION

This step generates analysis-ready BAM Files

```bash
gatk ApplyBQSR \
-R '/home/nacho/Desktop/Genomes/broad_hg19/hg19.genome/ucsc.hg19.fasta' \
-I Rokia-MS_mergebamalignment.bam \
--bqsr-recal-file Rokia-MS_mergebamalignment.table \
-O Rokia-MS_mergebamalignment.recal.bam
```



## 7. CALL VARIANTS WITH HaplotypeCaller

```bash
gatk HaplotypeCaller \
-R '/home/nacho/Desktop/Genomes/broad_hg19/hg19.genome/ucsc.hg19.fasta' \
-I Rokia-MS_mergebamalignment.recal.bam \
-L Covered.bed \
--genotyping-mode DISCOVERY \
-O raw_variants.vcf
```



## 8. APPLY HARD FILTERS TO CALL SET

The preferred method for filtering variants after the calling step is to use [VQSR](https://software.broadinstitute.org/gatk/documentation/article?id=11084), a.k.a. recalibration. However, it requires well-curated training/truth resources, which are typically not available for organisms other than humans, and it also requires a large amount of variant sites to operate properly, so it is **not suitable for some small-scale experiments such as targeted gene panels** or exome studies with fewer than 30 exomes. We will use the following [recommendations](https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set):

### 8.A Extract SNPs from call set

By using the tool [SelectVariants](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php) we create a VCF file called `raw_snps.vcf`, containing just the SNPs from the original file of raw variants.

```bash
gatk SelectVariants \
-R '/home/nacho/Desktop/Genomes/broad_hg19/hg19.genome/ucsc.hg19.fasta' \
-V raw_variants.vcf \
-select-type SNP \
-O raw_snps.vcf
```

### 8.B Determine Parameters for Filtering SNPs

SNPs matching any of these conditions will be considered bad and filtered out, *i.e.* marked `FILTER` in the output VCF file. The program will specify which parameter was chiefly responsible for the exclusion of the SNP using the culprit annotation. SNPs that do not match any of these conditions will be considered good and marked `PASS` in the output VCF file.

 Here are some recommended arguments to use with VariantFiltration when ALL other options are unavailable to you. Be sure to read the documentation explaining [how to understand and improve upon these recommendations](https://www.broadinstitute.org/gatk/guide/article?id=6925).

Note that these JEXL expressions will tag as filtered any sites where the annotation value **matches** the expression. So if you use the expression `QD < 2.0`, any site with a QD lower than 2 will be tagged as failing that filter.

 #### For SNPs:

- `QD < 2.0` This is the variant confidence (from the `QUAL` field) divided by the unfiltered depth of non-reference samples.

- `MQ < 40.0` This is the Root Mean Square of the mapping quality of the reads across all samples.

   

- `FS > 60.0` Phred-scaled p-value using Fisher’s Exact Test to detect strand bias (the variation being seen on only the forward or only the reverse strand) in the reads. More bias is indicative of false positive calls.

- `SOR > 3.0` The StrandOddsRatio annotation is one of several methods that aims to evaluate whether there is strand bias in the data. Higher values indicate more strand bias.

- `MQRankSum < -12.5` This is the u-based z-approximation from the Mann-Whitney Rank Sum Test for mapping qualities (reads with ref bases vs. those with the alternate allele). Note that the mapping quality rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles, *i.e.* this will only be applied to heterozygous calls.

- `ReadPosRankSum < -8.0` This is the u-based z-approximation from the Mann-Whitney Rank Sum Test for the distance from the end of the read for reads with the alternate allele. If the alternate allele is only seen near the ends of reads, this is indicative of error. Note that the read position rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles, *i.e.* this will only be applied to heterozygous calls.

   

If your callset was generated with UnifiedGenotyper for legacy reasons, you can add `HaplotypeScore > 13.0`.

#### For indels:

- `QD < 2.0`
- `ReadPosRankSum < -20.0`
- `InbreedingCoeff < -0.8`
- `FS > 200.0`
- `SOR > 10.0`

####  IMPORTANT caveats

- The `InbreedingCoeff` statistic is a population-level calculation that is only available with **10** or more samples. If you have fewer samples you will need to omit that particular filter statement.
- For shallow-coverage (**<10x**), it is virtually impossible to use manual filtering to reliably separate true positives from false positives. You **SHOULD** use the protocol involving variant quality score recalibration. If you can't do that, maybe you need to take a long hard look at your experimental design.
- The maximum DP (depth) filter only applies to whole genome data, where the probability of a site having exactly N reads given an average coverage of M is a well-behaved function. First principles suggest this should be a binomial sampling but in practice it is more a Gaussian distribution. Regardless, the DP threshold should be set a 5 or 6 sigma from the mean coverage across all samples, so that the DP > X threshold eliminates sites with excessive coverage caused by alignment artifacts. Note that **for exomes, a straight DP filter shouldn't be used** because the relationship between misalignments and depth isn't clear for capture data.



### 8.C Apply the Filter to the SNP call set

We will use the [VariantFiltration tool](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php):

```bash
gatk VariantFiltration \
-R '/home/nacho/Desktop/Genomes/broad_hg19/hg19.genome/ucsc.hg19.fasta' \
-V raw_snps.vcf \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filter-name "my_snp_filter" \
-O filtered_snps.vcf
```

This creates a VCF file called `filtered_snps.vcf`, containing all the original SNPs from the `raw_snps.vcf`file, but now the SNPs are annotated with either `PASS` or `FILTER` depending on whether or not they passed the filters.

For SNPs that failed the filter, the variant annotation also includes the name of the filter. That way, if you apply several different filters (simultaneously or sequentially), you can keep track of which filter(s) each SNP failed, and later you can retrieve specific subsets of your calls using the SelectVariants tool. 



### 8.D Extract the Indels from the call set

We will use again the tool [SelectVariants](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php) to generate a VCF file containing just the Indels from the original file of raw variants.

```bash
gatk SelectVariants \
-R '/home/nacho/Desktop/Genomes/broad_hg19/hg19.genome/ucsc.hg19.fasta' \
-V raw_variants.vcf \
-select-type INDEL \
-O raw_indels.vcf
```

### 8.E Determine Parameters for Filtering Indels

Indels matching any of these conditions will be considered bad and filtered out, *i.e.* marked `FILTER` in the output VCF file. The program will specify which parameter was chiefly responsible for the exclusion of the indel using the culprit annotation. Indels that do not match any of these conditions will be considered good and marked `PASS` in the output VCF file.

- QualByDepth (QD) 2.0

This is the variant confidence (from the `QUAL` field) divided by the unfiltered depth of non-reference samples.

- FisherStrand (FS) 200.0

Phred-scaled p-value using Fisher’s Exact Test to detect strand bias (the variation being seen on only the forward or only the reverse strand) in the reads. More bias is indicative of false positive calls.

- ReadPosRankSumTest (ReadPosRankSum) 20.0

This is the u-based z-approximation from the Mann-Whitney Rank Sum Test for the distance from the end of the read for reads with the alternate allele. If the alternate allele is only seen near the ends of reads, this is indicative of error. Note that the read position rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles, *i.e.* this will only be applied to heterozygous calls.

- StrandOddsRatio (SOR) 10.0

The StrandOddsRatio annotation is one of several methods that aims to evaluate whether there is strand bias in the data. Higher values indicate more strand bias.



### 8.F Apply the Filter to the Indel call set

We use again the [VariantFiltration tool](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php):

```bash
gatk VariantFiltration \
-R '/home/nacho/Desktop/Genomes/broad_hg19/hg19.genome/ucsc.hg19.fasta' \
-V raw_indels.vcf \
--filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
--filter-name "my_indel_filter" \
-O filtered_indels.vcf
```



 ## 9. PREPROCESSING AND LOADING A VCF FILE INTO GEMINI



### 9.A Split, Left-align and Trim Variants

Variants with multiple alternate alleles will not be handled correctly by gemini (or by the tools used to annotate the variants). As projects get more samples it is likely that a non-negligible percentage of site will have multiple alternate alleles.

In addition, variants that are not left-aligned and trimmed can be incorrectly (or not) annotated.

To reduce the number of false negatives, **we strongly recommend that gemini users split, left-align, and trim their variants**. The tool recommended by [Gemini](https://gemini.readthedocs.io/en/latest/content/preprocessing.html) is **vt**:

```bash
vt decompose -s $VCF | vt normalize -r $REFERENCE - > $NEW_VCF
```

GEMINI uses the allele depths from the **AD** tag. In order for **vt** to decompose correctly, users will have to change the #INFO field for AD in the header from Number=. to Number=R.

Then the $NEW_VCF can be annotated with snpEff or VEP.



### 9.B Annotate with snpEff or VEP

GEMINI supports gene/transcript level annotations from snpEff and VEP and hence we suggest that you first annotate your VCF with either of these tools, prior to loading it into GEMINI. 

The related database columns would be populated, which would otherwise be set to None if an unannotated VCF file is loaded into GEMINI. 

```bash
java -Xmx16G -jar path/to/snpEff/snpEff.jar -c path/to/snpEff/snpEff.config hg19 new.VCF \ -classic -formatEff > NEW_snpeff.vcf
```



```bash
vep -i '/home/nacho/Desktop/For Nacho/annotation/snpeff/filtered_snps.norm.vcf' \
-o filtered_snps.vep.vcf --vcf --offline --cache --sift b --polyphen b --symbol \
--numbers --biotype --total_length --canonical --ccds --fields \ Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,CANONICAL,CCDS,RadialSVM_score,RadialSVM_pred,LR_score,LR_pred,CADD_raw,CADD_phred,Reliability_index,LoF,LoF_filter,LoF_flags \
```





 

 









