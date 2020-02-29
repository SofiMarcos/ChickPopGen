## Bioinformatic steps for variant calling from FASTQ

### Softwares: 
- BWA
- SAMtools
- BCFtools
- Picard Tools
- BEDtools 
- GATK
- R

### Files:
- Data (raw_reads.fq)
- Reference genome (reference.fa)
- Database of known variants (dbssnp.vcf)

### Steps:

### 0- Preparing the reference genome
  #### 0.1- Generate the fasta file index: 
  ```
  samtools faid reference.fa
  ```
  This creates a collection of files used by BWA to perform the alignement. Among these files, there will be a file called  *reference.fa.fai* 
  This files contains one record per line for each of the contigs in the reference file
  Each record is composed of the contig name,  size, location, basePerLine and bytesPerLine. 

  #### 0.2- Generate the sequence dictionary: 
  ```
  picard CreateSequenceDictionary REFERENCE=reference.fa OUTPUT=reference.dict
  ```
  This creates a file called *reference.dic* formatted like a SAM header
  Describes the contents of your reference FASTA file

  #### 0.3- Preprare a read group information: 
  Here is where you enter the meta-data about your sample. This only will be visibly to analysis tools. 
  It will be used for BWA aligner
  ```
  @RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1
  ```

  #### 0.4- Create a text file with chromosome names for multi-sample variant calling by chromosome



### 1- Mapping to the reference genome
  #### 1.1- Generate a SAM file containing aligned reads:
  ```
  bwa mem -R '<read group info>' -p reference.fa raw_reads.fq > aligned_reads.sam
  ```
  Flags that might be interesting: 
  -t number of threads 
  -p assume that the first input query file is interleaved paired-end FASTQ
  -M mark BWA-MEM may produce multiple primary alignments for different part of a query sequence (crucial for long sequences)
  but some tools does not work with split alignments. Use option -M to flag shorter split hits as secondary (for picard compatibility)
  -R read group 

  This creates a file called *aligned_reads.sam* containing the aligned reads from all input files, combined, annotated and aligned to the reference
  
  more details: http://bio-bwa.sourceforge.net/bwa.shtml


  #### 1.2- Sorting SAM into coordinate order and save as BAM:
  Different approaches can be used: samtools or Picard. 
  ```
  samtools sort -o sorted.bam input.sam
  picard sortsam 
  ```
  more details: http://www.htslib.org/doc/samtools-sort.html#DESCRIPTION //
  https://broadinstitute.github.io/picard/command-line-overview.html#SortSam

  This creates a sorted BAM file, with the same content as the SAM file. A BAM file is the binary version of a SAM, BAM or CRAM file
  ```
  samtools view 
  ```  

  more details: http://www.htslib.org/doc/samtools-view.html


  #### 1.3- Mapping statistics: 
   a) Depth of coverage (number of times a nucleotide is read during sequencing)
   ```
   samtools depth
   gatk DepthOfCoverage
   ```   
   
   b) Breath of coverage (percentage of bases of my reference genome that are covered at a certain depth)
   ```
    bedtools genomecov -ibam sorted_alig_reads.bam -g reference.fa.fai -max 20 > sample_metrics.dcov 2> >(tee $logfile)
   ```
   Bedtools will compute a histogram of coverage for the genome. 
   -max flag to specify maximum depth of the histogram
   -d per-base
   The default ouput format is formed by 5 columns: chromosome, depth of coverage from features in input file, number of bases on chromosome (with depth equal to 2, size of chromosome in base pairs and fraction of bases on chromosome with depth equal to column 2.   
   more details: https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html
   
   c) Callable loci (total number of sites sequenced with an specified min number of reads)
   ```
   gatk CallableLoci
   ```
   d) Stats and flagstat (stats collects statistics from BAM files and outputs in a text format,
    flagstat counts 13 different  categories)
  ```
  samtools stats sorted_alig_reads.bam > sample_metrics.stats
  samtools flagstat sorted_alig_reads.bam > sample_metrics.flagstat
  ```

  #### 1.4- Local realignment around indels (?)
 
  #### 1.5- Base quality score recalibration (?)
 
 
### 2- Variant calling
## A) OPTION 
  #### 2.1- Call variants: 

HaplotypeCaller calls for two types of variation (via local de-novo assembly): 
- single nucleotide polymorphisms (SNPs)
- insertion-deletinos (indels)

Variant calling > specifying parameters

-genotyping_mode:
	> DISCOVERY: will choose the most likely alleles out of the data
	> GENOTYPE_GIVEN_ALLELES: only use the alleles passed in from a VCF file 
		(for this, specify the file using -alleles argument)

-output_mode:
	> EMIT_VARIANTS_ONLY: only calls in sites that appear to be variant (high or low confidence)
	> EMIT_ALL_CONFIDENT_SITES: sites with high confidence, whether the site is variant or non variant
	> EMIT_ALL_SITES: calls at any callable site (high or low confidence). Not okay for DISCOVERY mode
	
- stand_emit_conf:minimun confidence threshold 

-  stand_call_conf:
 ```
GenomeAnalysisTK.jar -T HaplotypeCaller -R reference.fa -I reduced_reads.bam -L 20 -- (specify parameters) -O raw_variants.vcf
 ```
This creates a file containing all the sites that the HaplotypeCaller evaluated to be potentially variant
It is a VCF file and contains SNPs and indels that should be filtered before they can be used

Consolidate GVCFs 
 ```
 Genomics DBImport 
 ```
 
Joint-cal cohort 
```
GenotypeGVCFs 
```
#### 2.2- option A) Variant quality score recalibration 
This option could not work with our samples. This is thought to be used with human analyses
This step requires high-quality sets of known variants to use as training and truth resources, which for many organisms are not yet available. It also requires quite a lot of data in order to learn the profiles of good vs bad variants
Better to use hard-filtering instead

2 steps:
- machine learning for well-calibrated probability to each variant call
- uses this score to filter raw calls


 #### 2.2- option B) Variant filtration 



  ### 3- Addition of variant annotations 
  In case we forget to specify an annotation, or you realize only later that a certain annotation would be useful, we can add variants later on. 


## B) OPTION
Using BCFtools from SAMtools. Can be single- and multi-sample calling 

1- Generate VCF containing genotype likelihoods for one or multiple alignement files: 
```
bcftools mpileup -C 50 -q -Q -Ou -f reference.fa.fai -r mapped_reads.bam mapped_reads2.bam mapped_reads3.bam
```
Flags that might be interesting: 
-C coefficient for downgrading mapping quality for reads containig excessive missmatches. The recommended value for BWA is 50
-q minimum mapping quality 
-Q base quality
-Ou output uncompressed vcf
-r specify region. Requires the alignement files to be indexed. 

2- Variant calling command from the output of mpileup command:
```
bcftools call -Ou -m -v
```
Flags that might be interesting:
-m multi allelic caller
-v variants only

3- Index:
```
tabix -p vcf 
```

4- Stats:
```
bcftools stats -F ref.fa 
```

3- Apply fixed-threshold filters:
```
bcftools filter 
```

4- zip:



 ### Done! 
We can start playing with the data! 
Good luck!






