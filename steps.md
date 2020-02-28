## Bioinformatic steps for variant calling from FASTQ

### Softwares: 
- BWA
- SAMtools
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
  #### 0.1- Generate the fasta gile index: 
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

  #### 0.4- Create a text file with chromosome names 



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

  This creates a file called *aligned_reads.sam* containing the aligned reads from all input files, combined, annotated and aligned to   the reference
  SAM files cannot be opened with basic unix commands. We need to use samtools:
  ```
  samtools view 
  ```
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

  #### 1.4- Local realignment around indels
 
  #### 1.5- Base quality score recalibration
 
 
 
 
### 2- Variant calling
  #### 2.1- Call variants: 
  ```
  gatk --java-options "-Xmx4g" HaplotypeCaller -R reference.fa -I input.bam -O output.g.vcf.gz -ERC GVCF
  ```
  There are two types of variation: 
  - single-nucleotide polymorphisms (SNPs)
  - insertion-deletios (indels)
  HaplotypeCaller is capable of calling both simultaneously
  The end product of the command will be a VCF file containing raw calls that should be filtered before they can be used in downstream  analyses 
  ```
  gatk --java-options "-Xmx4g -Xmx4g" GenomicsDBImport 
  ```
 
 
  ```
  gatk --java-options "-Xmx4g -Xmx4g" GatherVcfs 
  ```
  
 
  #### 2.2- option A) Variant quality score recalibration
  Separe for SNPs and indels
  For SNPs: 
  For indels:
 
 #### 2.2- option B) Variant filtration 
 

  ### 3- Addition of variant annotations 
  In case we forget to specify an annotation, or you realize only later that a certain annotation would be useful, we can add variants later on. 




 ### Done! 
We can start playing with the data! 
Good luck!






