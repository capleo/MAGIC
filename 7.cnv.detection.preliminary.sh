##### Detection of CNVs #####
# Author: Svenja Mager

# Note: "sample" in the commands stands for one of the following: A632, B73, F7, H99, HP301, Mo17, W153R

### prepare files for CNV detection based on read counts

## extract single chromosomes from bam files
# example given for chromosome 1 (to be repeated for all chromosomes and samples):
samtools view sample.bam -b -o sample_chr1.bam 1:0-307041717 # the numbers indicate the chromosome and the start- and end-position of the chromosome

## sort and index the produced bam files
# example given for chromosome 1 (to be repeated for all chromosomes and samples):
samtools sort sample_chr1.bam -o sample_chr1.sorted.bam
samtools index sample_chr1.sorted.bam
