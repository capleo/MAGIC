# Detection of Copy Number Variations (CNVs) #####
# Author: Svenja Mager
# email: s.mager@santannapisa.it

# Note: "sample" in the commands stands for one of the following: A632, B73, F7, H99, HP301, Mo17, W153R
## create a binned (0.1MB bin size) bed file (chr, start, end) for each chromomsome
# example given for chromosome 1 (to be repeated for all chromosomes and samples):
library(edgeR)
library(GenomicRanges)
library(cn.mops)
# Set number formatting options to prevent scientific notation
options(scipen = 999, digits = 20)
# chromosome 1
chromosome <- "1" # chromosome name
start_position <- 0 # chromosome start-position
end_position <- 307041717  # chromosome end-position
# Calculate the number of bins needed
bin_size <- 100000  # 0.1 MB in base pairs
num_bins <- ceiling((end_position - start_position + 1) / bin_size)
# Create a data frame to store the regions
regions_df <- data.frame(
  chr = character(num_bins),
  start = integer(num_bins),
  end = integer(num_bins))
# Fill in the data frame with the bins
for (i in 1:num_bins) {
  bin_start <- start_position + (i - 1) * bin_size
  bin_end <- bin_start + bin_size - 1
  regions_df[i, "chr"] <- chromosome
  regions_df[i, "start"] <- bin_start
  regions_df[i, "end"] <- bin_end}
# Write the data frame to the regions.bed file
write.table(regions_df, "regions_chr1.bed", sep = "\t", col.names = FALSE,
            row.names = FALSE, quote = FALSE)

# #####################################################################################################################################################
# #  need to use bedtools
# ## create count data for each chromosome: getting read counts for each bin with the help of the before-produced single-chromosome bam files and the binned bed files
# # example given for chromosome 1 (to be repeated for all chromosomes):
# bedtools multicov -bams A632_chr1.sorted.bam B73_chr1.sorted.bam F7_chr1.sorted.bam H99_chr1.sorted.bam HP301_chr1.sorted.bam Mo17_chr1.sorted.bam W153R_chr1.sorted.bam -bed regions_chr1.bed > counts_chr1.txt
# #####################################################################################################################################################

### detect CNVs based on read counts

# example given for chromosome 1 (to be repeated for all chromosomes):
# in R load the following packages
library(edgeR)
library(GenomicRanges)
library(cn.mops)
# prepare input count matrix for cn.mops
# Read the count data from the before-produced file
counts <- read.table("counts_chr1.txt", header = FALSE)
# Separate the genomic coordinates from the count data
genomic_coordinates <- counts[, 1:3]
sample_counts <- counts[, -(1:3)]
# Provide the sample names in the same order as in the count matrix
sample_names <- c("A632", "B73", "F7", "H99", "HP301", "Mo17", "W153R")
colnames(sample_counts) <- sample_names
# Combine the genomic coordinates into a single string column and set it as row names
genomic_coordinates_str <- paste(genomic_coordinates$V1, genomic_coordinates$V2, genomic_coordinates$V3, sep = ":")
rownames(sample_counts) <- genomic_coordinates_str
# Create the count_matrix with sample counts
count_matrix <- sample_counts
# find cnvs with cn.mops
inmatrix<-as.matrix(count_matrix)
res <- cn.mops(inmatrix)
res <- calcIntegerCopyNumbers(res)
cnvs <- data.frame(res@cnvs)
cnvr <- data.frame(res@cnvr)
# save tables with CNVs and CNV regions
write.table(cnvs, file = "cnvs_chr1.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
write.table(cnvr, file = "cnvr_chr1.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
