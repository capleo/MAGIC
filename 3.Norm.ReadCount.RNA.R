# The script does normalize Illumina RNA-seq read counts
# Author: Ettore Riccucci

#on B73 V4
setwd("/RNA/founders_onB73andMo17/mapping_STAR/leaf_on_B73")
library(edgeR)
library(tidyverse)
#Read in the counts table
rawdata <- read.delim("gene_count_leafonB73.txt", row.names = 1)
head(rawdata)
dim(rawdata)
genotype<-c(read.delim('samples.txt')$Genotype)
dim(genotype)
y<-DGEList(counts=rawdata[,1:54], genes=rownames(rawdata), group=genotype)

#We filter out lowly expressed genes
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
dim(y)
##We now perform normalization steps
#The calcNormFactors function normalizes the library sizes by finding a set of scaling factors
#for the library sizes that minimizes the log-fold changes between the samples for most genes.
#The default method for computing these scale factors uses a trimmed mean of M-values
#(TMM) between each pair of samples. We call the product of the original library size
#and the scaling factor the effective library size. The effective library size replaces the original
#library size in all downsteam analyses.
y <- calcNormFactors(y)
#Now we can see the scaling factors: these should be "reasonably" similar among all samples
#A normalization factor below one indicates that a small number of high count genes
#are monopolizing the sequencing, causing the counts for other genes to be lower than would
#be usual given the library size. As a result, the library size will be scaled down, analogous
#to scaling the counts upwards in that library. Conversely, a factor above one scales up the
#library size, analogous to downscaling the counts.
y$samples

#you can get normalized counts like this:
normcounts<-cpm(y, normalized.lib.sizes = T)
dim(normcounts)
write.table(normcounts, file = "normcounts_leafonB73.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
write.table(rownames(normcounts), file='normcountgene_lonB73.txt')
#You can also get it in log scale like this:
#normcounts<-cpm(y, normalized.lib.sizes = T, log = T)
