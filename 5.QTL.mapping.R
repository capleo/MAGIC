# ## QTL mapping

# Preliminary
rm(list=ls()) # clean workspace
options(stringsAsFactors = F)
wd <- "/home/l.caproni/Documents/MAGIC_CAM/output"
setwd(wd) # set wd

#load libraries
library(qtl2)
library(qtl2ggplot)
library(qtl2convert)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(patchwork)
library(ggplot2)

# Define color palette for founders
pal <- c("#F0A000", "#808080", "#F08080", "#00A0F0",
         "#00A000", "#F00000", "#9000E0", "#8E2323", "blue")

save(pal, file="../input/pal.founders.Rdata")

# load data
MAGIC<-read_cross2("MAGIC.CAM.yaml")
MAGIC <- drop_nullmarkers(MAGIC)  #drops non informative markers
summary(MAGIC)

# create directories to store results
subfolder_names <- c("F.contrib","CAM","geno.prob.RILs") 
 for (j in 1:length(subfolder_names)){
   folder<-dir.create(paste0(wd,"/",subfolder_names[j]))
 }

## calculate genotype probabilities with Hidden Markov Model HMM
pr <- calc_genoprob(MAGIC, 
                    map = MAGIC$gmap, # or map_est
                    error_prob = 0.01, #function allows genotyping error
                    quiet = T,
                    cores = 8
                    )

save(pr, file="genotype.probs.CAM.Rdata")
# load(file="genotype.probs.CAM.Rdata")


#Perform a GENOME SCAN with a linear model
kinship_loco <- calc_kinship(pr, "loco", cores = 16)
out <- scan1(genoprobs = pr, pheno = MAGIC$pheno, kinship = kinship_loco)

save(kinship_loco, out, file="genome.scan.CAM.Rdata")
# load(file="genome.scan.9.CAM.Rdata")

# Perform PERMUTATION TEST
operm<-scan1perm(genoprobs = pr, pheno = MAGIC$pheno, n_perm = 999)

thr1 <- summary_scan1perm(operm, alpha = 0.01)
thr5 = summary(operm)

save(operm, thr1, thr5, thr10, file ="permutation.test.CAM.Rdata")
# load(file="permutation.test.CAM.Rdata")

#Find LOD peaks in a set of LOD curves
peaks<-find_peaks(out, MAGIC$pmap, 
                  threshold = thr5, #manual or from summary permutations
                  expand2markers = T, #QTL intervals expanded to markers
                  prob = 0.99
                  #peakdrop = 1  #the amount to drop to get QTL intervals 
                  ) # #(can be vector for each LOD score column)


save(peaks, file="peaks_CAM.Rdata")
load(file = "peaks_CAM.Rdata")
write.table(peaks, file = "CAM/peaks.cam.joint.txt", quote = FALSE, sep="\t", row.names = FALSE)

# LOD score plots
setwd("CAM/")

for (p in 1:ncol(out)){ # 
  
  phenotype <- colnames(out)[p]
  # plot scan
  pdf(file=paste(phenotype, "LOD.pdf", sep="."), height = 7, width = 14)
  plot_scan1(out, map = MAGIC$pmap, lodcolumn = phenotype, 
             main= phenotype, ylim = c(0,(max(out, lodcolumn = p))*1.2),
             col = "black", lwd =0.8, cex =3)
  abline(h=thr1[,phenotype], col="red", lty=2, cex = 0.8)
  abline(h=thr5[,phenotype], col="orange", lty=2, cex = 0.8)
  dev.off()
  
} 

setwd("..")

