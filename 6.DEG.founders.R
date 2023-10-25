# Differential Expression analysis
# This script does identify Differentially Expressed Genes (DEGs) based on founders' haplotype effects
# and it is based on the use if RNA seq data on the founder lines
# Author: Leonardo Caproni 
# email: l.caproni@santannapisa.it

# Preliminary
rm(list=ls()) # clean workspace
options(stringsAsFactors = F)
setwd("~/MAGIC/output/") # 

#load libraries
library(qtl2)
library(qtl2convert)
library(RColorBrewer)
library(fpc)
library(pvclust)
library(tidyverse)
library(edgeR)
library(qtl2ggplot)
library(patchwork)
library(ggplot2)

# load data
MAGIC<-read_cross2("MAGIC.CAM.yaml") #cross2 object
MAGIC <- drop_nullmarkers(MAGIC)  #drops non informative markers
summary(MAGIC)

load(file="genotype.probs.map.inter.Rdata") # geno probabilities
load(file="genome.scan.CAM.Rdata") # scan out
load(file="../input/pal.founders.Rdata") # palette
load(file="permutation.test.CAM.Rdata") # thrs based on permutations
load(file = "peaks_CAM.Rdata") # your peaks

# Set Parameters
prob = 0.99 # Maximum marginal probability
FDR = 0.10 # multiple testing FDR threshold for DEG
max_span = 30 #max qtl span in Mb
min_span = 0.015 #min qtl span

#Make dirs to store 
dir.create(paste0(wd, "/founder.effect.CAM"))
dir.create(paste0(wd, "/DEG.movaps.fourth.leaf.CAM"))
dir.create(paste0(wd, "/QTL.RILs.state.CAM"))

# create query genes function
query_genes <- create_gene_query_func("../input/gff/maize_genes.sqlite")

#load normalized counts on fourth leaf
counts.n.leaf <- read.delim(file="../input/normalizedcounts_edgeR_leafonB73.txt", row.names = 1)
#names(counts.n.leaf) # check names
names(counts.n.leaf) <- gsub(x = names(counts.n.leaf), pattern = "W153", replacement = "W153R") 

newcounts.leaf<-list() #get averages across gene counts
for (i in unique(founders(MAGIC))){
  if(length(grep(i, colnames(counts.n.leaf)))>1){
    newcounts.leaf[[i]]<-rowMeans(counts.n.leaf[,grep(i, colnames(counts.n.leaf))])
  } else {
    newcounts.leaf[[i]]<-counts.n.leaf[,grep(i, colnames(counts.n.leaf))]
  }
}

# create norm counts leaf
counts.leaf<-do.call(cbind, newcounts.leaf) #get a df

# subset peaks
# test only QTL in which the interval is less than max_span Mb
peaks$qtl_span <- peaks$ci_hi-peaks$ci_lo
peaks_selected <- peaks %>% filter(qtl_span < max_span) 
peaks_selected <- peaks_selected %>% filter(qtl_span > min_span) 

#QTL effects BLUPs
for (p in 1:nrow(peaks_selected)){ # p=1
  
  phenotype <- peaks_selected$lodcolumn[p]
  # find position of max LOD score selectively
  peak_Mbp <- peaks_selected$pos[p]
  peak_round <- round(peak_Mbp)
  chr_tmp <- as.numeric(peaks_selected$chr[p])
  tmp_c2eff <- scan1coef (pr[,chr_tmp], 
                          MAGIC$pheno[,phenotype], 
                          kinship = kinship_loco[[chr_tmp]],
                          #cores = 8
  )
  
  # set physical window
  start<- peaks_selected$ci_lo[p] #- (peaks_selected$ci_hi[p]-peaks_selected$ci_lo[p])*0.10
  stop<- peaks_selected$ci_hi[p] #+ (peaks_selected$ci_hi[p]-peaks_selected$ci_lo[p])*0.10
  
  interval_tmp <- c(start, stop) #set the physical interval
  m_keep <- find_marker(MAGIC$pmap, chr_tmp, interval = interval_tmp) # grab the makers in the interval
  
  # Plot phenotype vs genotype for a single putative QTL and a single phenotype
  g.prob.x <- maxmarg(pr, MAGIC$pmap, chr=chr_tmp,
                      pos=peak_Mbp, #postion (gen map) # check
                      minprob = prob,
                      return_char=TRUE)
  # plot such probabilities
  pdf(paste( "founder.effect.CAM/pheno.vs.geno.at.QTL",phenotype, "Chr",chr_tmp, peak_round, "Mbp", "pdf", sep="."))
  plot_pxg(g.prob.x, MAGIC$pheno[,phenotype], swap_axes = F,
           SEmult = 2, xlab=phenotype, main=paste("Chromosome", chr_tmp, peak_round, "Mb p =", prob, sep=" " ), 
           #ylab=founders(MAGIC),
           sort = F)
  dev.off()
  
  # cut the coefficient to the confidence interval
  tmp_c2eff.red<-tmp_c2eff[m_keep,]
  tmp_c2eff.red <- as.matrix(tmp_c2eff.red)#[,c(1:8,10)] #exclude CML91
  
  ##################################################################################
  # CLUSTERING - Fourth Leaf
  ##################################################################################
  tmpc2eff.l<-tmp_c2eff.red[,c(1,2,4:8)]
  
  #get k means clusters
  tmpkm.l <- kmeansruns(t(tmpc2eff.l), 
                        krange = 1:ncol(tmpc2eff.l)-1, 
                        runs = 1000, scaledata = T,
                        criterion = "asw")
  #fviz_cluster(tmpkm, t(tmpc2eff))
  tmpclst.l <- tmpkm.l$cluster
  
  # ascribe letter factor groups to the k-means clusters, ordered from A is low to C is high
  km.order.l <- as.numeric(names(sort(tmpkm.l$centers[,1])))
  names(km.order.l) <- toupper(letters)[1:length(tmpkm.l$size)]
  km.order.l <- sort(km.order.l)
  clus.order.l <- factor(names(km.order.l[tmpkm.l$cluster]))
  
  ##################################################################################
  
  # now let's grab probabilities at the QTL for all RILs
  m_pos <- find_markerpos(MAGIC, m_keep, na.rm = TRUE) # find physical position of the markers
  
  # prepare a data.frame to store the allelic state of each RIL at each marker within the QTL
  out_tmp <- data.frame(matrix(NA, nrow=n_ind_geno(MAGIC), ncol=0))
  for(m in 1:nrow(m_pos)){ # m=1
    g.prob.tmp <- maxmarg(pr, MAGIC$pmap, chr=chr_tmp,
                          pos=m_pos$pmap [m], #postion (gen map) # check
                          minprob = prob,
                          return_char=T) # if false returns the numerical ID of the founder
    out_tmp_qtl <- data.frame(g.prob.tmp)
    colnames(out_tmp_qtl)[1] <- rownames(m_pos)[m]
    out_tmp <- cbind(out_tmp, out_tmp_qtl)
    
  }
  # fix the allelic state of each RIL at he QTL
  write.csv(out_tmp, paste("QTL.RILs.state.CAM/QTL", phenotype, chr_tmp, round(peak_Mbp),"Mb","csv", sep = "."), quote = F)
  
  ## now identify Differentially Expressed Genes (DEGs) within each QTL ci ##
  genes <- query_genes(chr_tmp, start, stop)
  #keep only genes
  genes <- genes[genes$type == "gene",]
  
  targetgene.leaf <- counts.leaf[rownames(counts.leaf) %in% genes$gene_id,] # fourth leaf
  targetgene.leaf <- targetgene.leaf[,c(1,2,3:7)] # with CML91 targetgene.leaf[,c(1,2,8,3:7)]
  
  ##################################################################################
  # DEG - Fourth Leaf
  ##################################################################################
  
  names(tmpclst.l)<-founders(MAGIC)[c(1,2,4:8)] #FOURTH LEAF no CML91
  
  tmpclst.ok.l <- rbind(tmpclst.l, clus.order.l) #Let's put in order
  tmpclst.ok.l <- tmpclst.ok.l [2,]
  #tmpclst.ok <- tmpclst.ok[-9]# and adapt to the counts table (no CML91) 
  
  d <- DGEList(counts = targetgene.leaf, group= tmpclst.ok.l)
  d <- calcNormFactors(d)
  #plotMDS(d,col=as.numeric(d$samples$group))
  #get dispersion
  d1 <- estimateCommonDisp(d, verbose=T)
  d1 <- estimateTagwiseDisp(d1)
  #plotBCV(d1)
  
  et31.l <- exactTest(d1, pair=factor(c(min(tmpclst.ok.l),max(tmpclst.ok.l)))) # compare groups max and 1
  degout.l<-data.frame(topTags(et31.l, n=10))
  write.csv(et31.l$table, file = paste("DEG.movaps.fourth.leaf.CAM/", phenotype,".chr.", chr_tmp,".all.GM.csv", sep = ""))
  
  #get threshold
  DEG.l<-degout.l[which(degout.l$FDR<FDR),]
  if (nrow(DEG.l)>0) {
    write.csv(DEG.l, file = paste("DEG.movaps.fourth.leaf.CAM/", phenotype,".chr.", chr_tmp,".csv", sep = ""))
    for(g in 1:nrow(DEG.l)){ # g=2
      toplot<-counts.n.leaf[which(rownames(counts.n.leaf) == row.names(DEG.l)[g]),]
      #make a boxplot
      toplot<-data.frame(t(toplot))
      toplot$geno<-rownames(toplot)
      
      toplot$geno<-gsub("\\..*", "",toplot$geno)
      toplot$geno<-as.factor(toplot$geno)
      
      toplot <- toplot[!toplot$geno=="CML91",] #exclude CML91
      
      bplot <- ggplot(toplot, aes(y=toplot[,1], fill=geno))+geom_boxplot()+
        scale_fill_manual(values = pal[c(1,2,4:8)]) +
        guides(fill=guide_legend("Founder")) +
        ggtitle(row.names(DEG.l)[g]) +
        labs(y="raw counts") +
        theme_bw()
      
      coef <- ggplot_coef(tmp_c2eff.red[,c(1,2,4:8)], map=MAGIC$pmap, columns = 1:7) +
        geom_line(linewidth=0.2) +
        scale_colour_manual(values = pal[c(1,2,4:8)], labels = founders(MAGIC)[c(1,2,4:8)]) +
        guides(col=guide_legend("Founder")) + 
        theme_bw()
      
      plot <- coef | bplot
      
      ggsave(plot, file=paste("DEG.movaps.fourth.leaf.CAM/", phenotype,".chr.", chr_tmp,".", 
                              colnames(toplot)[1], ".jpg", sep = ""), height = 5, width = 12 )
      
    }
    
  }
}
